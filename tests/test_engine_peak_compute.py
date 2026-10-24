"""Regression tests for the synced peak engine (``peak_compute``).

The DB repository used to vendor an older copy of ``peak_compute.py`` that
carried two real defects (both discovered in the 2026-10 acceptance run):

1. **2 theta = 180 deg boundary fake line.**  When a cell axis is an integer
   multiple of ``wavelength / 2`` (e.g. ``c = 7.703 = 10 * 0.7703``) the
   back-scattered reflection is kept by a naive ``2 * d > wavelength`` test
   through floating-point rounding; its Lorentz-polarisation factor diverges and
   the normalised pattern collapses to a single bogus line at d = 0.7703 A.
2. **gemmi placeholder-cell fake line.**  When gemmi cannot map the CIF cell it
   leaves ``(1, 1, 1, 90, 90, 90)`` in place, which yields one bogus reflection
   at d = 1.0 A.  Such a phase must be reported as uncomputable (or recovered
   from the CIF text tags) instead of written to the library.

Both fixtures come from ``tests/conftest.py`` so builder, validator and engine
tests exercise the same archive layout.
"""

from __future__ import annotations

from pathlib import Path

import gemmi
import numpy as np
import pytest

from tests.conftest import NO_CELL_CIF
from xinmatch.core.db import peak_compute as pc
from xinmatch.core.db.cif_meta import cell_parameters_from_block
from xinmatch.core.db.cif_scan import DEFAULT_INORG_COD_ROOT

WL = pc.DEFAULT_WAVELENGTH
D_BOUNDARY = WL / 2.0

#: Real inorganic COD archive (read-only).  Absent on a clean checkout, in
#: which case the real-file cases below are skipped.
SRC_COD = DEFAULT_INORG_COD_ROOT

#: gemmi cannot map this CIF (``_cell_length_a`` hides in a trailing
#: single-column ``loop_``), but all six parameters are in the text tags.
RECOVERABLE_CIF = """\
data_recoverable_cell
_space_group_name_H-M   'P 1'
_cell_length_b          4.0000
_cell_length_c          5.0000
_cell_angle_alpha       90.0
_cell_angle_beta        90.0
_cell_angle_gamma       90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Si1 Si 0.0 0.0 0.0 1.0
O1  O  0.25 0.25 0.25 1.0
loop_
_cell_length_a
4.0000
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    """Write *text* to ``tmp_path/name`` and return the path."""
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def _ss(text: str) -> gemmi.SmallStructure:
    """Build a small structure straight from CIF text."""
    return gemmi.make_small_structure_from_block(gemmi.cif.read_string(text).sole_block())


# ---------------------------------------------------------------------------
# 1. 2 theta = 180 deg boundary reflection (LP divergence / fake line)
# ---------------------------------------------------------------------------
def test_boundary_fixture_really_contains_180deg_reflection(boundary_cif: Path) -> None:
    """Prove the fixture hits the old defect: (0 0 +-10) passes the naive test."""
    ss = gemmi.read_small_structure(str(boundary_cif))
    enum = pc._enumerate_hkl(
        ss.cell,
        dmin=pc.DEFAULT_DMIN,
        dmax=pc.DEFAULT_DMAX,
        max_hkl=pc.DEFAULT_MAX_HKL,
    )
    assert enum is not None
    hkl, d = enum
    naive = (2.0 * d) > WL
    strict = (2.0 * d) > WL * (1.0 + pc._MEASURABLE_REL_EPS)
    leaked = np.nonzero(naive & ~strict)[0]
    assert leaked.size >= 1, "fixture no longer reproduces the boundary defect"
    leaked_hkl = {tuple(int(v) for v in row) for row in hkl[leaked]}
    assert (0, 0, 10) in leaked_hkl or (0, 0, -10) in leaked_hkl


def test_boundary_cif_returns_full_pattern(boundary_cif: Path) -> None:
    """The fixed engine returns the full pattern, not the single bogus line."""
    result = pc.compute_pattern(boundary_cif, max_peaks=0)
    assert result is not None
    d, i = result
    assert d.size > 100, "old engine collapsed this fixture to one fake peak"
    assert bool((2.0 * d > WL * (1.0 + pc._MEASURABLE_REL_EPS)).all())
    assert not np.any(np.isclose(d, D_BOUNDARY, rtol=0.0, atol=1e-9))
    assert float(d.min()) >= pc.DEFAULT_DMIN
    assert float(d.max()) <= pc.DEFAULT_DMAX + 1e-9
    assert float(i.max()) == pytest.approx(pc.INTENSITY_SCALE)


def test_lorentz_polarization_finite_for_any_angle() -> None:
    angles = np.array([0.0, 1e-6, 30.0, 90.0, 179.9, 179.9999, 180.0])
    lp = pc._lorentz_polarization(angles)
    assert bool(np.isfinite(lp).all())
    assert bool((lp > 0.0).all())


def test_lorentz_polarization_normal_angles_bit_identical() -> None:
    """The defensive floor must not change any normal-angle value."""
    normal = np.array([10.0, 30.0, 60.0, 90.0, 120.0, 150.0])
    theta = np.radians(normal) / 2.0
    expected = (1.0 + np.cos(2.0 * theta) ** 2) / (
        np.sin(theta) * np.sin(theta) * np.cos(theta)
    )
    np.testing.assert_array_equal(pc._lorentz_polarization(normal), expected)


# ---------------------------------------------------------------------------
# 2. gemmi placeholder cell (fake d = 1 A line)
# ---------------------------------------------------------------------------
def test_cell_is_usable_rejects_placeholder_and_degenerate_cells() -> None:
    class _FakeCell:
        """Duck-typed cell with an impossible angle (gemmi refuses to build it)."""

        a = 5.0
        b = 5.0
        c = 5.0
        alpha = 0.0
        beta = 90.0
        gamma = 90.0

    assert not pc._cell_is_usable(None)
    assert not pc._cell_is_usable(gemmi.UnitCell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0))
    assert not pc._cell_is_usable(gemmi.UnitCell(0.0, 5.0, 5.0, 90.0, 90.0, 90.0))
    assert not pc._cell_is_usable(_FakeCell())
    assert pc._cell_is_usable(gemmi.UnitCell(5.431, 5.431, 5.431, 90.0, 90.0, 90.0))


def test_gemmi_leaves_placeholder_for_recoverable_cif() -> None:
    """Fixture precondition: gemmi really cannot map this cell."""
    assert not pc._cell_is_usable(_ss(RECOVERABLE_CIF).cell)


def test_cell_recovered_from_text_tags(tmp_path: Path) -> None:
    path = _write(tmp_path, "recoverable.cif", RECOVERABLE_CIF)
    for result in (pc.compute_pattern_from_text(RECOVERABLE_CIF), pc.compute_pattern(path)):
        assert result is not None, "a text-recoverable cell must not be uncomputable"
        d, _i = result
        assert d.size > 1
        assert float(d.max()) == pytest.approx(5.0)  # c axis
        assert pc.DEFAULT_DMIN <= float(d.min()) <= float(d.max())


def test_no_cell_cif_is_uncomputable(no_cell_cif: Path) -> None:
    assert pc.compute_pattern(no_cell_cif) is None
    assert pc.compute_pattern_from_text(NO_CELL_CIF) is None
    # even when handed the placeholder-cell structure directly
    assert pc.compute_pattern_from_structure(_ss(NO_CELL_CIF)) is None


def test_cell_parameters_from_block_parses_uncertainty(boundary_cif: Path) -> None:
    text = boundary_cif.read_text(encoding="utf-8")
    text = text.replace("_cell_length_a 4.2606", "_cell_length_a 4.2606(10)")
    block = gemmi.cif.read_string(text).sole_block()
    assert cell_parameters_from_block(block) == (4.2606, 4.2606, 7.703, 90.0, 90.0, 90.0)
    partial = gemmi.cif.read_string(
        "data_x\n_cell_length_b 4\n_cell_length_c 5\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
    ).sole_block()
    assert cell_parameters_from_block(partial) is None
    assert cell_parameters_from_block(None) is None


# ---------------------------------------------------------------------------
# 3. Pattern semantics (truncation / normalisation / dmin)
# ---------------------------------------------------------------------------
def test_pattern_is_sorted_by_descending_d(write_cif) -> None:
    path = write_cif("1000001")
    result = pc.compute_pattern(path, max_peaks=0)
    assert result is not None
    d, i = result
    assert bool((np.diff(d) <= 0.0).all()), "peaks must be ordered by descending d"
    assert bool((i > 0.0).all())
    assert d.dtype == np.float64
    assert i.dtype == np.float64


def test_max_peaks_keeps_the_strongest_reflections(write_cif) -> None:
    path = write_cif("1000001")
    full = pc.compute_pattern(path, max_peaks=0)
    assert full is not None
    d_full, i_full = full
    assert d_full.size > 5, "fixture must have more than 5 reflections"
    truncated = pc.compute_pattern(path, max_peaks=5)
    assert truncated is not None
    d_trunc, i_trunc = truncated
    assert d_trunc.size == 5
    # every truncated line belongs to the untruncated pattern, bit for bit
    for d_value, i_value in zip(d_trunc, i_trunc, strict=True):
        match = np.nonzero(d_full == d_value)[0]
        assert match.size == 1
        assert i_full[match[0]] == i_value
    # and the retained set is exactly the 5 strongest
    strongest = np.argsort(-i_full, kind="stable")[:5]
    assert set(np.round(d_trunc, 12)) == set(np.round(d_full[strongest], 12))


def test_normalize_scales_strongest_peak(write_cif) -> None:
    path = write_cif("1000001")
    raw = pc.compute_pattern(path, max_peaks=0, normalize=False)
    normalised = pc.compute_pattern(path, max_peaks=0, normalize=True)
    assert raw is not None and normalised is not None
    _d_raw, i_raw = raw
    _d_norm, i_norm = normalised
    assert float(i_raw.max()) != pytest.approx(pc.INTENSITY_SCALE)
    assert float(i_norm.max()) == pytest.approx(pc.INTENSITY_SCALE)
    np.testing.assert_allclose(
        i_norm, i_raw * (pc.INTENSITY_SCALE / float(i_raw.max())), rtol=1e-12
    )


def test_custom_dmin_bounds_the_pattern(write_cif) -> None:
    path = write_cif("1000001")
    default = pc.compute_pattern(path, max_peaks=0)
    assert default is not None
    high_dmin = pc.compute_pattern(path, max_peaks=0, dmin=1.5)
    assert high_dmin is not None
    d_high, _i_high = high_dmin
    assert float(d_high.min()) >= 1.5
    assert d_high.size <= default[0].size


def test_new_path_reader_and_text_reader_agree(write_cif) -> None:
    """``read_small_structure`` and ``make_small_structure_from_block`` agree.

    ``build_pipeline`` relies on this: it parses the CIF once into a block and
    builds the structure from that block, while the old two-step path used
    ``gemmi.read_small_structure``.  The patterns must be bit-identical.
    """
    path = write_cif("1000001")
    text = path.read_text(encoding="utf-8")
    from_path = pc.compute_pattern(path, max_peaks=0)
    from_text = pc.compute_pattern_from_text(text, max_peaks=0)
    assert from_path is not None and from_text is not None
    np.testing.assert_array_equal(from_path[0], from_text[0])
    np.testing.assert_array_equal(from_path[1], from_text[1])


def test_unparseable_cif_is_none(tmp_path: Path) -> None:
    junk = _write(tmp_path, "junk.cif", "this file is not a CIF at all\n")
    assert pc.compute_pattern(junk) is None
    assert pc.compute_pattern_from_text("this file is not a CIF at all\n") is None


# ---------------------------------------------------------------------------
# 4. Real COD source CIFs (skipped when the archive is not on this machine)
# ---------------------------------------------------------------------------
def _cif_path(native_id: int) -> Path:
    text = str(native_id).zfill(7)
    return SRC_COD / text[0] / text[1:3] / text[3:5] / f"{text}.cif"


def _require_cif(native_id: int) -> Path:
    path = _cif_path(native_id)
    if not path.is_file():
        pytest.skip(f"COD source archive not on this machine: {path}")
    return path


@pytest.mark.parametrize(
    ("entry_id", "native_id", "expected_peaks"),
    [
        ("96-100-8057", 1008057, 122),
        ("96-151-1248", 1511248, 228),
        ("96-153-3251", 1533251, 152),
        ("96-154-5882", 1545882, 413),
        ("96-210-1175", 2101175, 79),
    ],
)
def test_real_180deg_phases_restore_full_pattern(
    entry_id: str, native_id: int, expected_peaks: int
) -> None:
    """These phases collapsed to one fake d = 0.7703 line in the old engine."""
    path = _require_cif(native_id)
    result = pc.compute_pattern(path, max_peaks=0)
    assert result is not None, entry_id
    d, i = result
    assert d.size == expected_peaks, entry_id
    assert bool((2.0 * d > WL * (1.0 + pc._MEASURABLE_REL_EPS)).all())
    assert not np.any(np.isclose(d, D_BOUNDARY, rtol=0.0, atol=1e-9))
    assert float(i.max()) == pytest.approx(pc.INTENSITY_SCALE)


@pytest.mark.parametrize(
    ("entry_id", "native_id", "expected_peaks", "expected_dmax"),
    [
        ("96-156-3735", 1563735, 150, 19.739),
        ("96-156-3736", 1563736, 114, 10.1704),
        ("96-210-4629", 2104629, None, None),  # no cell in the source data
    ],
)
def test_real_placeholder_cell_phases(
    entry_id: str,
    native_id: int,
    expected_peaks: int | None,
    expected_dmax: float | None,
) -> None:
    path = _require_cif(native_id)
    result = pc.compute_pattern(path, max_peaks=pc.DEFAULT_MAX_PEAKS)
    if expected_peaks is None:
        assert result is None, f"{entry_id}: a phase without a cell must be uncomputable"
        return
    assert result is not None, entry_id
    d, _i = result
    assert d.size == expected_peaks, entry_id
    assert float(d.max()) == pytest.approx(expected_dmax, abs=1e-3), entry_id
    assert not np.any(np.isclose(d, 1.0, rtol=0.0, atol=1e-9)), entry_id
