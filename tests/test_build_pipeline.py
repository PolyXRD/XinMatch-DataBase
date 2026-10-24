"""Parity tests for :func:`xinmatch.core.db.build_pipeline.compute_phase`.

The builder used to parse every CIF twice: once in ``parse_meta`` and once in
``compute_pattern``.  ``compute_phase`` produces both results from a single
parse and must stay a drop-in replacement for the old two-call sequence.  The
tests below define that sequence independently (``old_two_step``) and compare:

* the returned ``reason`` string, character for character;
* every :class:`PhaseMeta` field (via ``dataclasses.asdict``);
* the ``(d, i)`` arrays **bit for bit** (``ndarray.tobytes``), not just
  approximately.

Real COD CIFs are exercised when ``D:/Project/XRD/sync-COD`` is present; on a
machine without the archive those cases are skipped, never failed.
"""

from __future__ import annotations

from dataclasses import asdict
from pathlib import Path

import gemmi
import numpy as np
import pytest

from tests.conftest import render_cif
from xinmatch.core.db.build_pipeline import compute_phase
from xinmatch.core.db.cif_meta import parse_meta, parse_meta_from_structure
from xinmatch.core.db.cif_scan import (
    DEFAULT_INORG_COD_ROOT,
    entry_id_for_cod,
)
from xinmatch.core.db.peak_compute import (
    DEFAULT_DMIN,
    DEFAULT_MAX_PEAKS,
    compute_pattern,
    compute_pattern_from_block,
)

#: Real inorganic COD archive (read-only), skipped when absent.
SRC_COD = DEFAULT_INORG_COD_ROOT

#: ``reason`` values of the public contract.
REASONS = {"", "no_cif", "no_meta", "no_peaks"}

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

NO_SITES_CIF = """\
data_no_sites
_chemical_formula_sum 'Si O2'
_cell_length_a 5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
"""

NEGATIVE_CELL_CIF = """\
data_negative_cell
_chemical_formula_sum 'Si O2'
_cell_length_a -5.0
_cell_length_b 5.0
_cell_length_c 5.0
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Si1 Si 0.0 0.0 0.0
O1  O  0.25 0.25 0.25
"""


def _write(tmp_path: Path, name: str, text: str) -> Path:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return path


def old_two_step(
    cif_path: Path,
    *,
    entry_id: str,
    dmin: float,
    max_peaks: int,
    normalize: bool = True,
) -> tuple[object, object, str]:
    """Reference implementation: the old two-call builder sequence.

    This is intentionally written out rather than imported, so the parity test
    cannot silently drift with the implementation under test.
    """
    path = Path(cif_path)
    if not path.is_file():
        return None, None, "no_cif"
    meta = parse_meta(path, entry_id=entry_id)
    if meta is None:
        return None, None, "no_meta"
    peaks = compute_pattern(path, dmin=dmin, max_peaks=max_peaks, normalize=normalize)
    if peaks is None:
        return None, None, "no_peaks"
    return meta, (peaks[0], peaks[1]), ""


def assert_same_phase(result, reference) -> None:
    """Assert that ``compute_phase`` matched the old two-step result exactly."""
    meta, peaks, reason = result
    old_meta, old_peaks, old_reason = reference
    assert reason == old_reason
    assert reason in REASONS
    if old_meta is None:
        assert meta is None
    else:
        assert meta is not None
        assert asdict(meta) == asdict(old_meta)
    if old_peaks is None:
        assert peaks is None
    else:
        assert peaks is not None
        assert len(peaks) == len(old_peaks) == 2
        for new_arr, old_arr in zip(peaks, old_peaks, strict=True):
            assert isinstance(new_arr, np.ndarray)
            assert new_arr.dtype == old_arr.dtype == np.float64
            assert new_arr.shape == old_arr.shape
            assert new_arr.tobytes() == old_arr.tobytes(), "peak arrays must be bit-identical"


# ---------------------------------------------------------------------------
# 1. Fixture files covering success and every failure reason
# ---------------------------------------------------------------------------
@pytest.fixture
def parity_cases(tmp_path: Path, write_cif, boundary_cif: Path, no_cell_cif: Path) -> dict:
    return {
        "normal": write_cif("1000001"),
        "boundary": boundary_cif,
        "no_cell": no_cell_cif,
        "recoverable": _write(tmp_path, "recoverable.cif", RECOVERABLE_CIF),
        "no_sites": _write(tmp_path, "no_sites.cif", NO_SITES_CIF),
        "negative_cell": _write(tmp_path, "negative_cell.cif", NEGATIVE_CELL_CIF),
        "broken": _write(tmp_path, "broken.cif", "this file is not a CIF\n"),
        "missing": tmp_path / "missing.cif",
    }


#: Expected ``reason`` per fixture, documenting the old builder's check order.
EXPECTED_REASONS = {
    "normal": "",
    "boundary": "",
    "no_cell": "no_peaks",
    "recoverable": "",
    "no_sites": "no_peaks",
    "negative_cell": "no_peaks",
    "broken": "no_meta",
    "missing": "no_cif",
}


@pytest.mark.parametrize("case", sorted(EXPECTED_REASONS))
def test_parity_with_old_two_step(parity_cases: dict, case: str) -> None:
    path = parity_cases[case]
    kwargs = {
        "entry_id": "96-100-0001",
        "dmin": DEFAULT_DMIN,
        "max_peaks": DEFAULT_MAX_PEAKS,
    }
    result = compute_phase(path, **kwargs)
    reference = old_two_step(path, **kwargs)
    assert result[2] == EXPECTED_REASONS[case]
    assert_same_phase(result, reference)


@pytest.mark.parametrize(
    ("dmin", "max_peaks", "normalize"),
    [
        (DEFAULT_DMIN, DEFAULT_MAX_PEAKS, True),
        (DEFAULT_DMIN, 7, True),
        (DEFAULT_DMIN, 0, True),  # 0 disables the truncation
        (1.2, 5, True),
        (2.0, 20, False),
        (0.5, 10, True),  # explicit dmin below the back-scattering limit
    ],
)
@pytest.mark.parametrize("case", ["normal", "boundary"])
def test_parity_for_builder_parameters(
    parity_cases: dict, case: str, dmin: float, max_peaks: int, normalize: bool
) -> None:
    path = parity_cases[case]
    kwargs = {
        "entry_id": "96-100-0001",
        "dmin": dmin,
        "max_peaks": max_peaks,
        "normalize": normalize,
    }
    assert_same_phase(compute_phase(path, **kwargs), old_two_step(path, **kwargs))


def test_peak_arrays_are_shared_by_identity_within_result(parity_cases: dict) -> None:
    """The arrays returned by ``compute_phase`` are the computed float64 arrays."""
    meta, peaks, reason = compute_phase(
        parity_cases["normal"],
        entry_id="96-100-0001",
        dmin=DEFAULT_DMIN,
        max_peaks=DEFAULT_MAX_PEAKS,
    )
    assert meta is not None and peaks is not None and reason == ""
    d, i = peaks
    assert d.dtype == np.float64 and i.dtype == np.float64
    assert d.size == i.size > 0
    assert bool((np.diff(d) <= 0.0).all())


def test_success_returns_meta_and_peaks_together(parity_cases: dict) -> None:
    meta, peaks, reason = compute_phase(
        parity_cases["boundary"],
        entry_id="96-109-9999",
        dmin=DEFAULT_DMIN,
        max_peaks=0,
    )
    assert reason == ""
    assert meta is not None
    assert meta.entry_id == "96-109-9999"
    assert peaks is not None
    d, _i = peaks
    assert d.size > 100, "boundary fixture must keep its full fixed pattern"


def test_failure_returns_both_none(parity_cases: dict) -> None:
    for case in ("missing", "broken", "no_cell", "no_sites", "negative_cell"):
        meta, peaks, reason = compute_phase(
            parity_cases[case],
            entry_id="96-100-0001",
            dmin=DEFAULT_DMIN,
            max_peaks=DEFAULT_MAX_PEAKS,
        )
        assert (meta, peaks) == (None, None)
        assert reason == EXPECTED_REASONS[case]


# ---------------------------------------------------------------------------
# 2. Entry-id handling (COD 7-digit / TCOD 8-digit shifted ids)
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("native_id", "offset"),
    [
        (1_000_001, 0),  # COD
        (1_000_001, 20_000_000),  # PCOD
        (11_234_567, 100_000_000),  # TCOD (8-digit native id)
    ],
)
def test_entry_id_is_stored_verbatim(
    parity_cases: dict, native_id: int, offset: int
) -> None:
    entry_id = entry_id_for_cod(native_id + offset)
    kwargs = {"entry_id": entry_id, "dmin": DEFAULT_DMIN, "max_peaks": 10}
    result = compute_phase(parity_cases["normal"], **kwargs)
    assert result[0] is not None
    assert result[0].entry_id == entry_id
    assert_same_phase(result, old_two_step(parity_cases["normal"], **kwargs))


# ---------------------------------------------------------------------------
# 3. Whole tiny archive (every builder slice must stay in parity)
# ---------------------------------------------------------------------------
def test_parity_over_tiny_archive(tiny_archive: Path) -> None:
    from xinmatch.core.db.cif_scan import iter_cod_cifs

    paths = list(iter_cod_cifs(tiny_archive))
    assert len(paths) == 12
    for index, (cod_id, cif_path) in enumerate(paths):
        kwargs = {
            "entry_id": entry_id_for_cod(cod_id),
            "dmin": DEFAULT_DMIN,
            "max_peaks": DEFAULT_MAX_PEAKS,
        }
        assert_same_phase(
            compute_phase(cif_path, **kwargs),
            old_two_step(cif_path, **kwargs),
        ), f"parity broke on archive entry #{index} ({cif_path})"


# ---------------------------------------------------------------------------
# 4. Single-parse guards (the optimization itself)
# ---------------------------------------------------------------------------
def test_compute_phase_parses_the_cif_exactly_once(
    parity_cases: dict, monkeypatch: pytest.MonkeyPatch
) -> None:
    """One CIF parse + one structure build, no ``read_small_structure`` rerun.

    This is the throughput optimization under test: the old two-step path called
    ``gemmi.read_small_structure`` twice (plus 1-2 ``gemmi.cif.read`` calls).
    """
    calls = {"cif_read": 0, "read_small_structure": 0, "build_structure": 0}
    real_cif_read = gemmi.cif.read
    real_read_small_structure = gemmi.read_small_structure
    real_build_structure = gemmi.make_small_structure_from_block

    def counting_cif_read(*args: object, **kwargs: object):
        calls["cif_read"] += 1
        return real_cif_read(*args, **kwargs)

    def counting_read_small_structure(*args: object, **kwargs: object):
        calls["read_small_structure"] += 1
        return real_read_small_structure(*args, **kwargs)

    def counting_build_structure(*args: object, **kwargs: object):
        calls["build_structure"] += 1
        return real_build_structure(*args, **kwargs)

    monkeypatch.setattr(gemmi.cif, "read", counting_cif_read)
    monkeypatch.setattr(gemmi, "read_small_structure", counting_read_small_structure)
    monkeypatch.setattr(
        gemmi, "make_small_structure_from_block", counting_build_structure
    )

    meta, peaks, reason = compute_phase(
        parity_cases["normal"],
        entry_id="96-100-0001",
        dmin=DEFAULT_DMIN,
        max_peaks=DEFAULT_MAX_PEAKS,
    )
    assert reason == ""
    assert meta is not None and peaks is not None
    assert calls == {"cif_read": 1, "read_small_structure": 0, "build_structure": 1}


def test_block_entry_point_matches_path_entry_point(parity_cases: dict) -> None:
    """``compute_pattern_from_block`` is the block-based twin of ``compute_pattern``."""
    for case in ("normal", "boundary", "recoverable", "no_cell", "negative_cell"):
        path = parity_cases[case]
        block = gemmi.cif.read(str(path)).sole_block()
        from_block = compute_pattern_from_block(block, max_peaks=DEFAULT_MAX_PEAKS)
        from_path = compute_pattern(path, max_peaks=DEFAULT_MAX_PEAKS)
        if from_path is None:
            assert from_block is None, case
            continue
        assert from_block is not None, case
        for new_arr, old_arr in zip(from_block, from_path, strict=True):
            assert new_arr.tobytes() == old_arr.tobytes(), case


@pytest.mark.parametrize("case", ["normal", "boundary", "recoverable", "no_cell"])
def test_meta_helper_matches_public_parse_meta(parity_cases: dict, case: str) -> None:
    """The shared metadata core returns exactly what the public reader returns."""
    path = parity_cases[case]
    block = gemmi.cif.read(str(path)).sole_block()
    structure = gemmi.make_small_structure_from_block(block)
    from_helper = parse_meta_from_structure(structure, block, entry_id="96-100-0001")
    public = parse_meta(path, entry_id="96-100-0001")
    if public is None:
        assert from_helper is None
    else:
        assert from_helper is not None
        assert asdict(from_helper) == asdict(public)


# ---------------------------------------------------------------------------
# 5. Real COD CIFs (skipped when the archive is not on this machine)
# ---------------------------------------------------------------------------
#: Boundary-defect, placeholder-cell, no-cell and ordinary phases.
REAL_COD_IDS = [1008057, 1563735, 1563736, 2104629, 1511248, 1000001]


def _real_cif(native_id: int) -> Path:
    text = str(native_id).zfill(7)
    return SRC_COD / text[0] / text[1:3] / text[3:5] / f"{text}.cif"


@pytest.mark.parametrize("native_id", REAL_COD_IDS)
def test_parity_on_real_cod_cifs(native_id: int) -> None:
    path = _real_cif(native_id)
    if not path.is_file():
        pytest.skip(f"COD source archive not on this machine: {path}")
    for kwargs in (
        {"entry_id": entry_id_for_cod(native_id), "dmin": DEFAULT_DMIN, "max_peaks": DEFAULT_MAX_PEAKS},
        {"entry_id": entry_id_for_cod(native_id), "dmin": 1.2, "max_peaks": 0},
    ):
        result = compute_phase(path, **kwargs)
        reference = old_two_step(path, **kwargs)
        assert result[2] == reference[2], f"{path.name}: reason differs from two-step"
        assert_same_phase(result, reference)


def test_real_cod_phase_metadata_is_populated() -> None:
    """A real phase keeps every metadata field the old two-step produced."""
    path = _real_cif(1008057)
    if not path.is_file():
        pytest.skip(f"COD source archive not on this machine: {path}")
    meta, peaks, reason = compute_phase(
        path,
        entry_id=entry_id_for_cod(1008057),
        dmin=DEFAULT_DMIN,
        max_peaks=DEFAULT_MAX_PEAKS,
    )
    assert reason == ""
    assert meta is not None and peaks is not None
    assert meta.elements
    assert meta.cell is not None
    assert meta.volume is not None


def test_render_cif_helper_is_the_expected_fixture_format() -> None:
    """Guard against a conftest fixture change silently weakening these tests."""
    text = render_cif("1000001")
    assert "_cell_length_a 5.431" in text
    assert "_atom_site_type_symbol" in text
