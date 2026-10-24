"""Regression tests for the synced CIF metadata parser (``cif_meta``).

The DB repository used to vendor an older ``cif_meta.py``.  The synced version
adds two things the builder depends on:

* :func:`cell_parameters_from_block` -- the shared "recover the cell from CIF
  text tags" helper that ``peak_compute`` uses to avoid computing a pattern from
  gemmi's placeholder cell;
* a validated cell tuple in :class:`PhaseMeta` (positive lengths, angles in
  ``(0, 180)``) instead of only a ``None`` check.

Everything else (formula / element resolution, space-group normalisation,
volume fallback) must keep behaving exactly as before, because those values are
already stored in released libraries.
"""

from __future__ import annotations

from pathlib import Path

import gemmi
import pytest

from xinmatch.core.db import cif_meta
from xinmatch.core.db.cif_meta import (
    ELEMENT_BIT_TABLE_HASH,
    ELEMENT_MASK_SCHEME,
    IUPAC_ELEMENTS,
    PhaseMeta,
    cell_parameters_from_block,
    element_mask,
    elements_from_formula,
    elements_from_sites,
    parse_meta,
)
from xinmatch.core.db.cif_scan import DEFAULT_INORG_COD_ROOT

#: Real inorganic COD archive (read-only), skipped when absent.
SRC_COD = DEFAULT_INORG_COD_ROOT

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

#: No formula at all: the element list must fall back to the atom-site list.
NO_FORMULA_CIF = """\
data_no_formula
_cell_length_a 5.0
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
Fe1 Fe3+ 0.0 0.0 0.0
O1  O2-  0.25 0.25 0.25
"""

#: ``_chemical_formula_sum`` holds a base64 payload (real COD defect); the
#: parser must reject it and fall back to the structural formula / sites.
BASE64_FORMULA_CIF = """\
data_base64_formula
_chemical_formula_sum 'iVBORw0KGgoAAAANSUhEUg'
_chemical_formula_structural 'Fe2 O3'
_cell_length_a 5.0
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
Fe1 Fe 0.0 0.0 0.0
O1  O  0.25 0.25 0.25
"""

#: Invalid cell: gemmi accepts ``a = -5`` into a UnitCell, but the text-tag cell
#: parser must reject it (positive lengths are required).
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


# ---------------------------------------------------------------------------
# 1. parse_meta on a well-formed CIF
# ---------------------------------------------------------------------------
def test_parse_meta_reads_scalar_fields(write_cif) -> None:
    path = write_cif("1000001")
    meta = parse_meta(path, entry_id="96-100-0001")
    assert meta is not None
    assert isinstance(meta, PhaseMeta)
    assert meta.entry_id == "96-100-0001"
    assert meta.formula == "O2 Si"
    assert meta.elements == ["O", "Si"]
    assert meta.space_group == "P 1"
    assert meta.crystal_system == "triclinic"
    assert meta.cell == (5.431, 5.431, 5.431, 90.0, 90.0, 90.0)
    assert meta.volume == pytest.approx(5.431**3)
    assert meta.density is None
    assert meta.mineral_name is None
    assert meta.element_mask == element_mask(["O", "Si"])


def test_parse_meta_entry_id_defaults_to_empty(write_cif) -> None:
    meta = parse_meta(write_cif("1000001"))
    assert meta is not None
    assert meta.entry_id == ""


def test_parse_meta_reads_bibliographic_and_mineral_fields(tmp_path: Path) -> None:
    text = """\
data_biblio
_chemical_formula_sum 'Ca F2'
_chemical_name_mineral 'fluorite'
_cell_length_a 5.463
_cell_length_b 5.463
_cell_length_c 5.463
_cell_angle_alpha 90.0
_cell_angle_beta 90.0
_cell_angle_gamma 90.0
_symmetry_space_group_name_H-M 'F m -3 m'
_cell_density 3.18
_journal_year 1914
_publ_author_name 'Bragg W. H.'
_journal_name_full 'Proc. R. Soc. London'
_cod_database_code 1000001
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Ca1 Ca 0.0 0.0 0.0
F1  F  0.25 0.25 0.25
"""
    meta = parse_meta(_write(tmp_path, "biblio.cif", text), entry_id="96-100-0001")
    assert meta is not None
    assert meta.formula == "Ca F2"
    assert meta.mineral_name == "fluorite"
    assert meta.space_group == "F m -3 m"
    assert meta.crystal_system == "cubic"
    assert meta.density == pytest.approx(3.18)
    assert meta.year == 1914
    assert meta.author == "Bragg W. H."
    assert meta.journal == "Proc. R. Soc. London"
    assert meta.quality == "1000001"


# ---------------------------------------------------------------------------
# 2. parse_meta edge cases (missing cell / formula / unparseable file)
# ---------------------------------------------------------------------------
def test_parse_meta_keeps_a_phase_without_cell(no_cell_cif: Path) -> None:
    """No cell = no peaks, but the metadata itself is still readable."""
    meta = parse_meta(no_cell_cif, entry_id="96-109-9998")
    assert meta is not None
    assert meta.cell is None
    assert meta.volume is None
    assert meta.elements == ["La", "Mg", "Na", "O", "W"]


def test_parse_meta_recovers_cell_from_text_tags(tmp_path: Path) -> None:
    """gemmi cannot map the cell, the text tags can."""
    path = _write(tmp_path, "recoverable.cif", RECOVERABLE_CIF)
    meta = parse_meta(path)
    assert meta is not None
    assert meta.cell == (4.0, 4.0, 5.0, 90.0, 90.0, 90.0)
    assert meta.volume == pytest.approx(80.0)


def test_parse_meta_rejects_negative_cell_parameters(tmp_path: Path) -> None:
    meta = parse_meta(_write(tmp_path, "negative.cif", NEGATIVE_CELL_CIF))
    assert meta is not None
    assert meta.cell is None


def test_parse_meta_returns_none_for_unreadable_files(tmp_path: Path) -> None:
    assert parse_meta(tmp_path / "does_not_exist.cif") is None
    assert parse_meta(_write(tmp_path, "junk.cif", "not a CIF\n")) is None
    assert parse_meta(_write(tmp_path, "empty.cif", "")) is None


def test_parse_meta_falls_back_to_structural_formula(tmp_path: Path) -> None:
    """A base64 ``_chemical_formula_sum`` must not hide the chemistry."""
    meta = parse_meta(_write(tmp_path, "b64.cif", BASE64_FORMULA_CIF))
    assert meta is not None
    assert elements_from_formula(meta.formula) == ["Fe", "O"]
    assert meta.elements == ["Fe", "O"]


def test_parse_meta_falls_back_to_atom_sites(tmp_path: Path) -> None:
    meta = parse_meta(_write(tmp_path, "nosites.cif", NO_FORMULA_CIF))
    assert meta is not None
    assert meta.formula is None
    assert meta.elements == ["Fe", "O"]


# ---------------------------------------------------------------------------
# 3. Formula / site element extraction
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("O2 Si", ["O", "Si"]),
        ("C3 O6 Sr", ["C", "O", "Sr"]),
        ("Ca F2", ["Ca", "F"]),
        ("D2 O", ["H", "O"]),  # deuterium is hydrogen for X-ray scattering
        ("T2 O", ["H", "O"]),
        ("", []),
        (None, []),
        ("?", []),
        ("iVBORw0KGgoAAAANSUhEUg", []),  # base64 payload
    ],
)
def test_elements_from_formula(formula: str | None, expected: list[str]) -> None:
    assert elements_from_formula(formula) == expected


def test_elements_from_sites_reads_type_symbol_heads() -> None:
    class _Site:
        def __init__(self, type_symbol: str) -> None:
            self.type_symbol = type_symbol

    class _Structure:
        sites = [_Site("Fe3+"), _Site("O2-"), _Site("D1"), _Site("?")]

    assert elements_from_sites(_Structure()) == ["Fe", "H", "O"]


def test_elements_from_sites_tolerates_bad_input() -> None:
    assert elements_from_sites(object()) == []


# ---------------------------------------------------------------------------
# 4. Frozen element-mask contract
# ---------------------------------------------------------------------------
def test_element_bit_table_is_frozen() -> None:
    assert len(IUPAC_ELEMENTS) == 118
    assert len(set(IUPAC_ELEMENTS)) == 118
    assert ELEMENT_MASK_SCHEME == "iupac128-v1"
    assert len(ELEMENT_BIT_TABLE_HASH) == 16


def test_element_mask_low_and_high_halves() -> None:
    low, high = element_mask(["H", "O"])
    assert low == (1 << 0) | (1 << 7)
    assert high == 0
    # U has atomic number 92 -> zero-based bit 91 -> bit 27 of the high half
    low, high = element_mask(["U"])
    assert low == 0
    assert high == 1 << (91 - 64)
    # unknown symbols contribute no bit
    assert element_mask(["Xx"]) == (0, 0)


# ---------------------------------------------------------------------------
# 5. Shared cell-from-text helper
# ---------------------------------------------------------------------------
def test_cell_parameters_from_block_strips_uncertainty() -> None:
    block = gemmi.cif.read_string(
        "data_x\n"
        "_cell_length_a 4.913(2)\n_cell_length_b 5.404(2)\n_cell_length_c 5.404(2)\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 120\n"
    ).sole_block()
    assert cell_parameters_from_block(block) == (4.913, 5.404, 5.404, 90.0, 90.0, 120.0)


def test_cell_parameters_from_block_rejects_incomplete_or_absent_cell() -> None:
    partial = gemmi.cif.read_string(
        "data_x\n_cell_length_b 4\n_cell_length_c 5\n"
        "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
    ).sole_block()
    assert cell_parameters_from_block(partial) is None
    assert cell_parameters_from_block(None) is None


def test_cell_parameters_from_block_is_exported() -> None:
    """``peak_compute`` imports this helper, so ``__all__`` must advertise it."""
    assert "cell_parameters_from_block" in cif_meta.__all__


# ---------------------------------------------------------------------------
# 6. Real COD source CIFs (skipped when the archive is not on this machine)
# ---------------------------------------------------------------------------
def _cif_path(native_id: int) -> Path:
    text = str(native_id).zfill(7)
    return SRC_COD / text[0] / text[1:3] / text[3:5] / f"{text}.cif"


@pytest.mark.parametrize("native_id", [1008057, 1545882, 2104629])
def test_parse_meta_on_real_cod_cifs(native_id: int) -> None:
    path = _cif_path(native_id)
    if not path.is_file():
        pytest.skip(f"COD source archive not on this machine: {path}")
    meta = parse_meta(path, entry_id="96-100-0001")
    assert meta is not None
    assert meta.entry_id == "96-100-0001"
    assert meta.elements, "real COD phase must resolve at least one element"
    if native_id != 2104629:  # 96-210-4629 carries no cell in the source data
        assert meta.cell is not None
        assert meta.volume is not None
