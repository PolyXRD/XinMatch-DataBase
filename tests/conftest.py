"""Shared fixtures for the XinMatch-DataBase tool tests.

The repository deliberately has no test CIF archive checked in: the real data
lives in ``D:/Project/XRD/sync-COD`` (``XRD_DATA_ROOT`` override).  The fixtures
below synthesize tiny CIFs and a tiny archive in COD fan-out layout so the
builder / validator / engine can be exercised without touching the 1.8 GiB
source tree.
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import pytest

#: A minimal but complete P1 structure: gemmi parses it, it has a real cell,
#: two element species and a nonzero isotropic displacement (so Debye-Waller is
#: actually applied).
CIF_TEMPLATE = """\
data_{name}
_chemical_formula_sum '{formula}'
_cell_length_a {a}
_cell_length_b {b}
_cell_length_c {c}
_cell_angle_alpha {alpha}
_cell_angle_beta {beta}
_cell_angle_gamma {gamma}
_symmetry_space_group_name_H-M 'P 1'
loop_
_space_group_symop_operation_xyz
  'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
_atom_site_U_iso_or_equiv
Si1 Si 0.0 0.0 0.0 1.0 0.01
O1 O 0.25 0.25 0.25 1.0 0.02
O2 O 0.75 0.75 0.75 1.0 0.02
"""

#: Cell that reproduces the 2theta=180 deg defect of the pre-2026-10 engine:
#: c = 7.703 A = 10 * lambda/2, and with a = b = 4.2606 the reciprocal-matrix
#: rounding puts the (0,0,10) reflection at d = 0.7703000000000001 > lambda/2,
#: so the old ``2d > lambda`` test kept it, its LP factor exploded and the
#: normalised pattern collapsed to a single bogus line.  The fixed engine must
#: return the full pattern (>= 100 reflections) instead.
BOUNDARY_A = 4.2606
BOUNDARY_B = 4.2606
BOUNDARY_C = 7.703

#: A file gemmi cannot map: no ``_cell_*`` tags at all.  The engine must report
#: it as uncomputable instead of emitting one bogus d = 1 A reflection.
NO_CELL_CIF = """\
data_nocell
_chemical_formula_sum 'La Mg Na O W'
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
La1 La 0.0 0.0 0.0
Mg1 Mg 0.25 0.25 0.25
"""


def render_cif(
    name: str = "1000001",
    *,
    formula: str = "O2 Si",
    a: float = 5.431,
    b: float = 5.431,
    c: float = 5.431,
    alpha: float = 90.0,
    beta: float = 90.0,
    gamma: float = 90.0,
) -> str:
    """Return a valid P1 CIF document as text."""
    return CIF_TEMPLATE.format(
        name=name, formula=formula, a=a, b=b, c=c,
        alpha=alpha, beta=beta, gamma=gamma,
    )


@pytest.fixture
def write_cif(tmp_path: Path) -> Callable[..., Path]:
    """Return ``write(name, **cell_kwargs) -> path`` writing a valid CIF."""

    def _write(name: str = "1000001", **kwargs: float) -> Path:
        path = tmp_path / f"{name}.cif"
        path.write_text(render_cif(name, **kwargs), encoding="utf-8")
        return path

    return _write


@pytest.fixture
def boundary_cif(tmp_path: Path) -> Path:
    """A CIF that reproduces the old 2theta=180 deg single-bogus-line defect."""
    path = tmp_path / "1099999.cif"
    path.write_text(
        render_cif("1099999", a=BOUNDARY_A, b=BOUNDARY_B, c=BOUNDARY_C),
        encoding="utf-8",
    )
    return path


@pytest.fixture
def no_cell_cif(tmp_path: Path) -> Path:
    """A CIF without any cell tags (uncomputable)."""
    path = tmp_path / "1099998.cif"
    path.write_text(NO_CELL_CIF, encoding="utf-8")
    return path


@pytest.fixture
def tiny_archive(tmp_path: Path) -> Path:
    """A 12-phase archive in COD fan-out layout.

    Layout: ``<root>/1/00/00/10000NN.cif`` — the same shape
    :func:`xinmatch.core.db.cif_scan.iter_cod_cifs` walks, so builder tests can
    run a real build end-to-end in about a second.
    """
    root = tmp_path / "archive"
    for i in range(12):
        cod_id = 1_000_000 + i
        d = root / str(cod_id)[:1] / str(cod_id)[1:3] / str(cod_id)[3:5]
        d.mkdir(parents=True, exist_ok=True)
        cell = 5.431 + 0.01 * i
        (d / f"{cod_id}.cif").write_text(
            render_cif(str(cod_id), a=cell, b=cell, c=cell), encoding="utf-8"
        )
    return root
