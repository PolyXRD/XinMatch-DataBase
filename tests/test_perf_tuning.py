"""Regression tests for the adaptive structure-factor block size.

The block size only decides which reflections are processed together; each
reflection's sum is computed in the same order either way, so choosing it from
the atom count must not change a single bit of ``|F|**2``.  These tests pin that
down and make the tuning explicit rather than a hidden magic number.
"""
from __future__ import annotations

import gemmi
import numpy as np

from xinmatch.core.db import peak_compute as pc


def _structure_args(cif_path):
    ss = gemmi.read_small_structure(str(cif_path))
    expanded = pc._expand_sites(ss)
    recip = pc._reciprocal_matrix(ss.cell)
    assert expanded is not None and recip is not None
    enum = pc._enumerate_hkl(
        ss.cell,
        dmin=pc.DEFAULT_DMIN,
        dmax=pc.DEFAULT_DMAX,
        max_hkl=pc.DEFAULT_MAX_HKL,
    )
    assert enum is not None
    hkl, _d = enum
    frac, occ, b_iso, symbols = expanded
    return hkl, frac, occ, b_iso, symbols, recip


def test_adaptive_block_matches_historical_fixed_block(write_cif):
    cif = write_cif("1000001")
    args = _structure_args(cif)
    reference = pc._structure_factors(*args, chunk=8192)
    for chunk in (256, 512, 1024, 2048, 4096):
        got = pc._structure_factors(*args, chunk=chunk)
        assert np.array_equal(reference, got), f"chunk={chunk} is not bit-identical"


def test_default_block_is_bit_identical_to_fixed_8192(write_cif):
    cif = write_cif("1000001")
    args = _structure_args(cif)
    reference = pc._structure_factors(*args, chunk=8192)
    assert np.array_equal(reference, pc._structure_factors(*args))


def test_block_size_is_capped_for_large_cells():
    # 720 atoms -> 262144 // 720 = 364 rows, within [256, 8192]; 10 atoms -> capped
    # at the historical 8192 rows.
    def block(n_atoms: int) -> int:
        return min(
            pc._SF_BLOCK_MAX,
            max(pc._SF_BLOCK_MIN, pc._SF_BLOCK_ELEMENTS // max(n_atoms, 1)),
        )

    assert block(10) == 8192
    assert block(720) == 364
    assert block(100_000) == pc._SF_BLOCK_MIN
