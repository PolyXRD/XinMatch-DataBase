"""CIF -> d-I powder pattern computation (gemmi parser + vectorised physics).

This module is the **single source of truth** for turning a CIF file into a
simulated powder diffraction pattern.  It feeds the v3.0.0 ``.xinDB`` rebuild,
which regenerates the reference peak table from the COD CIF archive instead of
trusting previously computed values.

Peak convention (frozen -- see :data:`PEAK_CONVENTION`)
------------------------------------------------------
The convention reproduces ``pymatgen.analysis.diffraction.xrd.XRDCalculator``
exactly (verified numerically in ``tests/db/test_peak_compute.py``).  Its key
property is that **symmetry-equivalent reflections are merged and their
intensities summed**, i.e. the multiplicity of every reflection family is
included.  That is what makes the result a *powder* pattern:

    fluorite CaF2 (F m -3 m) -> the strongest line is (220) d = 1.931 A, not
    (111) d = 3.154 A; (111) carries ~82 % of the (220) height.

An asymmetric-unit engine ("one reflection per family, no multiplicity") gets
this backwards and must not be used to generate reference data.

Pipeline
--------
1. Parse with :mod:`gemmi` and expand the asymmetric unit to every atom of the
   unit cell (:meth:`gemmi.SmallStructure.get_all_unit_cell_sites`).
2. Enumerate **every** ``(h, k, l)`` with ``dmin <= d <= dmax``.  No
   asymmetric-unit filter: systematic absences come out of the physics, not out
   of a table (an extinct reflection simply gets ``|F|**2 == 0``).
3. Structure factor per reflection::

       s      = |g| / 2   with g = h a* + k b* + l c*     (no 2*pi)
       f(s)   = Z - 41.78214 * s**2 * sum_i a_i exp(-b_i s**2)
                (Cromer-Mann, the parametrisation pymatgen's XRDCalculator uses)
       T(s)   = exp(-B s**2),  B = 8*pi**2 * U_iso        (Debye-Waller;
                T == 1 when the CIF carries no temperature factor)
       F(hkl) = sum_j occ_j T_j(s) f_j(s) exp(2*pi*i*(h x_j + k y_j + l z_j))
       I(hkl) = |F(hkl)|**2 * LP(2 theta)
       LP     = (1 + cos**2 2 theta) / (sin**2 theta * cos theta)   (unpolarised)
4. Merge reflections whose 2-theta differ by less than :data:`TWO_THETA_TOL`
   degrees, **summing** their intensities -- this is where multiplicity enters.
5. Drop merged groups with ``I <= 0`` (systematically absent / non-finite).
6. Normalise so the strongest group equals :data:`INTENSITY_SCALE` (1000).
7. Truncate to the ``max_peaks`` **strongest** groups, then sort by descending
   d-spacing (= ascending 2-theta).

Deliberate differences from the legacy helper
---------------------------------------------
The legacy :func:`tools.build_cod_full_xinmatch.compute_peaks_from_cif` used a
gemmi asymmetric-unit enumeration, which (a) excluded multiplicity, (b) dropped
systematic absences inconsistently, and (c) produced all-NaN peak lists for CIFs
with a non-finite structure factor.  All three defects are gone.

Two behaviours are intentional and worth stating explicitly:

* **Unmeasurable reflections are dropped** rather than clamped.  ``2 d < lambda``
  means ``2 theta >= 180 deg``; clamping the ``asin`` argument there makes the
  Lorentz-polarisation factor diverge to ~``1e19``.
* **The Debye-Waller term is honoured** when the CIF provides
  ``_atom_site_B_iso_or_equiv`` / ``_atom_site_U_iso_or_equiv``.  Ignoring it
  inflates the high-angle intensities by up to an order of magnitude.

Why not call ``PolyXRD-DataBase/tools/fastxrd.py`` directly?
------------------------------------------------------------
``fastxrd`` is bit-equivalent to pymatgen and was used as the reference while
fixing this convention, but it needs ``fastcif`` for CIF parsing, and ``fastcif``
rejects a non-trivial share of real COD files (``_space_group_symop_*`` blocks,
quoted values, ...).  ``gemmi`` parses those same files correctly, so the physics
is re-implemented here on top of gemmi's parser.  On files both engines can read,
the numerical agreement with ``fastxrd``/pymatgen is exact.

The module is free of any database / Qt dependency: it only needs ``gemmi`` and
``numpy``, which keeps it usable from CLI and batch tooling.
"""

from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
from typing import TypeAlias

import gemmi
import numpy as np

from .cif_meta import cell_parameters_from_block

__all__ = [
    "DEFAULT_WAVELENGTH",
    "DEFAULT_DMIN",
    "DEFAULT_DMAX",
    "DEFAULT_MAX_PEAKS",
    "DEFAULT_MAX_HKL",
    "INTENSITY_SCALE",
    "TWO_THETA_TOL",
    "PEAK_CONVENTION",
    "PatternResult",
    "compute_pattern",
    "compute_pattern_from_block",
    "compute_pattern_from_structure",
    "compute_pattern_from_text",
]

#: Cu K-alpha-1 wavelength (A).
DEFAULT_WAVELENGTH = 1.5406

#: Lowest resolvable d-spacing (A).  ``lambda / 2`` = 0.7703 A is the physical
#: back-scattering limit (2 theta = 180 deg); nothing below it is measurable.
DEFAULT_DMIN = DEFAULT_WAVELENGTH / 2.0

#: Relative tolerance for the ``2 theta = 180 deg`` boundary.  A reflection whose
#: d is only "greater" than ``wavelength / 2`` because of floating-point
#: rounding (the 5 COD phases whose cell axis is an integer multiple of
#: 0.7703 A) is at exactly 180 deg: it is not measurable and the
#: Lorentz-polarisation factor diverges there.  Requiring
#: ``2 * d > wavelength * (1 + eps)`` drops those boundary reflections while
#: leaving every physically measurable reflection untouched.
_MEASURABLE_REL_EPS = 1e-9

#: Smallest ``sin(theta)`` / ``cos(theta)`` accepted by the Lorentz-polarisation
#: factor.  The classical LP expression diverges at both ends of the measurable
#: range (``1 / sin^2 theta`` as theta -> 0, ``1 / cos theta`` at 2 theta = 180
#: deg).  Flooring the two trig factors at the value for ``2 theta = 0.001`` /
#: ``179.999`` deg keeps the factor finite for *any* input -- including a caller
#: that passes a custom ``dmin`` below the back-scattering limit -- without
#: changing any value inside the normal peak-table range (2.94..180 deg).
_MIN_LP_TRIG = math.cos(math.radians(179.999) / 2.0)

#: Largest d-spacing kept (A).  Equals 2 theta = 2.94 deg for Cu K-alpha and
#: bounds the reflection count of very large cells.
DEFAULT_DMAX = 30.0

#: Number of strongest reflections kept per entry.  150 keeps every line that
#: can matter for identification (COD phases have 253 reflections on average and
#: 1187 at the 90th percentile) while bounding the library size: the intensity
#: truncation cuts the tail, which is exactly the part no matching algorithm
#: consumes.
DEFAULT_MAX_PEAKS = 150

#: Safety valve for pathological cells: a CIF that needs more reflections than
#: this is reported as uncomputable instead of exhausting memory.
DEFAULT_MAX_HKL = 400_000

#: Intensity normalisation target: the strongest reflection of every entry is
#: scaled to this value.  A single shared scale keeps cross-library comparisons
#: (and FoM thresholds) consistent.
INTENSITY_SCALE = 1000.0

#: Two reflections closer than this in 2 theta (degrees) belong to the same
#: merged group.  Symmetry-equivalent reflections agree analytically, so the
#: tolerance only has to absorb floating-point noise.
TWO_THETA_TOL = 1e-5

#: Relative intensity cut-off, as a fraction of the strongest reflection.
#: Systematically absent reflections do not evaluate to exactly zero in floating
#: point (``|F|**2`` lands around ``1e-32``), so a hard ``> 0`` test would keep
#: thousands of numerical ghosts.  0.001 % of the maximum is the same cut-off the
#: reference ``fastxrd`` pipeline uses, and it is far below any line that could
#: matter for phase identification.
REL_INTENSITY_CUTOFF = 1e-5

#: Human-readable identifier of the peak convention, meant to be stored in the
#: database ``meta`` table so a reader can detect a stale or foreign peak table.
PEAK_CONVENTION = "CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000"

#: Result of a peak computation: ``(d, i)`` or, with ``return_hkl=True``,
#: ``(d, i, hkl)``.
PatternResult: TypeAlias = tuple[np.ndarray, np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray]

#: Cromer-Mann coefficients keyed by element symbol.  Covers 99 elements; the
#: 20 missing ones are Es..Og, which do not occur in COD.
_PARAMS_PATH = Path(__file__).with_name("data") / "atomic_scattering_params.json"

#: The ``41.78214`` prefactor belongs to the pymatgen parametrisation: it folds
#: the wavelength factor and the 4-pi convention of ``s`` into one constant.
_CM_PREFACTOR = 41.78214


@lru_cache(maxsize=1)
def _scattering_table() -> dict[str, tuple[tuple[float, float], ...]]:
    """Load the Cromer-Mann coefficient table once per process."""
    with _PARAMS_PATH.open(encoding="utf-8") as fh:
        raw = json.load(fh)
    return {el: tuple((float(a), float(b)) for a, b in pairs) for el, pairs in raw.items()}


@lru_cache(maxsize=1)
def _atomic_numbers() -> dict[str, float]:
    """Atomic number per element symbol, as required by the Cromer-Mann form."""
    return {el: float(gemmi.Element(el).atomic_number) for el in _scattering_table()}


#: gemmi's default (placeholder) unit cell.  ``read_small_structure`` leaves it
#: in place when it cannot map the cell tags of a CIF; computing a pattern from
#: it would emit a single bogus reflection at d = 1 A, so such a phase must be
#: either recovered from the text tags or reported as uncomputable.
_DEFAULT_CELL_PARAMS = (1.0, 1.0, 1.0, 90.0, 90.0, 90.0)


def _cell_is_usable(cell: gemmi.UnitCell | None) -> bool:
    """True when *cell* is a real unit cell rather than gemmi's placeholder.

    ``gemmi.UnitCell.is_crystal()`` cannot be used as the "was the cell read?"
    test: it is ``False`` for any cell whose ``a`` is not greater than 1, which
    includes valid low-``a`` grids.  The placeholder is detected by its exact
    parameter tuple instead; a real crystal cell of (1, 1, 1, 90, 90, 90) does
    not exist (that is a 1 A^3 cell).
    """
    if cell is None:
        return False
    params = (cell.a, cell.b, cell.c, cell.alpha, cell.beta, cell.gamma)
    a, b, c, alpha, beta, gamma = params
    if not all(math.isfinite(value) for value in params):
        return False
    if params == _DEFAULT_CELL_PARAMS:
        return False
    if a <= 0.0 or b <= 0.0 or c <= 0.0:
        return False
    return 0.0 < alpha < 180.0 and 0.0 < beta < 180.0 and 0.0 < gamma < 180.0


def _read_cif_block(cif_path: str | Path) -> gemmi.cif.Block | None:
    """Return the sole block of a CIF file, or ``None`` when unparseable."""
    try:
        return gemmi.cif.read(str(cif_path)).sole_block()
    except Exception:
        return None


def _recover_cell(block: gemmi.cif.Block | None) -> gemmi.UnitCell | None:
    """Rebuild a unit cell from the CIF text tags (see ``cell_parameters_from_block``)."""
    params = cell_parameters_from_block(block)
    if params is None:
        return None
    try:
        cell = gemmi.UnitCell(*params)
    except Exception:  # pragma: no cover - gemmi validates its arguments
        return None
    return cell if _cell_is_usable(cell) else None


def _reciprocal_matrix(cell: gemmi.UnitCell) -> np.ndarray | None:
    """Rows are the reciprocal basis vectors ``a*, b*, c*`` (no 2*pi).

    ``|h a* + k b* + l c*| = 1 / d(hkl)``, so ``s = |g| / 2`` is directly
    ``sin(theta) / lambda``.  Returns ``None`` for a degenerate cell.
    """
    try:
        direct = np.array(
            [
                [getattr(cell.orthogonalize(gemmi.Fractional(*(1 if i == j else 0 for j in range(3)))), ax)
                 for ax in ("x", "y", "z")]
                for i in range(3)
            ],
            dtype=np.float64,
        )
    except Exception:
        return None
    det = float(np.linalg.det(direct))
    if not math.isfinite(det) or abs(det) < 1e-12:
        return None
    # The reciprocal basis is the inverse transpose of the direct basis.
    return np.linalg.inv(direct).T


def _expand_sites(
    ss: gemmi.SmallStructure,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[str, ...]] | None:
    """Expand the asymmetric unit to every atom of the unit cell.

    Returns ``(frac, occ, b_iso, symbols)`` or ``None`` when the structure
    carries no usable site.

    ``gemmi`` reports the isotropic displacement as ``U_iso`` in A^2, hence
    ``B_iso = 8 * pi**2 * U_iso``.  An occupancy of zero is treated as one: COD
    uses ``0`` for "not refined" rather than for "absent", and treating it as
    absent would silently delete atoms from the structure factor.
    """
    if not ss.sites:
        return None
    try:
        sites = ss.get_all_unit_cell_sites()
    except Exception:
        return None
    if not sites:
        return None

    frac = np.empty((len(sites), 3), dtype=np.float64)
    occ = np.ones(len(sites), dtype=np.float64)
    b_iso = np.zeros(len(sites), dtype=np.float64)
    symbols: list[str] = []
    for i, site in enumerate(sites):
        frac[i, 0] = site.fract.x
        frac[i, 1] = site.fract.y
        frac[i, 2] = site.fract.z
        o = float(site.occ)
        occ[i] = o if 0.0 < o <= 1.0 else 1.0
        u = float(site.u_iso)
        b_iso[i] = 8.0 * math.pi * math.pi * u if math.isfinite(u) and u > 0.0 else 0.0
        symbols.append(site.element.name)
    if not np.all(np.isfinite(frac)):
        return None
    return frac, occ, b_iso, tuple(symbols)


def _enumerate_hkl(
    cell: gemmi.UnitCell,
    *,
    dmin: float,
    dmax: float,
    max_hkl: int,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Enumerate every ``(h, k, l)`` with ``dmin <= d <= dmax``.

    A whole ``(k, l)`` plane is evaluated per ``h`` with
    :meth:`gemmi.UnitCell.calculate_d_array`, which is far cheaper than a Python
    loop over reflections.  Returns ``(hkl, d)``, or ``None`` when the cell needs
    more reflections than ``max_hkl`` or nothing falls in range.
    """
    try:
        hmax, kmax, lmax = cell.get_hkl_limits(dmin)
    except Exception:
        return None
    if hmax <= 0 or kmax <= 0 or lmax <= 0:
        return None

    k_arr = np.arange(-kmax, kmax + 1, dtype=np.int32)
    l_arr = np.arange(-lmax, lmax + 1, dtype=np.int32)
    plane = np.empty((k_arr.size * l_arr.size, 3), dtype=np.int32)
    plane[:, 1] = np.tile(k_arr, l_arr.size)
    plane[:, 2] = np.repeat(l_arr, k_arr.size)

    chunks: list[np.ndarray] = []
    d_chunks: list[np.ndarray] = []
    total = 0
    for h in range(-hmax, hmax + 1):
        plane[:, 0] = h
        try:
            d_plane = cell.calculate_d_array(plane)
        except Exception:
            continue
        # d is +inf for (0,0,0) and undefined for a degenerate cell.
        keep = np.isfinite(d_plane) & (d_plane >= dmin) & (d_plane <= dmax)
        if not keep.any():
            continue
        idx = np.nonzero(keep)[0]
        total += idx.size
        if total > max_hkl:
            return None
        chunks.append(plane[idx])
        d_chunks.append(d_plane[idx])
    if not chunks:
        return None
    return np.concatenate(chunks), np.concatenate(d_chunks)


#: Target number of ``(reflection x atom)`` phase values per structure-factor
#: block (~0.25 M doubles = 2 MiB of float64), i.e. the block size that keeps the
#: ``(m, n_atoms)`` temporaries (``cos``/``sin``/``w``) inside L2/L3 instead of
#: streaming through main memory.  This matters for the handful of large-cell
#: COD phases that dominate a library build: at 720 atoms the old fixed 8192-row
#: block allocated ~47 MiB *per temporary*, and capping by atom count measured
#: ~1.2x faster over a random real-CIF sample (~1.4x on the slowest phase).
#: The block boundary only decides which rows are processed together; every
#: row's sum is computed in the same order, so the result is bit-identical
#: (asserted by ``tests/test_perf_tuning.py``).
_SF_BLOCK_ELEMENTS = 262_144

#: Bounds of the adaptive block: never smaller than 256 rows (loop overhead) and
#: never larger than the historical 8192 rows (small cells are fine in one go).
_SF_BLOCK_MIN = 256
_SF_BLOCK_MAX = 8192


def _structure_factors(
    hkl: np.ndarray,
    frac: np.ndarray,
    occ: np.ndarray,
    b_iso: np.ndarray,
    symbols: tuple[str, ...],
    recip: np.ndarray,
    *,
    chunk: int | None = None,
) -> np.ndarray:
    """``|F(hkl)|**2`` for every reflection, vectorised over reflections.

    Reflections are processed in blocks so the ``(n_hkl, n_atoms)`` phase matrix
    never exceeds ``chunk * n_atoms`` elements.  ``chunk=None`` (the production
    default) sizes the block adaptively from the atom count (see
    :data:`_SF_BLOCK_ELEMENTS`); passing an explicit integer keeps the historical
    fixed-block behaviour (used by tests / benchmarks).
    """
    if chunk is None:
        chunk = min(
            _SF_BLOCK_MAX,
            max(_SF_BLOCK_MIN, _SF_BLOCK_ELEMENTS // max(len(symbols), 1)),
        )
    table = _scattering_table()
    znum = _atomic_numbers()

    unique = sorted(set(symbols))
    coeff = np.zeros((len(unique), 4, 2), dtype=np.float64)
    zarr = np.zeros(len(unique), dtype=np.float64)
    for i, el in enumerate(unique):
        # Deuterium and tritium are hydrogen isotopes: CIFs label them ``D`` and
        # ``T``, which are not element symbols, so map them onto H before the
        # lookup.  Without this they fall into the "unknown" branch below and
        # contribute no scattering power at all.
        lookup = {"D": "H", "T": "H"}.get(el, el)
        pairs = table.get(lookup)
        if pairs is None:
            # Only Es..Og land here, and they do not occur in COD.  Fall back to
            # the s -> 0 limit (f = Z) rather than silently dropping the atom.
            zarr[i] = float(gemmi.Element(el).atomic_number)
            continue
        for j, (a, b) in enumerate(pairs[:4]):
            coeff[i, j, 0] = a
            coeff[i, j, 1] = b
        zarr[i] = znum[lookup]

    el_of_atom = np.array([unique.index(s) for s in symbols], dtype=np.int64)
    hkl_f = hkl.astype(np.float64)
    s = np.linalg.norm(hkl_f @ recip, axis=1) / 2.0
    s2 = s * s

    # Occupancies and isotropic temperature factors almost always take a single
    # value across the cell (measured: 9 of the 10 largest-cell CIFs in COD have
    # exactly one unique value for each).  Computing ``exp(-B s^2)`` once *per
    # atom* is therefore mostly wasted work: deduplicating by unique value cuts
    # the exp count from ``n_hkl * n_atoms`` to ``n_hkl * n_unique``.  The result
    # is numerically identical -- the same value is simply gathered instead of
    # recomputed.
    use_dw = bool(np.any(b_iso > 0.0))
    u_occ = np.unique(occ)
    occ_uniform = u_occ.size == 1
    if use_dw:
        u_b, b_inv = np.unique(b_iso, return_inverse=True)
        dw_uniform = u_b.size == 1
    else:
        u_b = b_inv = None
        dw_uniform = True

    out = np.zeros(hkl.shape[0], dtype=np.float64)
    frac_t = frac.T
    for start in range(0, hkl.shape[0], chunk):
        sl = slice(start, min(start + chunk, hkl.shape[0]))
        s2c = s2[sl][:, None]                                             # (m, 1)
        # Cromer-Mann: f(s) = Z - 41.78214 * s^2 * sum_i a_i exp(-b_i s^2)
        terms = coeff[None, :, :, 0] * np.exp(-coeff[None, :, :, 1] * s2c[:, :, None])
        f_el = zarr[None, :] - _CM_PREFACTOR * s2c * terms.sum(axis=2)    # (m, n_uniq)
        if occ_uniform:
            # Fold the shared occupancy into the per-element factor, so the
            # (m, n_atoms) array is only ever materialised once below.
            f_el *= float(u_occ[0])
        w = f_el[:, el_of_atom]                                           # (m, n_atoms)
        if not occ_uniform:
            w *= occ[None, :]
        if use_dw:
            if dw_uniform:
                w *= np.exp(-float(u_b[0]) * s2c)                         # (m, 1)
            else:
                w *= np.exp(-u_b[None, :] * s2c)[:, b_inv]                # (m, n_uniqB)
        phase = (hkl_f[sl] @ frac_t) * (2.0 * math.pi)
        re = (w * np.cos(phase)).sum(axis=1)
        im = (w * np.sin(phase)).sum(axis=1)
        out[sl] = re * re + im * im
    return out


def _lorentz_polarization(two_theta_deg: np.ndarray) -> np.ndarray:
    """Unpolarised Lorentz-polarisation factor.

    ``LP = (1 + cos^2 2 theta) / (sin^2 theta * cos theta)``.  Reflections with
    ``2 theta >= 180 deg`` are unmeasurable and are never passed in.
    """
    theta = np.radians(two_theta_deg) / 2.0
    # A floor (not a clip of the angle) keeps every normal angle bit-identical:
    # np.maximum returns the original value whenever the trig factor is above
    # it.
    sin_t = np.maximum(np.sin(theta), _MIN_LP_TRIG)
    cos_t = np.maximum(np.cos(theta), _MIN_LP_TRIG)
    # ``np.asarray`` only re-establishes the concrete dtype for numpy's loose
    # ufunc stubs; the value is already a float64 array, so behaviour is
    # unchanged.
    return np.asarray(
        (1.0 + np.cos(2.0 * theta) ** 2) / (sin_t * sin_t * cos_t),
        dtype=np.float64,
    )


def _group_start(two_theta_sorted: np.ndarray, tol: float) -> np.ndarray:
    """Boolean mask marking the first element of every merged group."""
    new_group = np.empty(two_theta_sorted.size, dtype=bool)
    new_group[0] = True
    if two_theta_sorted.size > 1:
        np.greater(np.diff(two_theta_sorted), tol, out=new_group[1:])
    return new_group


def _merge_by_two_theta(
    two_theta: np.ndarray,
    intensity: np.ndarray,
    d: np.ndarray,
    tol: float = TWO_THETA_TOL,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Merge reflections closer than ``tol`` degrees, **summing** intensities.

    This is the step that injects multiplicity into the pattern.  Returns
    ``(d, i, rep_index)`` ordered by ascending 2 theta (= descending d); the d
    value of a group is the one of its first (lowest-angle) member and
    ``rep_index`` indexes that member in the input arrays.
    """
    order = np.argsort(two_theta, kind="stable")
    tt = two_theta[order]
    ii = intensity[order]
    dd = d[order]
    new_group = _group_start(tt, tol)
    group_id = np.cumsum(new_group) - 1
    sums = np.zeros(int(group_id[-1]) + 1, dtype=np.float64)
    np.add.at(sums, group_id, ii)
    return dd[new_group], sums, order[new_group]


def compute_pattern_from_structure(
    ss: gemmi.SmallStructure,
    *,
    dmin: float = DEFAULT_DMIN,
    dmax: float = DEFAULT_DMAX,
    wavelength: float = DEFAULT_WAVELENGTH,
    max_peaks: int = DEFAULT_MAX_PEAKS,
    max_hkl: int = DEFAULT_MAX_HKL,
    return_hkl: bool = False,
    normalize: bool = True,
) -> PatternResult | None:
    """Compute the d-I powder pattern of an already parsed small structure.

    Args:
        ss: Parsed small structure (cell and sites must be present).
        dmin: Lowest d-spacing kept, in A.
        dmax: Largest d-spacing kept, in A.
        wavelength: X-ray wavelength in A.
        max_peaks: Maximum number of merged reflections returned, keeping the
            **strongest** ones.  ``0`` disables the truncation.
        max_hkl: Give up (return ``None``) above this many reflections.
        return_hkl: When true, also return the representative Miller index of
            every merged group as an ``(N, 3)`` int64 array.
        normalize: When true, scale the intensities so the strongest reflection
            equals :data:`INTENSITY_SCALE`.

    Returns:
        ``(d_arr, i_arr)`` or ``(d_arr, i_arr, hkl_arr)``: float64 arrays of
        d-spacings (A) and intensities, ordered by **descending d** and
        truncated to ``max_peaks``.  ``None`` when the structure carries no
        cell / sites or yields no usable reflection.
    """
    if ss.cell is None or not _cell_is_usable(ss.cell) or not ss.sites:
        return None
    cell = ss.cell
    recip = _reciprocal_matrix(cell)
    if recip is None:
        return None
    expanded = _expand_sites(ss)
    if expanded is None:
        return None
    frac, occ, b_iso, symbols = expanded

    enum = _enumerate_hkl(cell, dmin=dmin, dmax=dmax, max_hkl=max_hkl)
    if enum is None:
        return None
    hkl, d = enum

    # 2 d < lambda  <=>  2 theta >= 180 deg: unmeasurable, and LP diverges there.
    # The relative epsilon is essential: for the COD phases whose cell axis is an
    # integer multiple of lambda/2 (e.g. c = 7.703 A = 10 * 0.7703) the boundary
    # reflection is kept by a naive ``2*d > lambda`` test through float rounding,
    # its LP factor explodes to 1e16..1e21, and the normalised pattern collapses
    # to that single bogus line at d = 0.7703 A.
    measurable = (2.0 * d) > wavelength * (1.0 + _MEASURABLE_REL_EPS)
    if not measurable.any():
        return None
    hkl = hkl[measurable]
    d = d[measurable]

    g_abs = np.linalg.norm(hkl.astype(np.float64) @ recip, axis=1)
    two_theta = np.degrees(2.0 * np.arcsin(np.clip(wavelength * g_abs / 2.0, -1.0, 1.0)))

    f2 = _structure_factors(hkl, frac, occ, b_iso, symbols, recip)
    intensity = f2 * _lorentz_polarization(two_theta)

    d_m, i_m, rep = _merge_by_two_theta(two_theta, intensity, d)
    if i_m.size == 0:
        return None
    imax = float(np.max(i_m))
    if not math.isfinite(imax) or imax <= 0.0:
        return None
    good = np.isfinite(i_m) & (i_m >= imax * REL_INTENSITY_CUTOFF)
    if not good.any():
        return None
    d_m = d_m[good]
    i_m = i_m[good]
    if return_hkl:
        rep_hkl = hkl[rep[good]]

    if normalize:
        peak = float(np.max(i_m))
        if peak > 0.0:
            i_m = i_m * (INTENSITY_SCALE / peak)

    if max_peaks and i_m.size > max_peaks:
        # Keep the strongest, then restore descending-d (= ascending 2 theta).
        sel = np.argsort(-i_m, kind="stable")[:max_peaks]
        sel.sort()
        d_m = d_m[sel]
        i_m = i_m[sel]
        if return_hkl:
            rep_hkl = rep_hkl[sel]

    if return_hkl:
        return d_m, i_m, rep_hkl.astype(np.int64)
    return d_m, i_m


def _compute_pattern_recovering_cell(
    ss: gemmi.SmallStructure,
    block: gemmi.cif.Block | None,
    *,
    dmin: float,
    dmax: float,
    wavelength: float,
    max_peaks: int,
    max_hkl: int,
    return_hkl: bool,
    normalize: bool,
) -> PatternResult | None:
    """Shared cell-recovery step of every file/text/block pattern entry point.

    gemmi leaves its default ``(1, 1, 1, 90, 90, 90)`` cell in place when it
    cannot map the CIF cell tags, which would produce one bogus d = 1 A line.
    The cell is therefore rebuilt from the CIF text tags (a few COD files split
    the six values over several items); when that fails too the phase is
    reported as uncomputable.  Keeping this in one function is what makes
    :func:`compute_pattern`, :func:`compute_pattern_from_text` and
    :func:`compute_pattern_from_block` agree bit for bit.
    """
    if not _cell_is_usable(ss.cell):
        recovered = _recover_cell(block)
        if recovered is None:
            return None
        ss.cell = recovered
    return compute_pattern_from_structure(
        ss,
        dmin=dmin,
        dmax=dmax,
        wavelength=wavelength,
        max_peaks=max_peaks,
        max_hkl=max_hkl,
        return_hkl=return_hkl,
        normalize=normalize,
    )


def compute_pattern_from_block(
    block: gemmi.cif.Block | None,
    *,
    structure: gemmi.SmallStructure | None = None,
    dmin: float = DEFAULT_DMIN,
    dmax: float = DEFAULT_DMAX,
    wavelength: float = DEFAULT_WAVELENGTH,
    max_peaks: int = DEFAULT_MAX_PEAKS,
    max_hkl: int = DEFAULT_MAX_HKL,
    return_hkl: bool = False,
    normalize: bool = True,
) -> PatternResult | None:
    """Compute the d-I powder pattern of a CIF block the caller already parsed.

    Block-based twin of :func:`compute_pattern`, used by the library builder so
    a single ``gemmi.cif.read`` feeds both the metadata and the peak table.  It
    applies exactly the same placeholder-cell recovery and physics as the
    path-based reader; :func:`compute_pattern_from_structure` documents the
    parameters.

    Args:
        block: The sole CIF block of the file (from ``gemmi.cif.read`` /
            ``gemmi.cif.read_string``), or ``None`` when the file could not be
            parsed.
        structure: The structure already built from *block*.  Passing it in
            avoids a second ``make_small_structure_from_block`` call; when
            omitted, it is built here.

    Returns:
        ``(d_arr, i_arr)`` (or the ``return_hkl`` variant), or ``None`` when the
        block yields no usable structure / reflection.
    """
    if structure is None:
        if block is None:
            return None
        try:
            structure = gemmi.make_small_structure_from_block(block)
        except Exception:
            return None
    return _compute_pattern_recovering_cell(
        structure,
        block,
        dmin=dmin,
        dmax=dmax,
        wavelength=wavelength,
        max_peaks=max_peaks,
        max_hkl=max_hkl,
        return_hkl=return_hkl,
        normalize=normalize,
    )


def compute_pattern_from_text(text: str, **kwargs) -> PatternResult | None:
    """Compute the pattern of CIF text held in memory (e.g. a ``cif_gz`` blob)."""
    try:
        block = gemmi.cif.read_string(text).sole_block()
    except Exception:
        return None
    return compute_pattern_from_block(block, **kwargs)


def compute_pattern(
    cif_path: str | Path,
    *,
    dmin: float = DEFAULT_DMIN,
    dmax: float = DEFAULT_DMAX,
    wavelength: float = DEFAULT_WAVELENGTH,
    max_peaks: int = DEFAULT_MAX_PEAKS,
    max_hkl: int = DEFAULT_MAX_HKL,
    return_hkl: bool = False,
    normalize: bool = True,
) -> PatternResult | None:
    """Compute the d-I powder pattern of a CIF file.

    A single unreadable file must never abort a whole library build, so every
    parse failure is reported as ``None``.  A phase whose cell gemmi cannot map
    is recovered from the CIF text tags; when even that fails the phase is
    uncomputable and ``None`` is returned instead of a pattern computed from
    gemmi's default (1 A) cell.  See :func:`compute_pattern_from_structure` for
    the parameters.
    """
    try:
        ss = gemmi.read_small_structure(str(cif_path))
    except Exception:
        return None
    # Only re-read the text when gemmi could not map the cell: the common case
    # must not pay for a second parse.
    block = None if _cell_is_usable(ss.cell) else _read_cif_block(cif_path)
    return _compute_pattern_recovering_cell(
        ss,
        block,
        dmin=dmin,
        dmax=dmax,
        wavelength=wavelength,
        max_peaks=max_peaks,
        max_hkl=max_hkl,
        return_hkl=return_hkl,
        normalize=normalize,
    )
