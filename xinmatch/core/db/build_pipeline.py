"""One-stop phase computation for the v3 ``.xinDB`` builder (single CIF read).

The builder needs two products per CIF file: the scalar :class:`PhaseMeta`
record and the ``(d, i)`` powder pattern.  Historically these were produced by
calling two independent helpers::

    meta  = parse_meta(cif_path, entry_id=...)
    peaks = compute_pattern(cif_path, dmin=..., max_peaks=...)

Each of those helpers re-opened and re-parsed the same file, so a library build
parsed every CIF two to four times.  :func:`compute_phase` is the builder-facing
entry point that produces both results from **one** parse -- the full six-library
rebuild touches half a million files per archive, so parsing once is the
cheapest throughput win available without touching the frozen physics.

How the single parse stays bit-identical
----------------------------------------
The two products are computed by the *same* functions the old two-step path
used, fed from one parse:

1. ``block = gemmi.cif.read(path).sole_block()`` -- the one and only CIF parse;
2. ``structure = gemmi.make_small_structure_from_block(block)`` -- the one and
   only structure build (``gemmi.read_small_structure`` is a thin wrapper around
   exactly these two calls, so the structure is identical);
3. :func:`~xinmatch.core.db.cif_meta.parse_meta_from_structure` -- the shared
   core of ``parse_meta``, so element fallbacks / cell recovery / space-group
   normalisation are reused, not re-implemented;
4. :func:`~xinmatch.core.db.peak_compute.compute_pattern_from_block` (with the
   already built structure) -- the block-based twin of ``compute_pattern``,
   including its placeholder-cell recovery.

Semantics: :func:`compute_phase` is a drop-in replacement for the two calls
above, in the same order and with the same parameters.  A caller that checks
``meta is None`` first and ``peaks is None`` second sees exactly the same
values -- metadata fields and peak arrays are bit-identical, including the
``reason`` string that mirrors the old ``_compute_one`` failure markers:

======================  =====================================================
reason                  meaning (checked in this order)
======================  =====================================================
``""``                  success: both ``meta`` and ``(d, i)`` are returned
``"no_cif"``            the path does not exist / is not a file
``"no_meta"``           the CIF could not be parsed (old ``parse_meta`` None)
``"no_peaks"``          metadata is fine but no usable reflection (old
                        ``compute_pattern`` None: no cell, no sites, ...)
======================  =====================================================

On every failure the function returns ``(None, None, reason)``, mirroring the
old two-call builder loop which discarded the metadata of a phase without peaks.

The module only depends on :mod:`numpy` / :mod:`gemmi` (through
``cif_meta`` and ``peak_compute``), so it is safe to import from a worker
process of the builder without dragging in Qt or scipy.
"""

from __future__ import annotations

from pathlib import Path
from typing import TypeAlias

import gemmi
import numpy as np

from .cif_meta import PhaseMeta, parse_meta_from_structure
from .peak_compute import compute_pattern_from_block

__all__ = [
    "PeakArrays",
    "REASON_NO_CIF",
    "REASON_NO_META",
    "REASON_NO_PEAKS",
    "REASON_OK",
    "compute_phase",
]

#: The ``(d, intensity)`` arrays of one phase: both ``float64``, ordered by
#: descending d-spacing.
PeakArrays: TypeAlias = tuple[np.ndarray, np.ndarray]

#: Why :func:`compute_phase` returned no phase.  ``""`` means success.
REASON_OK = ""
REASON_NO_CIF = "no_cif"
REASON_NO_META = "no_meta"
REASON_NO_PEAKS = "no_peaks"


def compute_phase(
    cif_path: Path,
    *,
    entry_id: str,
    dmin: float,
    max_peaks: int,
    normalize: bool = True,
) -> tuple[PhaseMeta | None, PeakArrays | None, str]:
    """Compute the metadata and the powder pattern of one CIF file.

    This is the builder-facing single entry point: it replaces the old
    ``parse_meta`` + ``compute_pattern`` pair and is bit-for-bit equivalent to it
    for the same arguments (see the module docstring), while reading the CIF only
    once.

    Args:
        cif_path: Path of the CIF file to read.
        entry_id: XinMatch entry id of the phase, recorded on ``PhaseMeta``.
            Callers own the ``library_kind`` id offset; this function only
            stores the string it is given.
        dmin: Lowest d-spacing (A) kept in the pattern.
        max_peaks: Strongest merged reflections kept per phase (``0`` disables
            the truncation).
        normalize: When true, scale the strongest reflection to
            :data:`~xinmatch.core.db.peak_compute.INTENSITY_SCALE`.

    Returns:
        ``(meta, (d, i), reason)``.  On success ``reason`` is ``""``; otherwise
        ``meta`` and the arrays are ``None`` and ``reason`` is one of
        ``"no_cif"`` / ``"no_meta"`` / ``"no_peaks"`` (checked in that order,
        matching the old builder loop).
    """
    path = Path(cif_path)
    if not path.is_file():
        return None, None, REASON_NO_CIF

    # ---- the single parse (contrast with the old two-step path) -------------
    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except Exception:
        # gemmi raises on malformed / empty CIFs; a single bad file must never
        # abort a whole library build.
        return None, None, REASON_NO_META
    try:
        structure = gemmi.make_small_structure_from_block(block)
    except Exception:
        return None, None, REASON_NO_META

    # ---- metadata: same shared core as the public parse_meta() --------------
    meta = parse_meta_from_structure(structure, block, entry_id=entry_id)
    if meta is None:
        return None, None, REASON_NO_META

    # ---- peaks: same physics (and cell recovery) as compute_pattern() ------
    peaks = compute_pattern_from_block(
        block,
        structure=structure,
        dmin=dmin,
        max_peaks=max_peaks,
        normalize=normalize,
    )
    if peaks is None:
        return None, None, REASON_NO_PEAKS

    return meta, (peaks[0], peaks[1]), REASON_OK
