"""Materials Project CIF archive scanning and the MP entry-id convention.

The Materials Project summary archive is a **flat** directory of CIF files named
after the formula (``cif/Ac.cif``, ``cif/Ac2CdSn.cif``, ...), so the file name
is not an id.  Every file starts with a comment block::

    # Materials Project summary 2025-09-25
    # mp-id: mp-862690 | spacegroup: P6_3/mmc | energy_above_hull: 0.000000 eV/atom
    # generated using pymatgen

The ``mp-id`` is stable and unique, so the library uses it as ``entry_id``
(``mp-862690``).  This differs from the COD-family ``96-XXX-XXXX`` ids on
purpose: MP is an external database and its native identifier is the more
useful provenance.

The archive is scanned **once** (155k small header reads) to build an ordered
``(mp_id, path)`` list; the builder and the acceptance tool share
:func:`scan_mp_archive` so they agree on which files are in scope.
"""

from __future__ import annotations

import os
import re
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

__all__ = [
    "DEFAULT_MP_ROOT",
    "MP_ID_RE",
    "MP_JOURNAL",
    "MP_QUALITY",
    "MpArchive",
    "MpHeader",
    "crystal_system_from_symbol",
    "entry_id_for_mp",
    "iter_mp_cifs",
    "parse_mp_header",
    "scan_mp_archive",
]

#: Materials Project CIF archive root (``MP_CIF_ROOT`` overrides it).
DEFAULT_MP_ROOT = Path(
    os.environ.get("MP_CIF_ROOT") or r"D:\Project\MP-summary\cif"
)

#: ``entry_id`` prefix used by this external library.
MP_PREFIX = "mp-"

#: Provenance strings written into ``entries.journal/author/quality``.
MP_JOURNAL = "Materials Project"
MP_AUTHOR = "Materials Project"
MP_QUALITY = "DFT"

MP_ID_RE = re.compile(r"mp-id:\s*mp-(\d+)", re.IGNORECASE)
_SPACEGROUP_RE = re.compile(r"spacegroup:\s*([^|]+)", re.IGNORECASE)
_ENERGY_RE = re.compile(r"energy_above_hull:\s*(-?[0-9.]+)", re.IGNORECASE)
_SUMMARY_RE = re.compile(r"Materials Project summary\s+(\d{4})-(\d{2})-(\d{2})", re.IGNORECASE)

#: Largest mp id accepted by :func:`entry_id_for_mp` (MP currently < 3.4e6).
_MAX_MP_ID = 10**9 - 1

#: International space-group number -> crystal system (same ranges as
#: :mod:`xinmatch.core.db.cif_meta`).
_CRYSTAL_SYSTEMS: tuple[tuple[int, str], ...] = (
    (2, "triclinic"),
    (15, "monoclinic"),
    (74, "orthorhombic"),
    (142, "tetragonal"),
    (167, "trigonal"),
    (194, "hexagonal"),
    (230, "cubic"),
)


@dataclass(frozen=True)
class MpHeader:
    """Parsed comment header of one MP CIF."""

    mp_id: int
    space_group: str | None = None
    energy_above_hull: float | None = None
    summary_date: str | None = None


@dataclass(frozen=True)
class MpArchive:
    """Result of :func:`scan_mp_archive`."""

    #: ``(mp_id, path)`` sorted by ``mp_id``.
    entries: tuple[tuple[int, Path], ...]
    #: Files that carry no ``mp-id`` comment (reported, never silently dropped).
    no_id: tuple[Path, ...]
    #: Total ``*.cif`` files seen (top level only).
    total: int


def parse_mp_header(path: str | Path) -> MpHeader | None:
    """Read the leading comment block of one MP CIF; ``None`` when no mp-id."""
    try:
        with Path(path).open("r", encoding="utf-8", errors="replace") as f:
            head = ""
            for _ in range(8):
                line = f.readline()
                if not line:
                    break
                head += line
                if "mp-id:" in line:
                    break
    except OSError:
        return None
    match = MP_ID_RE.search(head)
    if match is None:
        return None
    sg = _SPACEGROUP_RE.search(head)
    energy = _ENERGY_RE.search(head)
    summary = _SUMMARY_RE.search(head)
    date = f"{summary.group(1)}-{summary.group(2)}-{summary.group(3)}" if summary else None
    return MpHeader(
        mp_id=int(match.group(1)),
        space_group=sg.group(1).strip() if sg else None,
        energy_above_hull=float(energy.group(1)) if energy else None,
        summary_date=date,
    )


def scan_mp_archive(root: str | Path = DEFAULT_MP_ROOT) -> MpArchive:
    """Scan the flat MP CIF directory once and return id-sorted entries.

    Raises:
        FileNotFoundError: When *root* does not exist.
    """
    root_path = Path(root)
    if not root_path.is_dir():
        raise FileNotFoundError(f"MP CIF archive root not found: {root_path}")
    files = sorted(root_path.glob("*.cif"))
    entries: list[tuple[int, Path]] = []
    no_id: list[Path] = []
    for path in files:
        header = parse_mp_header(path)
        if header is None:
            no_id.append(path)
        else:
            entries.append((header.mp_id, path))
    entries.sort(key=lambda item: item[0])
    return MpArchive(entries=tuple(entries), no_id=tuple(no_id), total=len(files))


def iter_mp_cifs(
    root: str | Path = DEFAULT_MP_ROOT,
    *,
    start: int | None = None,
    limit: int | None = None,
    archive: MpArchive | None = None,
) -> Iterator[tuple[int, Path]]:
    """Yield ``(mp_id, path)`` in ascending id order.

    Args:
        root: MP CIF directory.
        start: Only yield ``mp_id >= start`` (resume support).
        limit: Stop after this many entries.
        archive: Pre-scanned :class:`MpArchive` (avoids a second 155k-file scan).
    """
    if limit is not None and limit <= 0:
        return
    scanned = archive if archive is not None else scan_mp_archive(root)
    if start is not None and limit is None:
        # Cheap binary search would also work; the list is already sorted.
        scanned_entries = scanned.entries
        lo, hi = 0, len(scanned_entries)
        while lo < hi:
            mid = (lo + hi) // 2
            if scanned_entries[mid][0] < start:
                lo = mid + 1
            else:
                hi = mid
        scanned_entries = scanned_entries[lo:]
    else:
        scanned_entries = scanned.entries
    for index, (mp_id, path) in enumerate(scanned_entries):
        yield mp_id, path
        if limit is not None and index + 1 >= limit:
            return


def entry_id_for_mp(mp_id: int) -> str:
    """Return the canonical entry id of an MP phase (``mp-862690``)."""
    value = int(mp_id)
    if value < 0 or value > _MAX_MP_ID:
        raise ValueError(f"MP id out of range: {mp_id} (expected 0..{_MAX_MP_ID})")
    return f"{MP_PREFIX}{value}"


def crystal_system_from_symbol(symbol: str | None) -> str | None:
    """Crystal system of an H-M symbol like ``P6_3/mmc`` (``None`` if unknown)."""
    if not symbol:
        return None
    import gemmi

    text = str(symbol).strip().strip("'\"")
    for candidate in (text, text.replace("_", ""), text.replace("_", " ")):
        try:
            space_group = gemmi.find_spacegroup_by_name(candidate)
        except Exception:
            space_group = None
        if space_group is not None:
            number = int(space_group.number)
            for upper, system in _CRYSTAL_SYSTEMS:
                if number <= upper:
                    return system
    return None
