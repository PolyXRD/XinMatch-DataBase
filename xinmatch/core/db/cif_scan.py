"""COD CIF archive scanner for the v3.0.0 database rebuild (M1-S1.1).

What this module does
---------------------
The v3.0.0 rebuild has to walk the local COD CIF archive twice: once to
count what is there, and once to feed every CIF into the parser.  Both jobs
need the same three facts:

* how a numeric COD id maps onto a path inside the archive;
* how a numeric COD id maps onto the reference-software entry id
  (``96-XXX-XXXX``) used throughout the existing XinMatch databases;
* how to enumerate the archive **in ascending COD id order** without
  materialising 535 K paths in memory.

Archive layout (verified against the real archive)
--------------------------------------------------
The archive is a fan-out tree whose directory names are slices of the
7-digit, zero-padded COD id::

    1/00/00/1000005.cif      <- cod_id 1000005
    9/01/02/9010203.cif      <- cod_id 9010203

The first path component is the first id digit, the second is digits 2-3,
the third is digits 4-5, and the file name carries the full 7-digit id.
Each leaf directory therefore holds exactly 100 consecutive ids.

Because every path component is fixed width and zero padded, sorting the
components lexicographically is identical to sorting the COD ids
numerically.  The scanner relies on that: it walks the tree level by level
and sorts each level by name, so the files come out in ascending COD id
order for free, without a global sort over half a million entries.

The archive root also contains a ``dictionaries/`` directory holding the
CIF *dictionary* definitions (``.dic`` files).  Those are schema
definitions, not crystal structures, and are therefore skipped by default:
they would otherwise pollute both the file count and the rebuild input.

Why not just glob ``**/*.cif``
------------------------------
``Path.rglob`` on this archive is both slow (it re-walks the same
directories) and unordered (it would require a full sort of 535 K paths to
restore the id order).  The directory names already encode the ordering
information, so walking the tree in name order is both faster and
order-correct by construction.
"""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from pathlib import Path

logger = logging.getLogger(__name__)

__all__ = [
    "DEFAULT_COD_ROOT",
    "DEFAULT_INORG_COD_ROOT",
    "DEFAULT_INORG_PCOD_ROOT",
    "DEFAULT_INORG_TCOD_ROOT",
    "cod_cif_relpath",
    "cod_id_from_stem",
    "entry_id_for_cod",
    "iter_cod_cifs",
    "count_cif_files",
]

#: 共享数据根：COD/PCOD/TCOD 原始归档、``inorg_cif`` 无机子集与
#: ``inorg_manifest.sqlite`` 清单都在这里，由 PolyXRD 与 XinMatch 两库共用。
#: 2026-10-23 起移出仓库、统一放 ``D:/Project/XRD/sync-COD``；
#: 换盘/换机器用环境变量 ``XRD_DATA_ROOT`` 覆盖。
DATA_ROOT = Path(os.environ.get("XRD_DATA_ROOT") or r"D:\Project\XRD\sync-COD")

#: Full COD CIF archive root (read-only, 535 166 files).
DEFAULT_COD_ROOT = DATA_ROOT / "COD"

#: Inorganic-only subset of the COD archive (read-only, 79 198 files).
DEFAULT_INORG_COD_ROOT = DATA_ROOT / "inorg_cif" / "COD"

#: Inorganic-only subset of the PCOD archive (predicted structures, 1 061 962 files).
DEFAULT_INORG_PCOD_ROOT = DATA_ROOT / "inorg_cif" / "PCOD"

#: Inorganic-only subset of the TCOD archive (theoretical structures, 46 943 files).
DEFAULT_INORG_TCOD_ROOT = DATA_ROOT / "inorg_cif" / "TCOD"

#: Number of digits in a COD/PCOD id (those archives zero-pad to 7).
_ID_DIGITS = 7

#: Largest id the 7-digit COD archive layout can express.
_MAX_COD_ID = 10**_ID_DIGITS - 1

#: Largest id :func:`entry_id_for_cod` accepts.  Deliberately wider than the
#: archive width: TCOD ids are 8 digits, and the 3COD convention shifts PCOD by
#: +20 000 000 and TCOD by +100 000 000, which pushes stored ids to 9 digits.
#: Pinning the formatter to 7 digits would reject every PCOD/TCOD phase.
_MAX_ENTRY_ID = 10**9 - 1

#: CIF file suffix (compared case-insensitively).
_CIF_SUFFIX = ".cif"

#: Directory holding CIF dictionary definitions; never a structure entry.
_DICTIONARIES_DIRNAME = "dictionaries"


# ---------------------------------------------------------------------------
# Pure id / path helpers
# ---------------------------------------------------------------------------
def _cod_id_text(cod_id: int, *, digits: int = _ID_DIGITS) -> str:
    """Return the zero-padded text of *cod_id*, at least *digits* wide.

    Args:
        cod_id: Numeric phase id.
        digits: Minimum width.  Defaults to the 7 digits the COD and PCOD
            archives use; longer ids are returned as-is rather than truncated,
            which is what lets TCOD's 8-digit ids and the shifted 9-digit entry
            ids go through the same path.

    Returns:
        str: The id zero-padded to *digits*, e.g. ``"1000005"``.

    Raises:
        ValueError: If *cod_id* is negative or wider than
            :data:`_MAX_ENTRY_ID`.  Failing loudly beats silently producing an
            id that cannot be mapped back to a file.
    """
    value = int(cod_id)
    if value < 0 or value > _MAX_ENTRY_ID:
        raise ValueError(
            f"phase id out of range: {cod_id} (expected 0..{_MAX_ENTRY_ID})"
        )
    return str(value).zfill(digits)


def cod_cif_relpath(cod_id: int) -> str:
    """Numeric COD id to its archive-relative CIF path.

    The layout splits the 7-digit, zero-padded id as
    ``{first}/{next two}/{next two}/{id:07d}.cif`` — e.g. id ``1000005``
    becomes ``"1/00/00/1000005.cif"``.

    Only the **COD/PCOD 7-digit layout** is described here.  TCOD uses an
    8-digit id with its own fan-out, so it is not addressable through this
    helper; :func:`iter_cod_cifs` discovers paths by walking the tree instead
    and works for all three archives.

    Args:
        cod_id: Numeric COD id.

    Returns:
        str: Relative path using POSIX ``/`` separators, independent of the
        host platform.

    Raises:
        ValueError: If *cod_id* does not fit the 7-digit layout.
    """
    s = _cod_id_text(cod_id)
    if len(s) != _ID_DIGITS:
        raise ValueError(
            f"cod_id {cod_id} does not fit the {_ID_DIGITS}-digit COD layout"
        )
    return f"{s[0]}/{s[1:3]}/{s[3:5]}/{s}.cif"


def cod_id_from_stem(stem: str) -> int | None:
    """Parse a COD id out of a CIF file name.

    Accepts both a bare stem (``"1000005"``) and a full file name
    (``"1000005.cif"``); the ``.cif`` suffix is optional.

    Args:
        stem: File name or stem to parse.

    Returns:
        int | None: The numeric COD id, or ``None`` when the name is not
        purely numeric.  Non-numeric names are expected in the archive's
        ``dictionaries/`` subtree and in hand-dropped files, so they are
        reported as "not a COD entry" instead of raising.
    """
    name = stem.strip()
    if name.lower().endswith(_CIF_SUFFIX):
        name = name[: -len(_CIF_SUFFIX)]
    # ``str.isdigit`` also accepts non-ASCII digits and superscripts that
    # ``int()`` would reject, so the ASCII check is load-bearing.
    if not name.isascii() or not name.isdigit():
        return None
    return int(name)


def entry_id_for_cod(cod_id: int) -> str:
    """Numeric phase id to the reference-software style entry id.

    The entry id keeps the PolyXRD ``ref_id`` convention
    (``"96-400-1849"``) that every existing XinMatch database and every
    baseline CSV already uses: the zero-padded id is split into its first
    three and the remaining digits and joined with ``-`` behind the ``96-``
    prefix.

    The split is what the whole ecosystem keys on, so it is applied verbatim
    to PCOD/TCOD ids as well — those are shifted into their own numeric range
    first (``+20 000 000`` / ``+100 000 000``, the 3COD convention), which
    keeps COD ``1000000`` and PCOD ``1000000`` apart:

    ========  ==========  ============= ==========================
    archive   numeric id  stored id     entry id
    ========  ==========  ============= ==========================
    COD       1000005     1000005       ``96-100-0005``
    PCOD      1000000     21000000      ``96-210-00000``
    TCOD      10000001    110000001     ``96-110-000001``
    ========  ==========  ============= ==========================

    Args:
        cod_id: Numeric id in the library's own (already shifted) range.

    Returns:
        str: Entry id, e.g. ``"96-100-0005"`` for id ``1000005``.

    Raises:
        ValueError: If *cod_id* is out of range (see
            :data:`_MAX_ENTRY_ID`).
    """
    s = _cod_id_text(cod_id)
    return f"96-{s[:3]}-{s[3:]}"


# ---------------------------------------------------------------------------
# Archive traversal
# ---------------------------------------------------------------------------
def _subtree_below_start(prefix: str | None, name: str, start: int | None) -> bool:
    """True when every id under the *name* subtree is certainly below *start*.

    Used only to prune the walk, so a wrong ``True`` silently drops phases --
    hence the decision is made from the leading digits alone and never from an
    assumed id width.  Because ids are stored zero-padded and each directory
    name is a slice of that spelling, comparing the accumulated prefix against
    the head of ``str(start)`` is exactly a numeric comparison.

    A width-based guess is what used to be here, and it was wrong for TCOD:
    with 8-digit ids the "remaining digits" arithmetic underflows and either
    prunes a populated subtree or fails to prune at all.

    Args:
        prefix: Id digits contributed by the ancestor directories, or ``None``
            when an ancestor name was not numeric.
        name: Name of the sub-directory.
        start: Requested lowest id, or ``None`` for "no lower bound".

    Returns:
        bool: ``True`` when the subtree can be skipped.
    """
    if start is None or prefix is None:
        return False
    if not name.isascii() or not name.isdigit():
        return False
    node = prefix + name
    start_text = str(int(start))
    head = start_text[: len(node)]
    if len(head) < len(node):
        # The subtree prefix is already wider than ``start`` itself, so the
        # comparison cannot narrow anything down -- walk it.
        return False
    return node < head


def _iter_cif_files(
    directory: Path,
    prefix: str | None,
    *,
    include_dictionaries: bool,
    start: int | None,
) -> Iterator[tuple[int, Path]]:
    """Recursively yield ``(cod_id, path)`` for the CIF files under *directory*.

    Args:
        directory: Directory to scan.
        prefix: Numeric id digits contributed by the ancestor directory
            names, or ``None`` when the prefix is unknown.
        include_dictionaries: Whether to descend into ``dictionaries/``.
        start: Only yield entries whose COD id is ``>= start``.

    Yields:
        tuple[int, Path]: The COD id and the path of one CIF file.
    """
    try:
        with os.scandir(directory) as scan:
            entries = list(scan)
    except OSError as exc:
        logger.warning("cannot scan CIF directory %s: %s", directory, exc)
        return

    # Sort by name: the archive's components are fixed-width and zero
    # padded, so lexicographic order equals ascending COD id order.
    entries.sort(key=lambda entry: entry.name)

    for entry in entries:
        try:
            if entry.is_dir(follow_symlinks=False):
                if (
                    not include_dictionaries
                    and entry.name == _DICTIONARIES_DIRNAME
                ):
                    continue
                if start is not None and _subtree_below_start(prefix, entry.name, start):
                    # The whole subtree lies below the requested start id;
                    # skipping it here is what makes ``start`` cheap.
                    continue
                name = entry.name
                child_prefix = (
                    prefix + name
                    if prefix is not None and name.isascii() and name.isdigit()
                    else None
                )
                yield from _iter_cif_files(
                    Path(entry.path),
                    child_prefix,
                    include_dictionaries=include_dictionaries,
                    start=start,
                )
            elif entry.is_file(follow_symlinks=False) and entry.name.lower().endswith(
                _CIF_SUFFIX
            ):
                cod_id = cod_id_from_stem(entry.name)
                if cod_id is None:
                    # Not a COD entry (e.g. a stray file dropped into the
                    # archive); it cannot be addressed by id, so skip it.
                    continue
                if start is not None and cod_id < start:
                    continue
                yield cod_id, Path(entry.path)
        except OSError as exc:
            logger.warning("cannot inspect %s: %s", entry.path, exc)


def iter_cod_cifs(
    root: str | Path,
    *,
    start: int | None = None,
    limit: int | None = None,
    include_dictionaries: bool = False,
) -> Iterator[tuple[int, Path]]:
    """Iterate over the CIF files of a COD archive in ascending COD id order.

    The archive is expected to follow the standard COD fan-out layout
    (``{d}/{dd}/{dd}/{id:07d}.cif``); the walk is ordered by the directory
    names, which are slices of the zero-padded id, so the yield order is
    ascending by COD id without a global sort.

    Args:
        root: Archive root, e.g. :data:`DEFAULT_COD_ROOT`.
        start: Only yield entries whose COD id is ``>= start``.  Useful to
            resume an interrupted rebuild without re-scanning the skipped
            part of the archive.
        limit: Stop after yielding this many entries.  ``None`` means "no
            limit"; a value ``<= 0`` yields nothing.
        include_dictionaries: When ``True``, descend into the archive's
            ``dictionaries/`` subtree as well.  Defaults to ``False``
            because those files are CIF *dictionary* definitions, not
            crystal structures.

    Yields:
        tuple[int, Path]: The COD id and the path of one CIF file.  The
        path is always an existing regular file.

    Raises:
        FileNotFoundError: If *root* does not exist.
    """
    root_path = Path(root)
    if not root_path.is_dir():
        raise FileNotFoundError(f"COD archive root not found: {root_path}")
    if limit is not None and limit <= 0:
        return

    for index, (cod_id, path) in enumerate(
        _iter_cif_files(
            root_path,
            "",
            include_dictionaries=include_dictionaries,
            start=start,
        )
    ):
        yield cod_id, path
        if limit is not None and index + 1 >= limit:
            return


def count_cif_files(root: str | Path) -> int:
    """Count the CIF files of a COD archive.

    Args:
        root: Archive root, e.g. :data:`DEFAULT_COD_ROOT`.

    Returns:
        int: Number of CIF files, excluding the ``dictionaries/`` subtree.
        This is exactly the number of entries :func:`iter_cod_cifs` yields
        with its default arguments, so the two can be used interchangeably
        for progress reporting.
    """
    return sum(1 for _ in iter_cod_cifs(root))
