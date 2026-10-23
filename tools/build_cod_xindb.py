"""Build a v3 ``.xinDB`` reference library from the COD CIF archive (M2 / S2.3 + S2.4).

Two products come out of this one implementation:

``COD-Full.xinDB``
    Built from the **full** COD CIF archive (535,166 files).  This is the
    exhaustive reference: every phase the archive contains.

``COD-Inorg.xinDB``
    Built from the **inorganic subset** archive (79,198 files).  This is
    the default interactive library — the full archive contains thousands
    of near-duplicate cards for the same phase, and those duplicates dilute
    the FoM vote (measured: top-1 agreement collapses as the library grows).

Both are **external** libraries: they are not shipped with the application.
The build writes a self-contained single file that the user can keep in any
folder and mount through *设置 → 参考数据库 → 添加外部库*.

Data flow
---------
::

    CIF archive ──► cif_scan (paths)  ──► cif_meta (scalars)  ──┐
                     │                                            ├──► XinDBWriter ──► .xinDB
                     └──► peak_compute (d, I, hkl) ─────────────┘

Nothing else is read: no legacy database, no intermediate index.  The CIF
archive is the single source of truth, so the result is reproducible from
the archive alone.

Resumability
------------
A full build parses half a million CIF files and computes a powder pattern
for each — hours of work.  The script therefore keeps a watermark table
inside the output file and supports ``--resume`` so an interrupted run
continues instead of starting over.  ``--start`` / ``--limit`` restrict the
work to a slice for verification.

Usage
-----
::

    # verification slice first (always do this before a full build)
    python -m tools.build_cod_xindb --limit 1000 --out /tmp/slice.xinDB

    # full build, resumable
    python -m tools.build_cod_xindb --out D:/XRD-DB/COD-Full.xinDB

    # the inorganic interactive library
    python -m tools.build_cod_inorg_xindb --out D:/XRD-DB/COD-Inorg.xinDB
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from xinmatch.core.db.cif_meta import PhaseMeta  # noqa: E402
from xinmatch.core.db.peak_compute import (  # noqa: E402
    DEFAULT_DMIN,
    DEFAULT_MAX_PEAKS,
    compute_pattern,
)
from xinmatch.core.db.writer_v3 import XinDBWriter  # noqa: E402

from xinmatch.core.db.cif_scan import (  # noqa: E402
    DEFAULT_COD_ROOT,
    DEFAULT_INORG_COD_ROOT,
    iter_cod_cifs,
)
from xinmatch.core.db.schema_v3 import (  # noqa: E402
    CHUNK_SIZE,
    TOP_K,
    XINDB_SUFFIX,
)

logger = logging.getLogger("tools.build_cod_xindb")


def default_out_dir() -> Path:
    """Directory the build scripts write to when ``--out`` is not given.

    ``XINMATCH_XINDB_DIR`` wins, which is the same variable
    :func:`xinmatch.core.db.multi_db.xindb_default_dir` consults, so one
    environment variable keeps the build output and the library picker pointing
    at the same folder -- typically a shared archive next to the CIF sources
    (``D:/Project/XRD/PolyXRD-DataBase/xindb``) instead of scattering one large
    ``.xinDB`` per checkout.  Without the variable the scripts keep writing to
    ``<repo>/data``, so a plain checkout stays self-contained.
    """
    override = os.environ.get("XINMATCH_XINDB_DIR", "").strip()
    if override:
        return Path(override).expanduser()
    return PROJECT_ROOT / "data"

#: The two libraries this module can produce.
LIBRARY_KINDS = {
    "full": (
        DEFAULT_COD_ROOT,
        "COD-Full",
        "COD full archive (all phases in the COD CIF collection)",
    ),
    "inorg": (
        DEFAULT_INORG_COD_ROOT,
        "COD-Inorg",
        "COD inorganic subset (default interactive library)",
    ),
}


@dataclass
class BuildStats:
    """Counters collected during a build."""

    entries: int = 0
    peaks: int = 0
    skipped_no_cif: int = 0
    skipped_no_meta: int = 0
    skipped_no_peaks: int = 0
    errors: int = 0
    elapsed: float = 0.0
    per_minute: float = 0.0
    #: Phases inherited from the interrupted run this build continued (0 = fresh).
    resumed_from: int = 0
    #: Phases the finished file holds (inherited + this run).
    entries_total: int = 0
    skipped_ids: list[int] = field(default_factory=list)

    def as_dict(self) -> dict[str, object]:
        """Return the stats as a plain dict (JSON-friendly)."""
        return {
            "entries": self.entries,
            "entries_total": self.entries_total,
            "resumed_from": self.resumed_from,
            "peaks": self.peaks,
            "skipped_no_cif": self.skipped_no_cif,
            "skipped_no_meta": self.skipped_no_meta,
            "skipped_no_peaks": self.skipped_no_peaks,
            "errors": self.errors,
            "elapsed_sec": round(self.elapsed, 1),
            "entries_per_minute": round(self.per_minute, 1),
        }


def _compute_one(task: tuple[str, int]) -> tuple[int, PhaseMeta | None, object, object, str]:
    """Compute one phase in a worker process.

    Kept at module level so :class:`~concurrent.futures.ProcessPoolExecutor` can
    pickle it.  Returns ``(cod_id, meta, d_arr, i_arr, reason)`` where *reason* is
    ``""`` on success and a short machine-readable marker otherwise.

    Peak computation is by far the most expensive step (a large-cell phase can
    take seconds), so it is what gets parallelised; the SQLite writer stays in the
    parent process because a database handle cannot be shared.
    """
    cif_path_str, cod_id = task
    cif_path = Path(cif_path_str)
    if not cif_path.is_file():
        return cod_id, None, None, None, "no_cif"
    meta = _parse_meta(cif_path, cod_id, _ID_OFFSET[0])
    if meta is None:
        return cod_id, None, None, None, "no_meta"
    peaks = compute_pattern(
        cif_path, dmin=_DMIN[0], max_peaks=_MAX_PEAKS[0], normalize=True
    )
    if peaks is None:
        return cod_id, None, None, None, "no_peaks"
    return cod_id, meta, peaks[0], peaks[1], ""


#: Worker-process globals: the parameter values `_compute_one` needs.  Set by
#: :func:`_init_worker` in every child, so they do not travel with each task.
_DMIN: list[float] = [DEFAULT_DMIN]
_MAX_PEAKS: list[int] = [DEFAULT_MAX_PEAKS]
#: Numeric shift that separates the archives in the entry-id space: COD keeps
#: its own ids, PCOD is shifted by +20 000 000 and TCOD by +100 000 000 (the
#: 3COD convention).  Without it PCOD ``1000000`` and COD ``1000000`` would
#: share the entry id ``96-100-0000`` and one library would shadow the other.
_ID_OFFSET: list[int] = [0]


def _init_worker(dmin: float, max_peaks: int, id_offset: int = 0) -> None:
    """Populate the worker-process globals (one call per child)."""
    _DMIN[0] = dmin
    _MAX_PEAKS[0] = max_peaks
    _ID_OFFSET[0] = int(id_offset)


def build(
    cif_root: Path,
    out: Path,
    *,
    library_kind: str,
    build_tool: str = "tools.build_cod_xindb",
    data_date: str = "",
    dmin: float = DEFAULT_DMIN,
    max_peaks: int = DEFAULT_MAX_PEAKS,
    start: int = 0,
    limit: int | None = None,
    resume: bool = False,
    tier: int = 1,
    log_every: int = 500,
    workers: int = 1,
    progress_cb=None,
    id_offset: int = 0,
) -> BuildStats:
    """Build a v3 library from the CIF archive at *cif_root*.

    Args:
        cif_root: CIF archive root (the fan-out directory).
        out: Output ``.xinDB`` path.
        library_kind: e.g. ``"COD-Inorg"`` / ``"PCOD-Inorg"`` — recorded in ``meta``.
        dmin: Minimum d-spacing (Å) kept in the pattern.
        max_peaks: Strongest reflections kept per phase.
        start: Lowest archive id to process (pre-offset, i.e. as found on disk).
        id_offset: Numeric shift applied when forming the stored entry id, so
            the three archives cannot collide (COD 0, PCOD +20 000 000,
            TCOD +100 000 000).  It only affects the id, never the file walk
            and never the peak data.
        limit: Maximum number of CIFs to process (``None`` = all).
        resume: Continue after the highest COD id already written.
        tier: Quality tier stored for every phase (0 = experimental first).
        log_every: Log a progress line every N phases.
        workers: Number of peak-computation processes (1 = in-process).  The
            SQLite writer always stays in the parent, and results are consumed in
            archive order, so the output file is identical for any worker count.
        progress_cb: Optional ``(done, total)`` callback.

    Returns:
        A :class:`BuildStats`.
    """
    t0 = time.perf_counter()

    out.parent.mkdir(parents=True, exist_ok=True)

    # Provenance travels through the writer, not through a separate
    # ``write_standard_meta`` call: ``writer.close()`` stamps the full meta set
    # (including ``element_bits``), so a second writer-side call would silently
    # overwrite whatever this function set.  That is exactly how every library
    # built before this fix ended up with ``library_kind = "unknown"``.
    writer = XinDBWriter(
        out,
        top_k=TOP_K,
        chunk_size=CHUNK_SIZE,
        tier=tier,
        cif_source=str(cif_root),
        library_kind=library_kind,
        build_tool=build_tool,
        data_date=data_date,
        resume=bool(resume),
    )

    # ``--resume`` continues an interrupted build instead of restarting it.
    # The writer has already wound the file back to a chunk boundary (dropping
    # any entries whose peak chunk never landed) and reloaded its counters, so
    # the only thing left is to point the CIF walk at the right phase.  Where
    # that is has to come from the writer, not from a pre-read of the file:
    #
    # * the trim usually *discards* the tail, and those phases were never
    #   written to a chunk, so continuing after the last surviving phase would
    #   drop them for good -- the walk has to restart *at* the first discarded
    #   phase;
    # * when nothing was discarded the walk continues after the last phase.
    if resume:
        if writer.resumed_from:
            first_redo = writer.resume_from_entry_id
            if first_redo is not None:
                resume_cod_id = _cod_id_of(first_redo)
                if resume_cod_id is None:
                    raise RuntimeError(
                        f"cannot parse the COD id of the first discarded phase "
                        f"{first_redo!r}; refusing to guess where to resume"
                    )
                # ``resume_from_entry_id`` is a *stored* id (offset applied);
                # the walk works in archive ids, so undo the shift here.
                start = max(int(start or 0), resume_cod_id - int(id_offset))
                logger.warning(
                    "resuming %s: %d phases already stored, restarting the walk at "
                    "the first discarded phase (archive id %d)",
                    out.name,
                    writer.resumed_from,
                    resume_cod_id - int(id_offset),
                )
            else:
                last_id = writer.conn.execute(
                    "SELECT entry_id FROM entries ORDER BY rowid DESC LIMIT 1"
                ).fetchone()
                next_cod_id = _cod_id_of(last_id[0]) if last_id else None
                if next_cod_id is None:
                    raise RuntimeError(
                        f"cannot parse the COD id of the last stored phase "
                        f"{last_id[0]!r}; refusing to guess where to resume"
                    )
                start = max(int(start or 0), next_cod_id + 1 - int(id_offset))
                logger.warning(
                    "resuming %s: %d phases already stored, walking the archive from "
                    "COD id %d",
                    out.name,
                    writer.resumed_from,
                    start,
                )
        else:
            logger.warning("--resume given but %s holds no phases; starting fresh", out.name)

    stats = BuildStats()
    stats.resumed_from = int(writer.resumed_from)
    last_done = start
    total = 0
    workers = max(1, int(workers))

    tasks = (
        (str(cif_path), cod_id)
        for cod_id, cif_path in iter_cod_cifs(cif_root, start=start, limit=limit)
    )

    def _consume(results) -> None:
        nonlocal last_done
        for cod_id, meta, d_arr, i_arr, reason in results:
            last_done = cod_id
            if reason == "no_cif":
                stats.skipped_no_cif += 1
                continue
            if reason == "no_meta":
                stats.skipped_no_meta += 1
                continue
            if reason == "no_peaks":
                stats.skipped_no_peaks += 1
                continue
            assert meta is not None
            writer.add(meta, d_arr, i_arr)
            stats.entries += 1
            stats.entries_total = stats.resumed_from + stats.entries
            stats.peaks += int(d_arr.size)  # type: ignore[union-attr]

            if progress_cb is not None and stats.entries % log_every == 0:
                progress_cb(stats.entries, total)
            if stats.entries % log_every == 0:
                rate = stats.entries / max(time.perf_counter() - t0, 1e-9) * 60.0
                logger.info(
                    "built %d phases this run (%d in file, %.0f/min) last_cod_id=%d",
                    stats.entries, stats.entries_total, rate, last_done,
                )

    try:
        if workers > 1:
            from concurrent.futures import ProcessPoolExecutor

            with ProcessPoolExecutor(
                max_workers=workers,
                initializer=_init_worker,
                initargs=(dmin, max_peaks, id_offset),
            ) as pool:
                # imap keeps archive order, so the file does not depend on how
                # the workers happen to interleave.
                _consume(pool.map(_compute_one, tasks, chunksize=8))
        else:
            _init_worker(dmin, max_peaks, id_offset)
            _consume(map(_compute_one, tasks))

        # ``writer.close()`` writes the meta table inside the same transaction as
        # the final chunk, so a file whose meta table exists is complete and no
        # separate flush is needed here.
    finally:
        closed = writer.close()
    # ``stats.entries`` counts what *this* run added; the file itself holds the
    # inherited phases too, and only the writer knows that total.
    stats.entries_total = int(closed.get("entries", stats.entries))

    stats.elapsed = time.perf_counter() - t0
    stats.per_minute = (
        stats.entries / stats.elapsed * 60.0 if stats.elapsed > 0 else 0.0
    )
    logger.info("build finished: %s", stats.as_dict())
    return stats


def _cod_id_of(entry_id: str) -> int | None:
    """Recover the numeric COD id from a stored ``96-XXX-XXXX`` entry id.

    The display id is built by ``entry_id_for_cod`` as ``"96-" + zfill(7)[:3] +
    "-" + zfill(7)[3:]``, so stripping the source prefix and the dash returns the
    original id.  ``None`` means the string is not in that shape and the caller
    must not guess.
    """
    parts = str(entry_id).split("-")
    if len(parts) < 2:
        return None
    digits = "".join(parts[1:])
    return int(digits) if digits.isdigit() and digits else None


def _parse_meta(cif_path: Path, cod_id: int, id_offset: int = 0) -> PhaseMeta | None:
    """Parse a CIF's metadata, tolerating malformed files.

    The ``entry_id`` must be the canonical display id (``96-XXX-XXXX``), the same
    form the v2 library and the GUI use, so it is derived with
    :func:`~xinmatch.core.db.cif_scan.entry_id_for_cod` rather than from the bare
    numeric id.  *id_offset* is added first so PCOD/TCOD phases land in their own
    id range instead of colliding with COD.
    """
    from xinmatch.core.db.cif_meta import parse_meta
    from xinmatch.core.db.cif_scan import entry_id_for_cod

    try:
        return parse_meta(cif_path, entry_id=entry_id_for_cod(cod_id + int(id_offset)))
    except Exception as exc:  # noqa: BLE001 - one bad CIF must not kill the build
        logger.debug("metadata parse failed for %s: %s", cif_path, exc)
        return None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--out", type=Path, default=None,
        help="output .xinDB path (may live in any folder); defaults to "
             "$XINMATCH_XINDB_DIR/<kind>-<data-date>.xinDB when that variable is "
             "set, otherwise <repo>/data/<kind>-<data-date>.xinDB",
    )
    parser.add_argument(
        "--kind", choices=sorted(LIBRARY_KINDS), default="full",
        help="which library to build",
    )
    parser.add_argument(
        "--data-date", default=None,
        help="COD snapshot date (YYYYMMDD) recorded in meta.data_date and used "
             "in the default output name; defaults to today",
    )
    parser.add_argument(
        "--cif-root", type=Path, default=None,
        help="CIF archive root (defaults to the archive for --kind)",
    )
    parser.add_argument("--dmin", type=float, default=DEFAULT_DMIN)
    parser.add_argument("--max-peaks", type=int, default=DEFAULT_MAX_PEAKS)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--tier", type=int, default=1)
    parser.add_argument("--log-every", type=int, default=500)
    parser.add_argument(
        "--workers", type=int, default=1,
        help="peak-computation processes (1 = in-process); the writer stays single "
             "and results are consumed in archive order, so the output is identical",
    )
    parser.add_argument(
        "--id-offset", type=int, default=0,
        help="numeric shift applied to the stored entry id so the archives do "
             "not collide (COD 0, PCOD +20000000, TCOD +100000000); it does not "
             "change which files are read",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    default_root, kind_label, _desc = LIBRARY_KINDS[args.kind]
    cif_root = args.cif_root or default_root
    if not cif_root.is_dir():
        print(f"[ERR] CIF root not found: {cif_root}")
        return 2

    # The COD snapshot date is part of the identity of the library: two builds
    # of the same kind from different snapshots must not be mistaken for each
    # other, so it lands in both the file name and ``meta.data_date``.
    data_date = (args.data_date or time.strftime("%Y%m%d")).strip()
    out = args.out or (default_out_dir() / f"{kind_label}-{data_date}{XINDB_SUFFIX}")

    print(f"kind:     {kind_label}")
    print(f"data:     {data_date}")
    print(f"cif root: {cif_root}")
    print(f"id offset:{args.id_offset:+,}")
    print(f"out:      {out}")
    print()

    stats = build(
        cif_root,
        out,
        library_kind=kind_label,
        data_date=data_date,
        dmin=args.dmin,
        max_peaks=args.max_peaks,
        start=args.start,
        limit=args.limit,
        resume=args.resume,
        tier=args.tier,
        log_every=args.log_every,
        workers=args.workers,
        id_offset=args.id_offset,
    )

    size_mb = out.stat().st_size / 1e6 if out.exists() else 0.0
    print("\nbuild summary:")
    for key, value in stats.as_dict().items():
        print(f"  {key}: {value}")
    print(f"  size_mb: {size_mb:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
