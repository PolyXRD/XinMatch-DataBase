"""Prove that ``--resume`` reproduces a fresh build, table for table.

``.xinDB`` is written by two independently committed streams: ``entries`` every
``_ENTRY_BATCH`` (2000) phases, ``peak_chunks`` every ``chunk_size`` (4096)
phases.  A killed process therefore leaves a file whose two streams do not line
up, and appending to such a file is only lossless if the writer winds both of
them back to the last boundary they agree on.  Reasoning about that on paper is
not enough: the streams feed an integer ``rowid`` / ``peak_off`` chain, so a
mistake shows up as silently shifted reflection offsets rather than an error.

This script settles it empirically.  It builds a **reference** library from *N*
phases in one uninterrupted run, then for each residue state:

1. writes the same phases in a child process that dies with ``os._exit``
   before ``close()`` -- no rollback, no flush, exactly like a kill;
2. reads back how far the file actually got (``entries`` vs ``peak_chunks``);
3. resumes it and writes the remaining phases;
4. compares **every table** against the reference, byte for byte.

Two residue states are distinguishable and both are exercised:

``entries ahead``
    the peak payload of the tail never landed.  The common case: ``entries``
    flushes twice as often as a chunk.
``chunks ahead``
    a chunk landed whose own tail entries were still buffered.  This happens in
    the window right after a chunk commit, which is not rare at all -- with
    ``chunk_size=4096`` and ``_ENTRY_BATCH=2000`` it covers roughly 44% of the
    phase counts a kill can land on.

A green run means the resumed file is indistinguishable from the one a
35-hour uninterrupted build would have produced.  ``rowid``, ``peak_off``,
chunk ids, chunk payloads and the element inverted index all have to match.

The phase data is synthetic by default: ``add()`` stores whatever it is handed,
so the physics is irrelevant to the question being asked, and a deterministic
generator makes the run seconds instead of half an hour.  Pass ``--real`` to
repeat it on genuine COD phases (slow, same conclusion).

Usage
-----
::

    python -m tools._verify_resume_equivalence
    python -m tools._verify_resume_equivalence --phases 12000 --keep
    python -m tools._verify_resume_equivalence --real --phases 10000

Exit code 0 means every case matched; 1 means at least one table differed.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import os
import pickle
import sqlite3
import subprocess
import sys
import types
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from xinmatch.core.db.schema_v3 import CHUNK_SIZE, TOP_K  # noqa: E402
from xinmatch.core.db.writer_v3 import XinDBWriter  # noqa: E402

logger = logging.getLogger("verify-resume")

WORK = REPO / "data" / "_resume_check"
CACHE = WORK / "phases.pkl"

#: Writer arguments shared by the reference and every resumed run.  ``build_time``
#: is deliberately left out of the comparison: it records *when* the file was
#: made and is expected to differ.
WRITER_KW = dict(
    top_k=TOP_K,
    chunk_size=CHUNK_SIZE,
    tier=1,
    cif_source="resume-equivalence-check",
    library_kind="full",
    build_tool="tools._verify_resume_equivalence",
    data_date="20261020",
)

#: ``(label, phases written before the child dies, residue we expect)``.  The
#: expected residue is checked against what actually happened, so a case that
#: silently stops exercising its state fails loudly instead of passing for the
#: wrong reason.
CASES = [
    ("A entries-ahead", 15_000, "entries"),
    ("B chunks-ahead", 13_000, "chunks"),
]

ELEMENTS = ["H", "C", "N", "O", "F", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
            "K", "Ca", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
            "Ga", "Ge", "As", "Se", "Br", "Sr", "Y", "Zr", "Nb", "Mo", "Ag",
            "Cd", "In", "Sn", "Sb", "Te", "I", "Ba", "La", "Ce", "W", "Pb"]


# ── phase data ───────────────────────────────────────────────────
def synth_phases(n: int):
    """Deterministic fake phases with varied peak counts and compositions.

    The variety matters: the offset chain is only stressed when phases carry
    different numbers of reflections, and the element inverted index only spans
    several symbols when the composition varies.  Both must come out identical
    from the reference and the resumed run.
    """
    import numpy as np

    rng = np.random.default_rng(20261020)
    out = []
    for k in range(n):
        n_peaks = int(rng.integers(8, 151))
        d = np.sort(rng.uniform(0.8, 12.0, n_peaks))[::-1]
        i = rng.uniform(1.0, 1000.0, n_peaks)
        n_el = int(rng.integers(1, 5))
        elements = list(rng.choice(ELEMENTS, size=n_el, replace=False))
        out.append(
            (
                types.SimpleNamespace(
                    entry_id=f"96-{k // 10000:03d}-{k % 10000:04d}",
                    formula="-".join(elements),
                    space_group=rng.choice(["P 21/c", "F m -3 m", "I 4/m m m", "R -3 c"]),
                    crystal_system=rng.choice(["monoclinic", "cubic", "tetragonal", "trigonal"]),
                    mineral_name=None,
                    elements=[str(e) for e in elements],
                    a=float(rng.uniform(3, 20)),
                    b=float(rng.uniform(3, 20)),
                    c=float(rng.uniform(3, 20)),
                    alpha=90.0,
                    beta=float(rng.uniform(80, 110)),
                    gamma=90.0,
                    volume=float(rng.uniform(100, 5000)),
                    density=float(rng.uniform(2, 9)),
                    year=int(rng.integers(1950, 2026)),
                    author="resume-check",
                    journal="J. Test",
                    quality="calculated",
                    tier=1,
                ),
                d,
                i,
            )
        )
    return out


def real_phases(n: int):
    """Genuine COD phases, cached on disk (slow: one powder pattern each)."""
    cached = None
    if CACHE.is_file():
        with CACHE.open("rb") as fh:
            cached = pickle.load(fh)
    if cached is not None and len(cached) >= n:
        logger.info("using %d cached real phases", len(cached))
        return cached[:n]

    from concurrent.futures import ProcessPoolExecutor

    from tools.build_cod_xindb import _compute_one, _init_worker
    from xinmatch.core.db.cif_scan import DEFAULT_COD_ROOT, iter_cod_cifs

    WORK.mkdir(parents=True, exist_ok=True)
    tasks = [(str(p), cid) for cid, p in iter_cod_cifs(DEFAULT_COD_ROOT, limit=n)]
    workers = max(1, min(14, (os.cpu_count() or 4)))
    logger.info("computing %d real phases with %d workers (one-off, then cached)",
                len(tasks), workers)
    out = []
    with ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(0.7703, 150)
    ) as pool:
        for cod_id, meta, d_arr, i_arr, reason in pool.map(_compute_one, tasks, chunksize=8):
            if reason:
                continue
            out.append((meta, d_arr, i_arr))
    with CACHE.open("wb") as fh:
        pickle.dump(out, fh, protocol=pickle.HIGHEST_PROTOCOL)
    if len(out) < n:
        logger.warning("archive yielded only %d of %d requested phases", len(out), n)
    return out[:n]


# ── writers ──────────────────────────────────────────────────────
def build_reference(path: Path, phases) -> dict[str, int]:
    """Write every phase in one uninterrupted run."""
    writer = XinDBWriter(path, resume=False, **WRITER_KW)
    for meta, d_arr, i_arr in phases:
        writer.add(meta, d_arr, i_arr)
    return writer.close()


def child_die(db_path: Path, n_add: int, real: bool) -> None:
    """Write *n_add* phases then vanish without finalising the database."""
    phases = real_phases(n_add) if real else synth_phases(n_add)
    writer = XinDBWriter(db_path, resume=False, **WRITER_KW)
    for meta, d_arr, i_arr in phases:
        writer.add(meta, d_arr, i_arr)
    # A hard exit skips close() *and* every pending buffered write, which is
    # exactly what an external kill does: committed pages are on disk, nothing
    # else is.
    sys.stdout.flush()
    os._exit(0)


# ── comparison ───────────────────────────────────────────────────
def _hash_rows(rows) -> str:
    h = hashlib.blake2b(digest_size=16)
    for row in rows:
        h.update(repr(row).encode("utf-8", "surrogatepass"))
        h.update(b"\x1f")
    return h.hexdigest()


def describe(path: Path) -> dict[str, object]:
    """Fingerprint every table that carries phase data."""
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    try:
        out: dict[str, object] = {}
        out["entries"] = conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
        out["peaks"] = conn.execute("SELECT COALESCE(SUM(peak_n),0) FROM entries").fetchone()[0]
        out["chunks"] = conn.execute("SELECT COUNT(*) FROM peak_chunks").fetchone()[0]
        out["elements"] = conn.execute("SELECT COUNT(*) FROM element_inv").fetchone()[0]
        # rowid + every column: catches a shifted index or a wrong offset even
        # when the counts happen to agree.
        out["entries_hash"] = _hash_rows(
            conn.execute("SELECT rowid,* FROM entries ORDER BY rowid")
        )
        # Raw payload bytes: catches a stale row left in the reused top_k buffer,
        # a mis-sliced chunk or a wrong chunk id.
        out["chunks_hash"] = _hash_rows(
            conn.execute(
                "SELECT chunk_id,n,top_d,top_i,all_d,all_i FROM peak_chunks ORDER BY chunk_id"
            )
        )
        out["elements_hash"] = _hash_rows(
            conn.execute("SELECT element,eids FROM element_inv ORDER BY element")
        )
        out["meta"] = dict(
            conn.execute("SELECT key,value FROM meta WHERE key != 'build_time' ORDER BY key")
        )
        return out
    finally:
        conn.close()


def compare(ref_path: Path, got_path: Path) -> list[str]:
    """Return the differences between two libraries (empty when they agree)."""
    ref, got = describe(ref_path), describe(got_path)
    diffs: list[str] = []
    for key in ("entries", "peaks", "chunks", "elements",
                "entries_hash", "chunks_hash", "elements_hash"):
        if ref[key] != got[key]:
            diffs.append(f"{key}: reference={ref[key]!r} resumed={got[key]!r}")
    if not ref["meta"]:
        diffs.append("reference has no meta table (it did not close?)")
    if ref["meta"] != got["meta"]:
        bad = {
            k: (ref["meta"].get(k), got["meta"].get(k))
            for k in set(ref["meta"]) | set(got["meta"])
            if ref["meta"].get(k) != got["meta"].get(k)
        }
        diffs.append(f"meta differs (reference, resumed): {bad}")
    return diffs


def residue(path: Path) -> tuple[str, int, int]:
    """Report which of the two streams is ahead, and by how much."""
    conn = sqlite3.connect(f"file:{path}?mode=ro&immutable=1", uri=True)
    try:
        n = conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
        covered = conn.execute("SELECT COALESCE(SUM(n),0) FROM peak_chunks").fetchone()[0]
    finally:
        conn.close()
    if n > covered:
        return "entries", n, covered
    if covered > n:
        return "chunks", n, covered
    return "aligned", n, covered


# ── driver ───────────────────────────────────────────────────────
def run(phases_n: int, keep: bool, real: bool) -> int:
    WORK.mkdir(parents=True, exist_ok=True)
    phases = real_phases(phases_n) if real else synth_phases(phases_n)
    n = len(phases)
    if n < CHUNK_SIZE * 3:
        raise SystemExit(f"need at least {CHUNK_SIZE * 3} phases, got {n}")

    ref_path = WORK / "reference.xinDB"
    ref_path.unlink(missing_ok=True)
    logger.info("building reference library (%d %s phases, one uninterrupted run)",
                n, "real" if real else "synthetic")
    ref_stats = build_reference(ref_path, phases)
    ref = describe(ref_path)
    logger.info("reference: %d phases / %d reflections / %d chunks / %d elements / %.1f MB",
                ref["entries"], ref["peaks"], ref["chunks"], ref["elements"],
                ref_path.stat().st_size / 1e6)

    failures = 0
    rows: list[tuple[str, str, int, str]] = []
    for label, cut, expect in CASES:
        cut = min(cut, n)
        db = WORK / f"interrupted-{label.split()[0]}.xinDB"
        db.unlink(missing_ok=True)

        proc = subprocess.run(
            [sys.executable, __file__, "--child-die", str(db), str(cut),
             "--real" if real else "--synth"],
            cwd=str(REPO),
            stdout=subprocess.DEVNULL,
        )
        if proc.returncode != 0 or not db.exists():
            logger.error("[%s] child died abnormally (rc=%d)", label, proc.returncode)
            failures += 1
            continue

        state, n_entries, n_covered = residue(db)
        logger.info("[%s] killed after %d phases -> residue '%s': %d entries vs %d in chunks",
                    label, cut, state, n_entries, n_covered)
        if state != expect:
            logger.error("[%s] this cut no longer produces the '%s' residue (got '%s'); "
                         "the case is not testing what it claims", label, expect, state)
            failures += 1

        writer = XinDBWriter(db, resume=True, **WRITER_KW)
        resumed_from = writer.resumed_from
        redo_id = writer.resume_from_entry_id
        logger.info("[%s] wound back to %d phases; recompute starts at %s",
                    label, resumed_from, redo_id or "<after the last stored phase>")
        if resumed_from % CHUNK_SIZE:
            logger.error("[%s] resume boundary %d is not chunk aligned", label, resumed_from)
            failures += 1
            writer.conn.close()
            continue
        for meta, d_arr, i_arr in phases[resumed_from:]:
            writer.add(meta, d_arr, i_arr)
        stats = writer.close()
        logger.info("[%s] resumed run wrote %d phases -> %d in file, %d reflections",
                    label, stats["entries"] - resumed_from, stats["entries"], stats["peaks"])

        # ``mask_truncated`` is a statistic, not a stored column: a resumed run
        # that quietly recounts only its own phases would under-report it.
        if stats["mask_truncated"] != ref_stats["mask_truncated"]:
            logger.error(
                "[%s] mask_truncated differs: reference=%d resumed=%d",
                label, ref_stats["mask_truncated"], stats["mask_truncated"],
            )
            failures += 1

        diffs = compare(ref_path, db)
        if diffs:
            failures += 1
            logger.error("[%s] FAILED -- %d difference(s)", label, len(diffs))
            for d in diffs:
                logger.error("    %s", d)
            verdict = "FAIL"
        else:
            logger.info("[%s] OK -- identical to the reference in every table", label)
            verdict = "OK"
        rows.append((label, f"{state} ({n_entries}/{n_covered})", resumed_from, verdict))
        if not keep and verdict == "OK":
            db.unlink(missing_ok=True)

    print()
    print(f"phases / case : {n:,}")
    print(f"reference     : {ref['entries']:,} phases, {ref['peaks']:,} reflections, "
          f"{ref['chunks']:,} chunks, {ref['elements']} elements")
    print()
    print(f"{'case':<20}{'residue (entries/chunks)':<30}{'resumed at':>12}  verdict")
    print("-" * 76)
    for label, state, resumed_from, verdict in rows:
        print(f"{label:<20}{state:<30}{resumed_from:>12,}  {verdict}")
    print("-" * 76)
    print(f"failures: {failures}")
    print("RESULT  : " + ("PASS -- resume is lossless" if failures == 0 else "FAIL"))
    return 0 if failures == 0 else 1


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--phases", type=int, default=16_000, help="phases per case (default 16000)")
    ap.add_argument("--real", action="store_true", help="use genuine COD phases instead of synthetic")
    ap.add_argument("--synth", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--keep", action="store_true", help="keep the interrupted databases")
    ap.add_argument("--child-die", nargs=2, metavar=("DB", "N"), help=argparse.SUPPRESS)
    args = ap.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s %(levelname)-7s %(message)s", datefmt="%H:%M:%S"
    )
    if args.child_die:
        child_die(Path(args.child_die[0]), int(args.child_die[1]), bool(args.real))
        return 0
    return run(args.phases, args.keep, bool(args.real))


if __name__ == "__main__":
    raise SystemExit(main())
