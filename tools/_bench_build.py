"""Throughput benchmark for the single-parse build hot path (task-5 diagnostic).

The library builder used to compute every phase in two steps::

    meta  = parse_meta(cif_path, entry_id=...)
    peaks = compute_pattern(cif_path, dmin=..., max_peaks=...)

which parsed each CIF three times (``read_small_structure`` + ``cif.read`` in
``parse_meta``, then another ``read_small_structure`` in ``compute_pattern``).
``build_pipeline.compute_phase`` does the same work from a single
``gemmi.cif.read``.  This script measures that change on a deterministic sample
of the **real** COD archive (read-only) and proves the two paths agree bit for
bit:

* sequential (single-process) per-file time, measured **paired** (old then new
  on the same file, repeated) so machine load cancels out;
* the broken-down parse cost, to show how much of a phase's time is parsing;
* ``ProcessPoolExecutor`` throughput for the ``--workers`` / ``--chunksize``
  grids (the builder parallelises exactly like this);
* parity fingerprint per file: the full ``PhaseMeta`` (pickled ``asdict``) and
  the raw ``(d, i)`` array bytes must be identical.

Nothing outside the output JSON (default ``tmp/bench/``) is written; no library
and no CIF is modified.

Usage::

    python -m tools._bench_build                    # 800 CIFs, default grid
    python -m tools._bench_build --limit 500
    python -m tools._bench_build --cif-root D:/Project/XRD/sync-COD/inorg_cif/COD

The JSON result is written to ``tmp/bench/bench_build.json`` (not tracked by git).
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import pickle
import platform
import random
import statistics
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import gemmi  # noqa: E402
import numpy  # noqa: E402

from xinmatch.core.db.build_pipeline import compute_phase  # noqa: E402
from xinmatch.core.db.cif_meta import parse_meta  # noqa: E402
from xinmatch.core.db.cif_scan import (  # noqa: E402
    DEFAULT_INORG_COD_ROOT,
    entry_id_for_cod,
    iter_cod_cifs,
)
from xinmatch.core.db.peak_compute import (  # noqa: E402
    DEFAULT_DMIN,
    DEFAULT_MAX_PEAKS,
    compute_pattern,
)

#: One benchmark unit: CIF path as ``str`` (picklable), entry id and engine args.
Task = tuple[str, str, float, int]

#: Default result path (``tmp/`` is scratch space, never committed).
DEFAULT_OUT = PROJECT_ROOT / "tmp" / "bench" / "bench_build.json"

#: Seed of the deterministic CIF sample, so two runs benchmark the same files.
DEFAULT_SEED = 20261024

#: A phase slower than this (ms) is a "large-cell / high-symmetry" outlier.
SLOW_PHASE_MS = 100.0


# ---------------------------------------------------------------------------
# The two implementations under test
# ---------------------------------------------------------------------------
def _old_two_step(task: Task) -> tuple[str, Any, Any]:
    """The shipped two-step sequence: ``parse_meta`` then ``compute_pattern``."""
    path_str, entry_id, dmin, max_peaks = task
    path = Path(path_str)
    if not path.is_file():
        return "no_cif", None, None
    meta = parse_meta(path, entry_id=entry_id)
    if meta is None:
        return "no_meta", None, None
    peaks = compute_pattern(path, dmin=dmin, max_peaks=max_peaks, normalize=True)
    if peaks is None:
        return "no_peaks", None, None
    return "", meta, (peaks[0], peaks[1])


def _new_single_parse(task: Task) -> tuple[str, Any, Any]:
    """The new single-parse path: ``build_pipeline.compute_phase``."""
    path_str, entry_id, dmin, max_peaks = task
    meta, peaks, reason = compute_phase(
        Path(path_str),
        entry_id=entry_id,
        dmin=dmin,
        max_peaks=max_peaks,
        normalize=True,
    )
    return reason, meta, peaks


def _fingerprint(
    result: tuple[str, Any, Any],
) -> tuple[str, bytes | None, tuple[bytes, bytes] | None]:
    """Bit-exact fingerprint of one phase result (metadata + peak payload).

    Metadata is pickled instead of compared field by field so that ``-0.0`` vs
    ``0.0`` and NaN payloads cannot slip through a plain ``==``.  Peak arrays are
    compared as raw bytes.
    """
    reason, meta, peaks = result
    meta_fp = None if meta is None else pickle.dumps(dataclasses.asdict(meta), protocol=4)
    peaks_fp = None if peaks is None else (peaks[0].tobytes(), peaks[1].tobytes())
    return reason, meta_fp, peaks_fp


def _worker_old(task: Task) -> tuple[str, int]:
    """Parallel worker for the old path (returns a small summary only)."""
    reason, _meta, peaks = _old_two_step(task)
    return reason, 0 if peaks is None else int(peaks[0].size)


def _worker_new(task: Task) -> tuple[str, int]:
    """Parallel worker for the new path (returns a small summary only)."""
    reason, _meta, peaks = _new_single_parse(task)
    return reason, 0 if peaks is None else int(peaks[0].size)


# ---------------------------------------------------------------------------
# Statistics helpers
# ---------------------------------------------------------------------------
def _percentiles(samples: list[float]) -> dict[str, float]:
    """Summarise per-file timings (ms)."""
    ordered = sorted(samples)
    return {
        "min": ordered[0],
        "p50": statistics.median(ordered),
        "p90": ordered[min(len(ordered) - 1, int(0.90 * len(ordered)))],
        "p99": ordered[min(len(ordered) - 1, int(0.99 * len(ordered)))],
        "max": ordered[-1],
        "mean": statistics.fmean(ordered),
        "sum_ms": sum(ordered),
    }


def _speedup(old_value: float, new_value: float) -> float:
    """Times faster of *new_value* (guard against a zero divisor)."""
    return old_value / new_value if new_value > 0.0 else float("inf")


# ---------------------------------------------------------------------------
# Sequential (paired) benchmark
# ---------------------------------------------------------------------------
def bench_sequential(tasks: list[Task], *, repeats: int, check_parity: bool) -> dict[str, Any]:
    """Compare both implementations in one process with paired per-file timing.

    Old and new are run back-to-back on the same file, and every repetition
    keeps the per-file minimum.  That removes most machine-load noise without
    hiding the real per-phase distribution.
    """
    count = len(tasks)
    old_ms = [float("inf")] * count
    new_ms = [float("inf")] * count
    checked = 0
    mismatches = 0
    first_mismatch: dict[str, Any] | None = None
    old_reasons: dict[str, int] = {}
    new_reasons: dict[str, int] = {}

    for repeat_index in range(repeats):
        for index, task in enumerate(tasks):
            start = time.perf_counter()
            old_result = _old_two_step(task)
            middle = time.perf_counter()
            new_result = _new_single_parse(task)
            end = time.perf_counter()
            old_ms[index] = min(old_ms[index], (middle - start) * 1e3)
            new_ms[index] = min(new_ms[index], (end - middle) * 1e3)
            if not check_parity:
                continue
            checked += 1
            old_fp = _fingerprint(old_result)
            new_fp = _fingerprint(new_result)
            if repeat_index == repeats - 1:
                old_reasons[old_fp[0]] = old_reasons.get(old_fp[0], 0) + 1
                new_reasons[new_fp[0]] = new_reasons.get(new_fp[0], 0) + 1
            if old_fp != new_fp:
                mismatches += 1
                if first_mismatch is None:
                    first_mismatch = {
                        "cif": task[0],
                        "entry_id": task[1],
                        "reason_old": old_fp[0],
                        "reason_new": new_fp[0],
                        "meta_equal": old_fp[1] == new_fp[1],
                        "peaks_equal": old_fp[2] == new_fp[2],
                    }
        print(f"  sequential repeat {repeat_index + 1}/{repeats} done", flush=True)

    # ---- parse-cost breakdown on the same files (second pass) --------------
    old_parse_ms: list[float] = []
    new_parse_ms: list[float] = []
    for task in tasks:
        path = Path(task[0])
        entry_id = task[1]
        # old: read_small_structure + cif.read (parse_meta) + read_small_structure
        # (compute_pattern) -- i.e. every parse the old hot path performed.
        start = time.perf_counter()
        parse_meta(path, entry_id=entry_id)
        gemmi.read_small_structure(str(path))
        old_parse_ms.append((time.perf_counter() - start) * 1e3)
        # new: one cif.read + one make_small_structure_from_block.
        start = time.perf_counter()
        block = gemmi.cif.read(str(path)).sole_block()
        gemmi.make_small_structure_from_block(block)
        new_parse_ms.append((time.perf_counter() - start) * 1e3)

    old_total = sum(old_ms)
    new_total = sum(new_ms)
    old_parse_total = sum(old_parse_ms)
    new_parse_total = sum(new_parse_ms)

    ratios = [_speedup(old, new) for old, new in zip(old_ms, new_ms, strict=True)]
    ratios.sort()
    slow_old = [index for index in range(count) if old_ms[index] >= SLOW_PHASE_MS]
    fast_old = [index for index in range(count) if old_ms[index] < SLOW_PHASE_MS]
    fast_old_total = sum(old_ms[index] for index in fast_old)
    fast_new_total = sum(new_ms[index] for index in fast_old)

    return {
        "files": count,
        "repeats": repeats,
        "time_weighted": {
            "old_two_step": {
                "sum_seconds": old_total / 1e3,
                "files_per_s": count / (old_total / 1e3),
                "ms_per_file": old_total / count,
            },
            "compute_phase": {
                "sum_seconds": new_total / 1e3,
                "files_per_s": count / (new_total / 1e3),
                "ms_per_file": new_total / count,
            },
            "speedup_x": _speedup(old_total, new_total),
            "speedup_pct": (_speedup(old_total, new_total) - 1.0) * 100.0,
            "meets_5pct": _speedup(old_total, new_total) >= 1.05,
        },
        "per_file": {
            "old_two_step": _percentiles(old_ms),
            "compute_phase": _percentiles(new_ms),
            "speedup": {
                "p10": ratios[int(0.10 * count)],
                "p50": ratios[int(0.50 * count)],
                "p90": ratios[int(0.90 * count)],
                "mean": statistics.fmean(ratios),
            },
        },
        "fast_phases": {
            "threshold_ms": SLOW_PHASE_MS,
            "files": len(fast_old),
            "share_of_files_pct": len(fast_old) / count * 100.0,
            "old_seconds": fast_old_total / 1e3,
            "new_seconds": fast_new_total / 1e3,
            "speedup_pct": (_speedup(fast_old_total, fast_new_total) - 1.0) * 100.0,
        },
        "slow_phases": {
            "threshold_ms": SLOW_PHASE_MS,
            "files": len(slow_old),
            "share_of_files_pct": len(slow_old) / count * 100.0,
            "share_of_old_time_pct": (
                sum(old_ms[index] for index in slow_old) / old_total * 100.0
            ),
            "speedup_pct": (
                _speedup(
                    sum(old_ms[index] for index in slow_old),
                    sum(new_ms[index] for index in slow_old),
                )
                - 1.0
            )
            * 100.0,
        },
        "parse_cost": {
            "old_seconds": old_parse_total / 1e3,
            "new_seconds": new_parse_total / 1e3,
            "old_ms_per_file": old_parse_total / count,
            "new_ms_per_file": new_parse_total / count,
            "speedup_x": _speedup(old_parse_total, new_parse_total),
            "old_share_of_total_pct": old_parse_total / old_total * 100.0,
            "new_share_of_total_pct": new_parse_total / new_total * 100.0,
        },
        "parity": {
            "unique_phases": count,
            "checked": checked,
            "mismatches": mismatches,
            "first_mismatch": first_mismatch,
            "old_reasons": old_reasons,
            "new_reasons": new_reasons,
        },
    }


# ---------------------------------------------------------------------------
# Parallel throughput
# ---------------------------------------------------------------------------
def bench_parallel(
    worker,
    tasks: list[Task],
    *,
    workers: int,
    chunksize: int,
) -> dict[str, Any]:
    """Time one ``ProcessPoolExecutor`` configuration (startup included)."""
    start = time.perf_counter()
    if workers <= 1:
        results = [worker(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = list(pool.map(worker, tasks, chunksize=chunksize))
    elapsed = time.perf_counter() - start
    succeeded = sum(1 for reason, _n in results if reason == "")
    return {
        "workers": workers,
        "chunksize": chunksize,
        "seconds": elapsed,
        "files_per_s": len(tasks) / elapsed,
        "ok": succeeded,
        "skipped": len(tasks) - succeeded,
    }


# ---------------------------------------------------------------------------
# Sampling / environment / CLI
# ---------------------------------------------------------------------------
def sample_entries(
    cif_root: Path,
    *,
    limit: int,
    seed: int,
) -> tuple[int, list[tuple[int, Path]]]:
    """Deterministically sample *limit* ``(cod_id, path)`` entries from the archive."""
    all_entries = list(iter_cod_cifs(cif_root))
    total = len(all_entries)
    if limit >= total:
        return total, all_entries
    rng = random.Random(seed)
    selected = rng.sample(all_entries, limit)
    selected.sort(key=lambda item: item[0])
    return total, selected


def build_tasks(entries: list[tuple[int, Path]], *, dmin: float, max_peaks: int) -> list[Task]:
    """Turn archive entries into benchmark tasks with the builder's entry ids."""
    return [(str(path), entry_id_for_cod(cod_id), dmin, max_peaks) for cod_id, path in entries]


def machine_info() -> dict[str, Any]:
    """Record the environment the numbers were measured in."""
    return {
        "cpu_count": os.cpu_count(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "platform": platform.platform(),
        "python": sys.version.split()[0],
        "gemmi": gemmi.__version__,
        "numpy": numpy.__version__,
    }


def run(args: argparse.Namespace) -> int:
    """Execute the whole benchmark and return a process exit code."""
    cif_root = Path(args.cif_root)
    if not cif_root.is_dir():
        print(f"error: CIF archive root not found: {cif_root}", file=sys.stderr)
        return 2

    print(f"scanning {cif_root} ...")
    total, entries = sample_entries(cif_root, limit=args.limit, seed=args.seed)
    tasks = build_tasks(entries, dmin=args.dmin, max_peaks=args.max_peaks)
    print(
        f"sampled {len(tasks)}/{total} CIFs "
        f"(seed={args.seed}, dmin={args.dmin}, max_peaks={args.max_peaks})"
    )

    result: dict[str, Any] = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "machine": machine_info(),
        "slice": {
            "cif_root": str(cif_root),
            "archive_total": total,
            "limit": len(tasks),
            "seed": args.seed,
            "dmin": args.dmin,
            "max_peaks": args.max_peaks,
        },
        "notes": [
            "old_two_step = parse_meta(path) + compute_pattern(path) "
            "(the shipped builder hot path)",
            "compute_phase = build_pipeline.compute_phase (single gemmi.cif.read)",
            "sequential timings are paired per file (old then new) and keep the "
            "per-file minimum over repeats, so machine load cancels out",
            "parallel timings include ProcessPoolExecutor start-up; "
            "chunksize is the pool.map chunk size",
        ],
    }

    # ---- sequential comparison (the headline number) ----------------------
    sequential: dict[str, Any] | None = None
    if args.skip_sequential:
        print("\n[1/2] sequential comparison skipped (--skip-sequential)")
    else:
        print("\n[1/2] sequential: old two-step vs compute_phase (paired per file)")
        sequential = bench_sequential(
            tasks, repeats=args.repeats, check_parity=not args.no_parity
        )
        result["sequential"] = sequential
        tw = sequential["time_weighted"]
        print(
            f"  old two-step : {tw['old_two_step']['sum_seconds']:.2f} s  "
            f"({tw['old_two_step']['ms_per_file']:.2f} ms/CIF, "
            f"{tw['old_two_step']['files_per_s']:.2f} CIF/s)"
        )
        print(
            f"  compute_phase: {tw['compute_phase']['sum_seconds']:.2f} s  "
            f"({tw['compute_phase']['ms_per_file']:.2f} ms/CIF, "
            f"{tw['compute_phase']['files_per_s']:.2f} CIF/s)"
        )
        print(
            f"  throughput   : {tw['speedup_x']:.4f}x ({tw['speedup_pct']:+.2f} %)"
            f"  -> {'PASS >=5%' if tw['meets_5pct'] else 'BELOW 5% TARGET'}"
        )
        parse = sequential["parse_cost"]
        print(
            f"  parse only   : old {parse['old_ms_per_file']:.3f} ms/CIF "
            f"({parse['old_share_of_total_pct']:.3f} % of total) -> "
            f"new {parse['new_ms_per_file']:.3f} ms/CIF "
            f"({parse['new_share_of_total_pct']:.3f} %), {parse['speedup_x']:.2f}x faster"
        )
        fast = sequential["fast_phases"]
        slow = sequential["slow_phases"]
        print(
            f"  fast <{SLOW_PHASE_MS:.0f}ms: {fast['files']}/{sequential['files']} files "
            f"({fast['share_of_files_pct']:.1f} %), speedup {fast['speedup_pct']:+.2f} %"
        )
        print(
            f"  slow >={SLOW_PHASE_MS:.0f}ms: {slow['files']}/{sequential['files']} files "
            f"({slow['share_of_files_pct']:.1f} %) but "
            f"{slow['share_of_old_time_pct']:.1f} % of total time, "
            f"speedup {slow['speedup_pct']:+.2f} %"
        )
        parity = sequential["parity"]
        if parity["checked"]:
            print(
                f"  parity       : {parity['checked']} phase comparisons, "
                f"{parity['mismatches']} mismatches"
            )
            if parity["first_mismatch"]:
                print(f"  first mismatch: {parity['first_mismatch']}")

    # ---- parallel grid -----------------------------------------------------
    print("\n[2/2] parallel throughput (old vs compute_phase)")
    parallel: list[dict[str, Any]] = []
    baseline_seconds: dict[tuple[int, int], float] = {}

    old_worker_counts = [args.workers, *(w for w in args.extra_workers if w != args.workers)]
    old_chunk_sizes = list(args.old_chunk_sizes) or [args.chunksize]
    for old_chunksize in old_chunk_sizes:
        for workers in old_worker_counts:
            row = bench_parallel(
                _worker_old, tasks, workers=workers, chunksize=old_chunksize
            )
            row["impl"] = "old_two_step"
            parallel.append(row)
            baseline_seconds[(workers, old_chunksize)] = row["seconds"]
            print(
                f"  old   workers={workers:>2} chunksize={old_chunksize:>2}: "
                f"{row['seconds']:8.2f} s  ({row['files_per_s']:.2f} CIF/s)"
            )

    new_chunk_sizes = [
        args.chunksize,
        *[chunk for chunk in args.chunk_sizes if chunk != args.chunksize],
    ]
    for workers in [args.workers, *args.extra_workers]:
        for chunksize in new_chunk_sizes:
            if workers == args.workers and chunksize == args.chunksize:
                continue
            row = bench_parallel(_worker_new, tasks, workers=workers, chunksize=chunksize)
            row["impl"] = "compute_phase"
            baseline = baseline_seconds.get((workers, args.chunksize))
            row["speedup_vs_old_same_workers_pct"] = (
                None if baseline is None else (baseline / row["seconds"] - 1.0) * 100.0
            )
            parallel.append(row)
            print(
                f"  new   workers={workers:>2} chunksize={chunksize:>2}: "
                f"{row['seconds']:8.2f} s  ({row['files_per_s']:.2f} CIF/s)"
            )

    primary_old = next(
        row
        for row in parallel
        if row["impl"] == "old_two_step" and row["workers"] == args.workers
    )
    primary_new = bench_parallel(
        _worker_new, tasks, workers=args.workers, chunksize=args.chunksize
    )
    primary_new["impl"] = "compute_phase"
    primary_speedup = (primary_old["seconds"] / primary_new["seconds"] - 1.0) * 100.0
    primary_new["speedup_vs_old_same_workers_pct"] = primary_speedup
    parallel.append(primary_new)
    print(
        f"  primary workers={args.workers} chunksize={args.chunksize}: "
        f"old {primary_old['seconds']:.2f} s -> new {primary_new['seconds']:.2f} s "
        f"({primary_speedup:+.2f} %)"
    )
    result["parallel"] = {
        "primary": {
            "workers": args.workers,
            "chunksize": args.chunksize,
            "old_two_step": primary_old,
            "compute_phase": primary_new,
            "speedup_pct": primary_speedup,
        },
        "rows": parallel,
    }

    # ---- persist + verdict -------------------------------------------------
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\nJSON written to {out}")

    if sequential is not None:
        if sequential["parity"]["mismatches"]:
            print(
                f"FAIL: {sequential['parity']['mismatches']} parity mismatches "
                "between old and new paths"
            )
            return 1
        if not sequential["time_weighted"]["meets_5pct"]:
            print(
                "NOTE: measured whole-slice throughput speedup is below the 5 % "
                "target; the single-parse code is kept because it is semantically "
                "cleaner (and measurably faster per phase on typical files)."
            )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse the CLI and run the benchmark."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--cif-root",
        default=str(DEFAULT_INORG_COD_ROOT),
        help="CIF archive root (default: the real inorganic COD subset)",
    )
    ap.add_argument("--limit", type=int, default=800, help="CIFs sampled (500-1000 recommended)")
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED, help="deterministic sample seed")
    ap.add_argument("--dmin", type=float, default=DEFAULT_DMIN, help="peak engine dmin (A)")
    ap.add_argument("--max-peaks", type=int, default=DEFAULT_MAX_PEAKS, help="peak engine max_peaks")
    ap.add_argument("--repeats", type=int, default=2, help="paired sequential repetitions")
    ap.add_argument("--workers", type=int, default=4, help="primary worker count (task target 4)")
    ap.add_argument(
        "--extra-workers",
        type=int,
        nargs="*",
        default=[8, 14],
        help="additional worker counts for the parallel grid",
    )
    ap.add_argument(
        "--chunksize",
        type=int,
        default=16,
        help="primary ProcessPool chunksize (builder default 8; 8/16/32 all measured)",
    )
    ap.add_argument(
        "--chunk-sizes",
        type=int,
        nargs="*",
        default=[8, 16, 32],
        help="additional chunksizes for the parallel grid",
    )
    ap.add_argument(
        "--old-chunk-sizes",
        type=int,
        nargs="*",
        default=[],
        help="chunksizes the old two-step path is measured at "
        "(default: only --chunksize); pass 8 16 32 for a same-config comparison",
    )
    ap.add_argument(
        "--skip-sequential",
        action="store_true",
        help="skip the sequential/parity pass and benchmark only the parallel grid",
    )
    ap.add_argument(
        "--out", default=str(DEFAULT_OUT), help="JSON output path (under tmp/ by default)"
    )
    ap.add_argument(
        "--no-parity", action="store_true", help="skip the bit-exact parity fingerprint check"
    )
    args = ap.parse_args(argv)
    return run(args)


if __name__ == "__main__":
    raise SystemExit(main())
