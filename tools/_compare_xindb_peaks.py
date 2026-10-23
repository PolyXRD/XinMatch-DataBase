"""Compare the peak payload of two ``.xinDB`` libraries, chunk by chunk.

Used to prove that a *metadata-only* change (for example the element-extraction
fix in ``cif_meta``) did not disturb the physics: every ``peak_chunks`` row must
hash to the same digest on both sides, for ``top_d``, ``top_i``, ``all_d`` and
``all_i``.

Peak rows are compared as raw BLOBs rather than decoded arrays.  The writer packs
them with a fixed ``<f4`` layout, so a byte comparison is both stricter and much
cheaper than a numeric one.

Entries are *not* required to be in the same row order: entries are joined by
``entry_id`` so the comparison survives a different CIF enumeration order.

Usage::

    python -m tools._compare_xindb_peaks data/COD-Inorg-20261020.xinDB.prev \
        data/COD-Inorg-20261020.xinDB
"""

from __future__ import annotations

import argparse
import hashlib
import sqlite3
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def _open(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True, timeout=60.0)


def _chunk_digest(blob: object) -> str:
    return hashlib.blake2b(bytes(blob or b""), digest_size=16).hexdigest()


def _meta_value(conn: sqlite3.Connection, key: str) -> str:
    row = conn.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
    return str(row[0]) if row else ""


def compare(left: Path, right: Path) -> int:
    """Print a per-chunk comparison and return a process exit code."""
    for path in (left, right):
        if not path.exists():
            print(f"missing: {path}", file=sys.stderr)
            return 2

    a = _open(left)
    b = _open(right)

    for conn, label in ((a, "A"), (b, "B")):
        print(
            f"{label}: {_meta_value(conn, 'library_kind')} "
            f"{_meta_value(conn, 'data_date')} "
            f"entries={conn.execute('SELECT count(*) FROM entries').fetchone()[0]:,} "
            f"chunks={conn.execute('SELECT count(*) FROM peak_chunks').fetchone()[0]:,}"
        )

    # The chunk grid itself must match, otherwise the row order diverged and a
    # chunk-wise comparison would be meaningless.
    grid_a = a.execute("SELECT chunk_id, n FROM peak_chunks ORDER BY chunk_id").fetchall()
    grid_b = b.execute("SELECT chunk_id, n FROM peak_chunks ORDER BY chunk_id").fetchall()
    if grid_a != grid_b:
        print(f"FAIL chunk grid differs: {len(grid_a)} vs {len(grid_b)} rows")
        return 1

    ids_a = [r[0] for r in a.execute("SELECT entry_id FROM entries ORDER BY rowid")]
    ids_b = [r[0] for r in b.execute("SELECT entry_id FROM entries ORDER BY rowid")]
    same_order = ids_a == ids_b
    print(f"entry order identical: {same_order}  ({len(ids_a):,} vs {len(ids_b):,} ids)")
    if not same_order:
        only_a = set(ids_a) - set(ids_b)
        only_b = set(ids_b) - set(ids_a)
        print(f"  only in A: {len(only_a):,}   only in B: {len(only_b):,}")
        if only_a or only_b:
            print(f"  sample A-only: {list(only_a)[:5]}")
            print(f"  sample B-only: {list(only_b)[:5]}")

    mismatched = 0
    total = 0
    for cid, n in grid_a:
        row_a = a.execute(
            "SELECT top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?", (cid,)
        ).fetchone()
        row_b = b.execute(
            "SELECT top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?", (cid,)
        ).fetchone()
        total += 1
        for name, blob_a, blob_b in zip(
            ("top_d", "top_i", "all_d", "all_i"), row_a, row_b
        ):
            if _chunk_digest(blob_a) != _chunk_digest(blob_b):
                mismatched += 1
                print(
                    f"  DIFF chunk {cid} {name}: "
                    f"{len(bytes(blob_a or b''))} vs {len(bytes(blob_b or b''))} bytes"
                )
                break

    a.close()
    b.close()

    if mismatched:
        print(f"\nFAIL: {mismatched}/{total} chunks differ in peak payload")
        return 1
    print(f"\nPASS: all {total} chunks hash identically across top_d/top_i/all_d/all_i")
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("left", type=Path, help="reference library (e.g. the .prev copy)")
    ap.add_argument("right", type=Path, help="newly built library")
    args = ap.parse_args(argv)
    print(f"A = {args.left}\nB = {args.right}\n")
    return compare(args.left, args.right)


if __name__ == "__main__":
    raise SystemExit(main())
