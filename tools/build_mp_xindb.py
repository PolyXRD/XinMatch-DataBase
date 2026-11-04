"""Build ``MP-<date>.xinDB`` — the Materials Project external reference library.

Unlike the COD/PCOD/TCOD archives, the Materials Project summary is a flat
directory of formula-named CIFs whose stable identifier lives in a leading
comment (``# mp-id: mp-862690``).  This tool reuses the six-library build
pipeline (:func:`tools.build_cod_xindb.build`) with ``id_mode="mp"``:

* ``entry_id`` = ``mp-<number>`` (the MP native id);
* real space group / crystal system taken from the CIF comment header
  (pymatgen writes the structure as P1 and keeps the original group only
  there), with the crystal system derived from the H-M symbol;
* ``tier = 2`` (computed / theoretical structure), journal/author/quality
  recorded as ``Materials Project`` / ``Materials Project`` / ``DFT``;
* ``meta.library_kind = "MP"``, ``meta.data_date = <CIF download date>``;
* skip the 3 archive files that carry no ``mp-id`` comment.

Usage::

    python -m tools.build_mp_xindb --limit 1000 --out tmp/accept/MP-slice.xinDB
    python -m tools.build_mp_xindb --data-date 20261031 --workers 14 -v
"""

from __future__ import annotations

import argparse
import logging
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.build_cod_xindb import build, default_out_dir  # noqa: E402
from xinmatch.core.db.mp_scan import (  # noqa: E402
    DEFAULT_MP_ROOT,
    scan_mp_archive,
)
from xinmatch.core.db.schema_v3 import XINDB_SUFFIX  # noqa: E402

logger = logging.getLogger("tools.build_mp_xindb")

#: CIF download date recorded in ``meta.data_date`` and the file name.
DEFAULT_DATA_DATE = "20261031"


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--cif-root",
        type=Path,
        default=DEFAULT_MP_ROOT,
        help=f"Materials Project CIF directory (default {DEFAULT_MP_ROOT})",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="output .xinDB path; defaults to <repo>/data/MP-<data-date>.xinDB",
    )
    parser.add_argument(
        "--data-date",
        default=DEFAULT_DATA_DATE,
        help=f"CIF download date (default {DEFAULT_DATA_DATE})",
    )
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--workers",
        type=int,
        default=14,
        help="peak-computation processes (default 14)",
    )
    parser.add_argument("--log-every", type=int, default=5000)
    parser.add_argument("-v", "--verbose", action="store_true")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if not args.cif_root.is_dir():
        print(f"[ERR] MP CIF root not found: {args.cif_root}")
        return 2

    data_date = str(args.data_date).strip()
    out = args.out or (default_out_dir() / f"MP-{data_date}{XINDB_SUFFIX}")

    print("kind:     MP")
    print(f"data:     {data_date}")
    print(f"cif root: {args.cif_root}")
    print(f"out:      {out}")

    t0 = time.perf_counter()
    archive = scan_mp_archive(args.cif_root)
    print(
        f"archive:  {archive.total:,} CIF, {len(archive.entries):,} with mp-id, "
        f"{len(archive.no_id):,} without (skipped)"
    )
    if archive.no_id:
        for path in archive.no_id:
            print(f"  [skip] {path.name}  (no '# mp-id:' header)")
    print()

    stats = build(
        args.cif_root,
        Path(out),
        library_kind="MP",
        build_tool="tools.build_mp_xindb",
        data_date=data_date,
        start=args.start,
        limit=args.limit,
        resume=args.resume,
        tier=2,
        log_every=args.log_every,
        workers=args.workers,
        id_mode="mp",
        extra_meta={
            "dataset": "MP",
            "mp_archive_cifs": str(archive.total),
            "mp_with_id": str(len(archive.entries)),
            "mp_without_id": str(len(archive.no_id)),
            "mp_archive_root": str(args.cif_root),
        },
    )

    size_mb = Path(out).stat().st_size / 1e6 if Path(out).exists() else 0.0
    print("\nbuild summary:")
    for key, value in stats.as_dict().items():
        print(f"  {key}: {value}")
    print(f"  size_mb: {size_mb:.1f}")
    print(f"  wall_sec: {time.perf_counter() - t0:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
