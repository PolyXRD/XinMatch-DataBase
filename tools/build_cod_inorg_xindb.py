"""Build a v3 **inorganic** ``.xinDB`` library from COD, PCOD or TCOD (M2 / S2.4).

Thin wrapper around :mod:`tools.build_cod_xindb` that pins the CIF source to one
of the three *inorganic subset* archives and stamps the matching library kind.

Why the inorganic subset is the interactive library
---------------------------------------------------
The full COD archive holds thousands of near-duplicate cards for the same
phase, and those duplicates spread the FoM vote so thin that the top-1 hit rate
collapses as the library grows.  The subset keeps the ranking honest while the
full archive stays available as the exhaustive reference.

The three sources
-----------------
================  =======================================  =========  =====================
``--source``      archive                                  CIF count  entry-id shift
================  =======================================  =========  =====================
``COD``           ``inorg_cif/COD``  (experimental)          79 198   ``0``
``PCOD``          ``inorg_cif/PCOD`` (predicted, GRINSP)  1 061 962  ``+20 000 000``
``TCOD``          ``inorg_cif/TCOD`` (theoretical)            46 943  ``+100 000 000``
================  =======================================  =========  =====================

The shift is not cosmetic.  PCOD ids occupy the same 7-digit numeric space as
COD ids, so PCOD ``1000000`` and COD ``1000000`` would both become
``96-100-0000`` — one library would silently shadow the other in any merged
view.  Shifting keeps every entry id unique while preserving the ``96-XXX-XXXX``
shape the GUI, the baselines and the FoM tables already key on.  TCOD ids are 8
digits to begin with, which is why the shift is applied before formatting and
the formatter accepts wider ids.

Usage
-----
::

    # verification slice (always do this before a full build)
    python -m tools.build_cod_inorg_xindb --source PCOD --limit 500 --out /tmp/pcod_slice.xinDB

    # the three interactive libraries
    python -m tools.build_cod_inorg_xindb --source COD  --data-date 20261020 --workers 14 -v
    python -m tools.build_cod_inorg_xindb --source PCOD --data-date 20261020 --workers 14 -v
    python -m tools.build_cod_inorg_xindb --source TCOD --data-date 20261020 --workers 14 -v

Output defaults to ``$XINMATCH_XINDB_DIR/<kind>-<data-date>.xinDB`` when that
variable is set, otherwise ``<repo>/data/...``.
"""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from xinmatch.core.db.peak_compute import DEFAULT_DMIN, DEFAULT_MAX_PEAKS  # noqa: E402

from tools.build_cod_xindb import build, default_out_dir  # noqa: E402
from xinmatch.core.db.cif_scan import (  # noqa: E402
    DEFAULT_INORG_COD_ROOT,
    DEFAULT_INORG_PCOD_ROOT,
    DEFAULT_INORG_TCOD_ROOT,
)

logger = logging.getLogger("tools.build_cod_inorg_xindb")

#: ``source -> (inorganic archive root, library kind, entry-id shift)``.
#: The shift follows the 3COD convention so entry ids stay compatible with the
#: historical baselines: COD keeps its own ids, PCOD +20 000 000, TCOD
#: +100 000 000.
SOURCES: dict[str, tuple[Path, str, int]] = {
    "COD": (DEFAULT_INORG_COD_ROOT, "COD-Inorg", 0),
    "PCOD": (DEFAULT_INORG_PCOD_ROOT, "PCOD-Inorg", 20_000_000),
    "TCOD": (DEFAULT_INORG_TCOD_ROOT, "TCOD-Inorg", 100_000_000),
}

#: Default source, kept as ``COD`` so existing invocations behave unchanged.
DEFAULT_SOURCE = "COD"


def main(argv: list[str] | None = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--source", choices=sorted(SOURCES), default=DEFAULT_SOURCE,
        help=f"which archive to build from (default {DEFAULT_SOURCE})",
    )
    parser.add_argument(
        "--out", type=Path, default=None,
        help="output .xinDB path (may live in any folder); defaults to "
             "$XINMATCH_XINDB_DIR/<kind>-<data-date>.xinDB when that variable "
             "is set, otherwise <repo>/data/<kind>-<data-date>.xinDB",
    )
    parser.add_argument(
        "--data-date", default=None,
        help="source snapshot date (YYYYMMDD) recorded in meta.data_date and used "
             "in the default output name; defaults to today",
    )
    parser.add_argument(
        "--cif-root", type=Path, default=None,
        help="CIF archive root; defaults to the inorganic subset of --source",
    )
    parser.add_argument("--dmin", type=float, default=DEFAULT_DMIN)
    parser.add_argument("--max-peaks", type=int, default=DEFAULT_MAX_PEAKS)
    parser.add_argument("--start", type=int, default=0,
                        help="lowest archive id to process (pre-offset)")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--tier", type=int, default=1)
    parser.add_argument("--log-every", type=int, default=500)
    parser.add_argument("--workers", type=int, default=1,
                        help="peak-computation processes (1 = in-process)")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args(argv)

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    default_root, library_kind, id_offset = SOURCES[args.source]
    cif_root = args.cif_root or default_root
    if not cif_root.is_dir():
        print(f"[ERR] CIF root not found: {cif_root}")
        return 2

    data_date = (args.data_date or time.strftime("%Y%m%d")).strip()
    out = args.out or (default_out_dir() / f"{library_kind}-{data_date}.xinDB")

    print(f"kind:     {library_kind}")
    print(f"source:   {args.source}")
    print(f"data:     {data_date}")
    print(f"cif root: {cif_root}")
    print(f"id offset:{id_offset:+,}")
    print(f"out:      {out}")
    print()

    stats = build(
        cif_root,
        out,
        library_kind=library_kind,
        build_tool="tools.build_cod_inorg_xindb",
        data_date=data_date,
        dmin=args.dmin,
        max_peaks=args.max_peaks,
        start=args.start,
        limit=args.limit,
        resume=args.resume,
        tier=args.tier,
        log_every=args.log_every,
        workers=args.workers,
        id_offset=id_offset,
    )

    size_mb = out.stat().st_size / 1e6 if out.exists() else 0.0
    print("\nbuild summary:")
    for key, value in stats.as_dict().items():
        print(f"  {key}: {value}")
    print(f"  size_mb: {size_mb:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
