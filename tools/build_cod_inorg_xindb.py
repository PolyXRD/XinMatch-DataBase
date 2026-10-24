"""Build a v3 **inorganic** ``.xinDB`` library from COD, PCOD or TCOD (M2 / S2.4).

Thin compatibility wrapper around :mod:`tools.build_cod_xindb`.  The unified
builder already understands ``--source {COD,PCOD,TCOD}`` + ``--set inorg``; this
script only pins the default subset to ``inorg`` and records
``tools.build_cod_inorg_xindb`` in ``meta.build_tool``, so existing invocations
(and the published libraries' provenance) stay valid.  There is exactly one
build implementation and one source registry -- no second code path.

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
``96-100-0000`` -- one library would silently shadow the other in any merged
view.  Shifting keeps every entry id unique while preserving the ``96-XXX-XXXX``
shape the GUI, the baselines and the FoM tables already key on.  TCOD ids are 8
digits to begin with, which is why the shift is applied before formatting and
the formatter accepts wider ids.

Usage
-----
::

    # verification slice (always do this before a full build)
    python -m tools.build_cod_inorg_xindb --source PCOD --limit 500 --out tmp/accept/pcod_slice.xinDB

    # the three interactive libraries
    python -m tools.build_cod_inorg_xindb --source COD  --data-date 20261020 --workers 14 -v
    python -m tools.build_cod_inorg_xindb --source PCOD --data-date 20261020 --workers 14 -v
    python -m tools.build_cod_inorg_xindb --source TCOD --data-date 20261020 --workers 14 -v

The same builder also produces the three *full* libraries through
``python -m tools.build_cod_xindb --source {...} --set full``.

Output defaults to ``<repo>/data/<kind>-<data-date>.xinDB``; ``XINMATCH_XINDB_DIR``
is an explicit override of that directory.  Every other flag (``--out``,
``--limit``, ``--start``, ``--resume``, ``--workers``, ``--tier``,
``--log-every``, ``--data-date``, ``-v``, ``--cif-root``, ``--id-offset``) is
the unified builder's.
"""

from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.build_cod_xindb import main as _build_main  # noqa: E402

#: ``meta.build_tool`` recorded for libraries built through this wrapper.  Kept
#: stable so the published inorganic libraries and their rebuilds agree.
BUILD_TOOL = "tools.build_cod_inorg_xindb"


def main(argv: list[str] | None = None) -> int:
    """Delegate to the unified builder, defaulting ``--set`` to ``inorg``."""
    return _build_main(argv, default_subset="inorg", build_tool=BUILD_TOOL)


if __name__ == "__main__":
    sys.exit(main())
