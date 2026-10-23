"""``.xinDB`` v3 schema — columnar reference-database layout (M2 / S2.1).

Why a new format at all
-----------------------
The v2 reference databases store every reflection as **one row** in a
``peaks`` table plus three B-tree indexes over it.  Measured on the shipped
``cod_inorganics.db`` (745 MB, 71,206 entries, 6.21 M peaks):

===========================  ==========  =====
object                       size        share
===========================  ==========  =====
``peaks`` row table          283.8 MB    38 %
``idx_peaks_entry_intensity`` 176.7 MB   24 %
``idx_peaks_entry_id``       125.9 MB    17 %
``idx_peaks_d``              119.9 MB    16 %
everything else              ~39 MB      5 %
===========================  ==========  =====

So 95 % of the file is the peak table and its indexes — the *same* data
stored three times in three different shapes.  The v3 layout stores it
**once**, in column-major binary blocks, and keeps only the two indexes that
actually pay for themselves.

Layout
------
::

    entries      one row per phase (scalars + CSR offsets into peak_chunks)
    peak_chunks  fixed-size blocks of columnar peak data
    meta         self-describing key/value table

``peak_chunks`` holds two representations of the same reflections:

``top_d`` / ``top_i``
    A dense ``(n, TOP_K)`` float32 matrix of each phase's *strongest*
    reflections.  This is the **prescreen** matrix: the search loads it once
    and intersects measurement windows against it with plain NumPy, so
    candidate generation never touches the full peak data and never needs a
    runtime index build.

``all_d`` / ``all_i``
    The complete reflection list, stored as a **CSR-style** float32 array:
    phase *i* owns the slice ``all_d[peak_off[i] : peak_off[i] + peak_n[i]]``.
    Only the (few hundred) surviving candidates are expanded from it.

Both use little-endian float32 (``<f4``).  The d-spacing precision this
implies is ~1e-7 Å at the relevant magnitudes — three orders of magnitude
below the 0.01–0.02 Å matching windows — so the float32 rounding cannot
change a match decision.  This is verified by the parity test.

Design rules that protect accuracy
----------------------------------
The database layout only decides *where the numbers live*; it never decides
*how a candidate is scored*.  Ranking stays in
:mod:`xinmatch.core.search` untouched.  Two rules follow from that:

* The prescreen matrix is a **pruning** aid, not a ranking input — the final
  ordering always uses the full reflections from ``all_d`` / ``all_i``.
* The prescreen candidate cap is generous by design; ``prescreen()`` reports
  when it truncates so the caller can widen it rather than silently lose a
  true phase.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path

from xinmatch.core.db.peak_compute import PEAK_CONVENTION

#: Schema version written to ``meta.schema_version``.
SCHEMA_VERSION = 3

#: File extension of a v3 database.  Deliberately distinct from the v2
#: ``.db`` / ``.sqlite`` spellings so the two can never be confused.
XINDB_SUFFIX = ".xinDB"

#: ``XINDB_SUFFIX`` in lower case, for comparing against ``Path.suffix.lower()``.
#: Comparing the mixed-case constant directly would be wrong: the left side is
#: always lower-cased, so ``suffix.lower() == XINDB_SUFFIX`` is false for every
#: file and no ``.xinDB`` is ever recognised.  Keep the display spelling in
#: ``XINDB_SUFFIX`` (it is used to *build* file names) and compare with this.
XINDB_SUFFIX_LC = XINDB_SUFFIX.lower()

#: Number of strongest reflections kept in the prescreen matrix per phase.
#: 32 is comfortably above the 12–40 strongest reflections the scoring
#: actually consumes, so the prescreen never truncates a phase's
#: discriminating lines.
TOP_K = 32

#: Number of phases per ``peak_chunks`` row.  4096 keeps each block in the
#: low-megabyte range: large enough that the per-block overhead disappears,
#: small enough that a reader only has to mmap a few blocks for the
#: candidates it actually needs.
CHUNK_SIZE = 4096

#: Peak-intensity normalisation used at build time (the project convention
#: inherited from the reference software: I/I0 with I0 = 1000).
INTENSITY_SCALE = 1000.0

#: The peak-calculation convention recorded in ``meta.peak_convention``.
#: A reader that finds an unknown convention must refuse the file rather
#: than silently mix incompatible data.
#:
#: Re-exported from :mod:`peak_compute` on purpose.  The v3.0.0 rebuild
#: originally carried a second, stale literal here
#: (``CuKa1-1.5406-dmin-1.0-norm-1000``); because the writer stamps the value
#: from this module, every built library claimed a ``dmin = 1.0`` cutoff while
#: the engine actually computed down to 0.7703 A.  One definition, one truth.

#: The canonical DDL.  Kept as a module-level string (rather than a
#: ``.sql`` resource) so the schema travels with the code that writes it —
#: the v3 format is new, and a resource file would add a packaging
#: failure mode for no benefit.
SCHEMA_V3_SQL = """
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);

CREATE TABLE IF NOT EXISTS entries (
    entry_id       TEXT PRIMARY KEY,
    formula        TEXT,
    space_group    TEXT,
    crystal_system TEXT,
    mineral_name   TEXT,
    elements       TEXT,
    elem_mask      INTEGER,
    a REAL, b REAL, c REAL,
    alpha REAL, beta REAL, gamma REAL,
    volume  REAL,
    density REAL,
    year    INTEGER,
    author  TEXT,
    journal TEXT,
    quality TEXT,
    n_peaks   INTEGER,
    peak_off  INTEGER,
    peak_n    INTEGER,
    tier      INTEGER DEFAULT 1
);

CREATE TABLE IF NOT EXISTS peak_chunks (
    chunk_id INTEGER PRIMARY KEY,
    n        INTEGER,
    top_d    BLOB,
    top_i    BLOB,
    all_d    BLOB,
    all_i    BLOB
);

CREATE TABLE IF NOT EXISTS element_inv (
    element TEXT PRIMARY KEY,
    eids    BLOB
);

CREATE INDEX IF NOT EXISTS idx_entries_space_group
    ON entries(space_group);
CREATE INDEX IF NOT EXISTS idx_entries_crystal_system
    ON entries(crystal_system);
CREATE INDEX IF NOT EXISTS idx_entries_formula
    ON entries(formula);
CREATE INDEX IF NOT EXISTS idx_entries_tier
    ON entries(tier);
"""

#: Tables a v3 database must have for a reader to accept it.
REQUIRED_TABLES: tuple[str, ...] = ("meta", "entries", "peak_chunks")

#: ``meta`` keys written by the builder.
META_KEYS: tuple[str, ...] = (
    "schema_version",
    "peak_convention",
    "top_k",
    "chunk_size",
    "dtype",
    "entry_count",
    "build_time",
    "cif_source",
    "build_tool",
    "library_kind",
)


def create_tables(conn: sqlite3.Connection) -> None:
    """Create the v3 tables on *conn* (idempotent).

    Every statement uses ``IF NOT EXISTS``, so calling this on an existing
    database is harmless and can be used to upgrade a partially built file.
    """
    conn.executescript(SCHEMA_V3_SQL)


def set_meta(conn: sqlite3.Connection, key: str, value: object) -> None:
    """Write one ``meta`` key/value pair."""
    conn.execute(
        "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
        (str(key), str(value)),
    )


def set_meta_many(conn: sqlite3.Connection, values: dict[str, object]) -> None:
    """Write several ``meta`` key/value pairs in one transaction."""
    with conn:
        conn.executemany(
            "INSERT OR REPLACE INTO meta (key, value) VALUES (?, ?)",
            [(str(k), str(v)) for k, v in sorted(values.items())],
        )


def read_meta(conn: sqlite3.Connection) -> dict[str, str]:
    """Return the whole ``meta`` table as a dict."""
    return {
        str(k): str(v)
        for k, v in conn.execute("SELECT key, value FROM meta")
    }


def write_standard_meta(
    conn: sqlite3.Connection,
    *,
    entry_count: int,
    cif_source: str,
    build_time: str,
    library_kind: str,
    build_tool: str = "tools.build_cod_xindb",
    data_date: str = "",
    top_k: int = TOP_K,
    chunk_size: int = CHUNK_SIZE,
) -> None:
    """Write the standard ``meta`` set describing a freshly built library.

    ``library_kind`` distinguishes the two shipped products (``"COD-Full"`` /
    ``"COD-Inorg"``) so a reader can tell an inorganic subset from the full
    archive without inspecting the contents.

    ``data_date`` is the **COD snapshot date** the CIFs came from (``YYYYMMDD``),
    which is *not* the same thing as ``build_time`` — the same snapshot may be
    rebuilt with a different peak engine, and a later snapshot may be built the
    same day.  Keeping them separate is what makes "which data is this?" and
    "when was it computed?" independently answerable.
    """
    values = {
        "schema_version": SCHEMA_VERSION,
        "peak_convention": PEAK_CONVENTION,
        "top_k": top_k,
        "chunk_size": chunk_size,
        "dtype": "f32",
        "entry_count": entry_count,
        "build_time": build_time,
        "cif_source": cif_source,
        "build_tool": build_tool,
        "library_kind": library_kind,
    }
    if data_date:
        values["data_date"] = str(data_date)
    set_meta_many(conn, values)


def is_xindb_path(path: str | Path) -> bool:
    """True when *path* carries the v3 extension (cheap pre-check)."""
    return Path(path).suffix.lower() == XINDB_SUFFIX_LC


__all__ = [
    "CHUNK_SIZE",
    "INTENSITY_SCALE",
    "META_KEYS",
    "PEAK_CONVENTION",
    "REQUIRED_TABLES",
    "SCHEMA_VERSION",
    "SCHEMA_V3_SQL",
    "TOP_K",
    "XINDB_SUFFIX",
    "XINDB_SUFFIX_LC",
    "create_tables",
    "is_xindb_path",
    "read_meta",
    "set_meta",
    "set_meta_many",
    "write_standard_meta",
]
