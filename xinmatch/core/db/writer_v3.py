"""Columnar writer for ``.xinDB`` v3 reference databases (M2 / S2.2).

What this module does
---------------------
The v3 rebuild stores every reflection **once**, in column-major binary
blocks, instead of one SQLite row per reflection.  This module is the
producer side of that format: it takes phases one at a time and packs them
into a ``.xinDB`` file that :mod:`xinmatch.core.db.reader_v3` can read back
without ever building a runtime index.

Storage layout produced here
----------------------------
::

    meta(key, value)          self-describing build information
    entries(...)              one row per phase (scalars + CSR offsets)
    peak_chunks(chunk_id, n,
                top_d, top_i, all_d, all_i)
    element_inv(element, eids)

``peak_chunks`` carries two representations of the same reflections:

``top_d`` / ``top_i``
    A dense ``(n, top_k)`` float32 matrix of each phase's *strongest*
    reflections — the **prescreen** matrix.  Phases with fewer than
    ``top_k`` reflections are right-padded with ``0.0`` so the block stays
    rectangular and can be decoded with a single ``reshape``.

``all_d`` / ``all_i``
    The complete reflection list, stored CSR-style: phase *i* owns the
    slice ``all_d[peak_off[i] : peak_off[i] + peak_n[i]]`` of the
    *conceptual* global array.  The global array is not stored as one
    blob; it is cut into per-chunk pieces, and a chunk's ``all_d`` holds
    the reflections of the phases belonging to that chunk, in phase order.

Why ``peak_off`` is a *peak index* and not a byte offset
---------------------------------------------------------
``peak_off`` counts **reflections** in the conceptual global array, not
bytes.  A byte offset would be ambiguous: the reader would have to know
the chunk base to turn it into a slice, which is exactly the information
the offset was supposed to provide.  With a reflection index the reader
walks the chunks in ``chunk_id`` order, accumulates the reflection count
of each chunk, and the first chunk whose range covers ``peak_off`` holds
the phase; the phase's slice then starts at ``peak_off - chunk_base``.
That derivation needs no extra bookkeeping on the write side, and the
reflection index is also the natural unit for "how many peaks are stored".

Consequences the writer guarantees
-----------------------------------
* Reflections are appended to the current chunk in phase order, and a
  chunk is written to SQLite as soon as it is full — the peak data of a
  whole 535 k-phase build never sits in memory at once.
* ``entries`` rows are inserted in the same order the phases were added,
  and the integer "entry index" used by ``element_inv`` is that same
  position.  The index is written as the SQLite ``rowid`` explicitly, so
  ``SELECT rowid FROM entries WHERE entry_id = ?`` returns the stored
  number unchanged and the two representations can never drift apart.
* ``elem_mask`` and ``element_inv`` are produced by the *same* element
  mask function, so the fast bitmask filter and the explicit inverted
  index can never disagree about which phases contain an element.  The
  bit table itself comes from :mod:`xinmatch.core.db.cif_meta` alone —
  the bit layout is baked into the stored data, so the writer and the
  reader must share one table.  The table that was used is recorded in
  ``meta.element_bits``, which makes every file self-describing.

How wide the element mask can be
--------------------------------
``entries.elem_mask`` is a single SQLite ``INTEGER`` column, and SQLite
integers are *signed* 64-bit values (the Python driver refuses anything
above ``2**63 - 1`` outright).  A bitmask stored there can therefore
carry bits 0..62 and nothing more — 63 elements.

The canonical element table in :mod:`cif_meta` is wider than that (it
covers all 118 IUPAC symbols in 128 bits), so a phase may contain
elements whose bit does not fit into the column.  Those elements are
**not** silently mangled: the bitmask is truncated, while

* ``entries.elements`` keeps the exact element list, and
* ``element_inv`` keeps the exact inverted index over *all* elements,

so the file never loses the information — only the fast bitmask path
cannot express it.  A file states what it did in ``meta``:

``meta.element_mask_scheme``
    Version tag of the bit table (from :mod:`cif_meta`).
``meta.element_mask_bits``
    Number of bits the stored mask actually uses (63).
``meta.element_mask_full_bits``
    Width of the table the bit numbers come from (128).
``meta.element_mask_unmapped``
    The elements whose bit does not fit into the column, as
    ``"Symbol:bit;..."``.

The truncation is *sound* for the "required elements" semantics the
reader implements (a phase keeps a required element only if its bit is
present), but it is **not** sound for "subset of an allowed pool"
semantics, where the exact ``element_inv`` index has to be used.

Precision and scale
-------------------
Both peak arrays are stored as little-endian ``float32``.  The rounding
error of a float32 grows with the magnitude of the value, so it is
quoted relative to the d-spacing: ``|error| <= d * 6e-8``.  Measured on
the shipped COD inorganic subset, the largest d-spacing is 42.2 Å and
the largest absolute error is 1.7e-6 Å — four orders of magnitude below
the 0.01–0.02 Å matching windows, so the rounding cannot change a match
decision.

Intensities are written on the ``max = 1000`` scale that the file
declares in ``meta.peak_convention``.  The upstream peak calculator
already normalises, which makes the writer's own normalisation a no-op
there, but enforcing it here means the file can never contradict the
convention it advertises.  The d-cutoff (>= 1.0 Å) belongs to the peak
calculator and is not touched here.
"""

from __future__ import annotations

import contextlib
import logging
import sqlite3
from array import array
from collections.abc import Callable, Iterable, Sequence
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, NamedTuple

import numpy as np

from .schema_v3 import (
    CHUNK_SIZE,
    INTENSITY_SCALE,
    TOP_K,
    create_tables,
    set_meta,
    write_standard_meta,
)

if TYPE_CHECKING:  # pragma: no cover - import cycle-free typing aid
    from .cif_meta import PhaseMeta

logger = logging.getLogger(__name__)

__all__ = [
    "XinDBWriter",
    "pack_float32",
    "pack_uint32",
]

#: Number of ``entries`` rows buffered before an ``executemany`` flush.
#: 2000 rows is small enough to keep the Python-side memory flat and big
#: enough that the per-statement SQLite overhead stays negligible.
_ENTRY_BATCH = 2000

#: Largest bitmask a single SQLite ``INTEGER`` column can hold *exactly*.
#: SQLite stores integers as signed 64-bit values and the Python driver
#: refuses to bind anything above ``2**63 - 1`` (it raises
#: ``OverflowError``), so bits 0..62 are usable and bit 63 is not.
_MAX_STORABLE_MASK = (1 << 63) - 1

#: Number of mask bits that fit into ``entries.elem_mask``.
_STORABLE_MASK_BITS = 63


def pack_float32(values: Iterable[float] | np.ndarray) -> bytes:
    """Pack *values* into little-endian float32 bytes.

    The on-disk representation of every float array in a ``.xinDB`` file
    is ``<f4``; producing it through this single helper keeps the
    byte order and the item size identical everywhere.

    Parameters
    ----------
    values:
        Any 1-D sequence of numbers, or a NumPy array.

    Returns
    -------
    bytes
        The packed buffer (``4 * len(values)`` bytes).
    """
    return np.asarray(values, dtype="<f4").tobytes()


def pack_uint32(values: Iterable[int] | np.ndarray) -> bytes:
    """Pack *values* into little-endian uint32 bytes.

    Used for the ``element_inv`` inverted index, where the payload is a
    list of integer entry indices.

    Parameters
    ----------
    values:
        Any 1-D sequence of integers, or a NumPy array.

    Returns
    -------
    bytes
        The packed buffer (``4 * len(values)`` bytes).
    """
    return np.asarray(values, dtype="<u4").tobytes()


class _ElementBits(NamedTuple):
    """The element bit table shared by the writer and the reader.

    Attributes
    ----------
    mask_fn:
        Maps a sequence of element symbols onto a bitmask.  The
        canonical implementation in :mod:`cif_meta` returns a 128-bit
        ``(low, high)`` pair; a 64-bit implementation may return a
        plain integer.  Both shapes are accepted.
    positions:
        Element symbol -> bit number.
    scheme:
        Version tag of the table, recorded in ``meta.element_mask_scheme``.
    full_bits:
        Width of the table the bit numbers come from.
    """

    mask_fn: Callable[[Sequence[str]], int | tuple[int, int]]
    positions: dict[str, int]
    scheme: str
    full_bits: int


def _resolve_element_mask() -> _ElementBits:
    """Load the element mask function and bit table from :mod:`cif_meta`.

    :mod:`xinmatch.core.db.cif_meta` is the **single source of truth**
    for the element bit layout.  The reader imports the very same
    ``element_bit_positions`` table to translate a user's element
    restriction into a bitmask, so the writer and the reader must agree
    bit for bit: the bit layout is *baked into the stored data*
    (``entries.elem_mask`` holds the bits), and a library written with
    one table and read with another yields silently wrong filter
    results.

    For that reason this module keeps **no** fallback table of its own.
    An earlier revision carried one (ordered by atomic number, 56
    elements) and it was a liability rather than a convenience: a
    library built before :mod:`cif_meta` existed would have been
    misread forever after, with no error raised anywhere.  Failing
    loudly is the only safe behaviour.

    Returns
    -------
    _ElementBits
        The mask function, the bit table, its version tag and its width.

    Raises
    ------
    ImportError
        When :mod:`cif_meta` is unavailable or does not expose the
        element bit table.
    """
    try:
        from . import cif_meta
    except ImportError as exc:  # pragma: no cover - depends on deployment
        raise ImportError(
            "xinmatch.core.db.cif_meta provides the element bit table and is "
            "required by the v3 writer: the bit layout is baked into the "
            "stored data, so the writer and the reader must share one "
            "table.  See xinmatch/core/db/cif_meta.py."
        ) from exc

    mask_fn = getattr(cif_meta, "element_mask", None)
    positions = getattr(cif_meta, "element_bit_positions", None)
    if mask_fn is None or positions is None:  # pragma: no cover - defensive
        raise ImportError(
            "xinmatch.core.db.cif_meta must expose 'element_mask' and "
            "'element_bit_positions'; the v3 writer stores the resulting "
            "bits in entries.elem_mask and cannot build it on its own."
        )
    table = dict(positions)
    scheme = str(getattr(cif_meta, "ELEMENT_MASK_SCHEME", "unknown"))
    full_bits = int(getattr(cif_meta, "ELEMENT_MASK_BITS", 0) or 0)
    if full_bits <= 0:
        full_bits = (max(table.values()) + 1) if table else 0
    return _ElementBits(mask_fn, table, scheme, full_bits)


def _split_mask(raw: int | tuple[int, int]) -> tuple[int, int]:
    """Split a raw element mask into its low and high 64-bit halves.

    :mod:`cif_meta` encodes a phase's elements as a 128-bit
    ``(low, high)`` pair; a 64-bit implementation returns a single
    integer.  Both are normalised here so the writer does not have to
    care which scheme the bit table uses.

    Parameters
    ----------
    raw:
        The value returned by the mask function.

    Returns
    -------
    tuple[int, int]
        ``(low, high)``; ``high`` is ``0`` for a 64-bit scheme.
    """
    if isinstance(raw, tuple):
        low, high = raw
        return int(low), int(high)
    value = int(raw)
    return value & 0xFFFFFFFFFFFFFFFF, value >> 64


def _normalize_elements(raw: Any, known_bits: dict[str, int]) -> list[str]:
    """Normalise the ``elements`` field of a phase into a list of symbols.

    Accepts either a sequence of symbols or a delimited string (``","``
    or ``";"``), because both shapes occur in CIF-derived metadata.
    Symbols are canonicalised to title case (``"SI"`` -> ``"Si"``) so the
    stored text, the bitmask and the reader's element filter all speak
    the same language, and duplicates are dropped so a phase contributes
    exactly one bit and one inverted-index entry per element.

    Parameters
    ----------
    raw:
        The value taken from the phase metadata.
    known_bits:
        The element bit table (element symbol -> bit position).  Symbols
        outside it are **dropped**: an element that has no bit cannot be
        filtered on, and keeping it in the stored text would make
        ``entries.elements`` disagree with ``entries.elem_mask``.  The
        whitelist also keeps formula-parsing noise (a base64-encoded
        formula is a classic source of bogus "elements" such as ``"Q"``
        or ``"Zz"``) out of the database entirely.

    Returns
    -------
    list[str]
        Cleaned, sorted element symbols that have a bit.
    """
    if raw is None:
        return []
    if isinstance(raw, str):
        items: list[Any] = raw.replace(";", ",").split(",")
    else:
        items = list(raw)
    out: set[str] = set()
    for item in items:
        symbol = str(item).strip().title()
        if symbol and symbol in known_bits:
            out.add(symbol)
    return sorted(out)


_CELL_NAMES = ("a", "b", "c", "alpha", "beta", "gamma")


def _cell_values(meta: Any) -> tuple[float | None, float | None, float | None,
                                     float | None, float | None, float | None]:
    """Extract ``(a, b, c, alpha, beta, gamma)`` from phase metadata.

    Three shapes of the cell are accepted, because the metadata
    dataclass has changed shape over the life of the project and the
    writer must not silently drop the cell parameters:

    * flat attributes on the metadata object (``meta.a``, ``meta.b`` …),
    * a nested cell *object* with those attributes (``meta.cell.a`` …),
    * a nested cell *sequence* of six numbers in the canonical order
      ``(a, b, c, alpha, beta, gamma)`` — this is what
      :class:`xinmatch.core.db.cif_meta.PhaseMeta` carries.

    Parameters
    ----------
    meta:
        Phase metadata object.

    Returns
    -------
    tuple
        The six cell values, ``None`` where unknown.
    """
    cell = getattr(meta, "cell", None)
    values: list[float | None] = []
    for name in _CELL_NAMES:
        value = getattr(meta, name, None)
        if value is None and cell is not None:
            value = getattr(cell, name, None)
        values.append(None if value is None else float(value))

    if cell is not None and any(value is None for value in values):
        # A plain sequence carries no attributes, so it is read
        # positionally.  A cell object that is not iterable raises here
        # and is simply ignored.
        try:
            sequence = [
                None if item is None else float(item) for item in cell
            ]
        except (TypeError, ValueError):
            sequence = []
        if len(sequence) == len(_CELL_NAMES):
            values = [
                value if value is not None else fallback
                for value, fallback in zip(values, sequence, strict=False)
            ]

    return (
        values[0], values[1], values[2], values[3], values[4], values[5]
    )


def _cell_volume(
    a: float | None,
    b: float | None,
    c: float | None,
    alpha: float | None,
    beta: float | None,
    gamma: float | None,
) -> float | None:
    """Compute the unit-cell volume (Å³) from the cell parameters.

    Uses ``V = abc * sqrt(1 - cos²α - cos²β - cos²γ + 2 cosα cosβ cosγ)``.

    Parameters
    ----------
    a, b, c:
        Cell edge lengths in Å.
    alpha, beta, gamma:
        Cell angles in degrees.

    Returns
    -------
    float | None
        The volume, or ``None`` when the cell is incomplete or degenerate.
    """
    if None in (a, b, c, alpha, beta, gamma):
        return None
    assert a is not None and b is not None and c is not None
    assert alpha is not None and beta is not None and gamma is not None
    ca = np.cos(np.radians(alpha))
    cb = np.cos(np.radians(beta))
    cg = np.cos(np.radians(gamma))
    radicand = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if radicand <= 0.0:
        return None
    return float(a * b * c * np.sqrt(radicand))


class XinDBWriter:
    """Incrementally write phases into a ``.xinDB`` (v3 columnar) file.

    The writer is *streaming*: phases are handed over one at a time with
    :meth:`add` and are packed into chunk rows as soon as a chunk is
    full, so the peak data of a full build never has to be held in
    memory.  :meth:`close` flushes the last chunk, materialises the
    element inverted index and writes the ``meta`` table.

    Args:
        path: Target ``.xinDB`` file.  Parent directories are created.
        top_k: Number of strongest reflections per phase kept in the
            prescreen matrix.  Must match the ``top_k`` recorded in
            ``meta``, because the reader uses it to decode the block.
        chunk_size: Number of phases per ``peak_chunks`` row.
        tier: Default tier written to the ``entries`` table
            (``0`` = experimentally determined, ``1`` = predicted).
        cif_source: Provenance of the input CIFs, stored in ``meta``.
        library_kind: Kind of library being built, stored in ``meta``.
        extra_meta: Additional ``meta`` key/value pairs.
        build_time: ISO timestamp written to ``meta``; defaults to now.

    Raises
    ------
    ImportError
        When :mod:`xinmatch.core.db.cif_meta` is unavailable.  The element
        bit table lives there because the bit layout is baked into the
        stored data, so the writer and the reader must share one table.

    Example
    -------
    >>> with XinDBWriter("COD-Inorg.xinDB") as writer:      # doctest: +SKIP
    ...     for meta, d, i in phases:
    ...         writer.add(meta, d, i)
    """

    def __init__(
        self,
        path: str | Path,
        *,
        top_k: int = TOP_K,
        chunk_size: int = CHUNK_SIZE,
        tier: int = 1,
        cif_source: str = "",
        library_kind: str = "unknown",
        build_tool: str = "tools.build_cod_xindb",
        data_date: str = "",
        extra_meta: dict[str, object] | None = None,
        build_time: str | None = None,
        resume: bool = False,
    ) -> None:
        if int(top_k) < 1:
            raise ValueError("top_k must be >= 1")
        if int(chunk_size) < 1:
            raise ValueError("chunk_size must be >= 1")

        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.top_k = int(top_k)
        self.chunk_size = int(chunk_size)
        self.tier = int(tier)
        self.cif_source = str(cif_source)
        self.library_kind = str(library_kind)
        self.build_tool = str(build_tool)
        self.data_date = str(data_date)
        self.extra_meta = dict(extra_meta or {})
        self.build_time = build_time or datetime.now().astimezone().isoformat(
            timespec="seconds"
        )

        # The element bit table comes from cif_meta alone: the bit layout
        # is baked into the stored data, so the writer and the reader must
        # share one table (see _resolve_element_mask).
        self._bits = _resolve_element_mask()
        self._element_bits = self._bits.positions
        # Elements whose bit does not fit into the single ``INTEGER``
        # column; they stay in ``entries.elements`` and in
        # ``element_inv`` but cannot be expressed as a bit.
        self._unmapped_elements = sorted(
            symbol
            for symbol, bit in self._element_bits.items()
            if bit >= _STORABLE_MASK_BITS
        )
        self._mask_truncated = 0

        self.conn = sqlite3.connect(str(self.path))
        self._apply_pragmas()
        create_tables(self.conn)

        # -- entries bookkeeping -----------------------------------------
        # Initialised *before* the guard / resume step below, because resuming
        # loads the counters that live here.
        self._entry_batch: list[tuple[Any, ...]] = []
        self._seen_ids: set[str] = set()
        self._entry_count = 0
        self._next_index = 0
        self._skipped = 0

        # -- chunk bookkeeping -------------------------------------------
        self._chunk_id = 0
        self._chunk_count = 0
        self._global_peak = 0
        self._top_d_buf = np.zeros((self.chunk_size, self.top_k), dtype="<f4")
        self._top_i_buf = np.zeros((self.chunk_size, self.top_k), dtype="<f4")
        self._chunk_all_d: list[np.ndarray] = []
        self._chunk_all_i: list[np.ndarray] = []

        #: Number of phases inherited from an earlier run (0 for a fresh build).
        self._resumed_from = 0
        #: Entry id of the first inherited phase that has to be recomputed;
        #: ``None`` means the build simply appends after the last stored phase.
        self._resume_from_entry: str | None = None
        if resume:
            self._load_resume_state()
        else:
            self._guard_existing_database()

        self._closed = False
        self._stats: dict[str, int] = {}

    # ── setup helpers ───────────────────────────────────────────
    def _apply_pragmas(self) -> None:
        """Tune the connection for a write-once build.

        The database is written once and then only ever read, so the
        durability machinery is pure overhead: skipping the journal and
        the fsyncs roughly halves the build time.  The trade-off is that
        an interrupted build cannot be rolled back — it has to be redone
        from the start (or resumed by the build script, which tracks its
        own progress).
        """
        self.conn.execute("PRAGMA journal_mode = OFF")
        self.conn.execute("PRAGMA synchronous = OFF")
        self.conn.execute("PRAGMA temp_store = MEMORY")
        self.conn.execute("PRAGMA cache_size = -131072")

    def _guard_existing_database(self) -> None:
        """Refuse to write into a database that already holds phases.

        Appending to a finished file is only safe when the previous run can be
        wound back to a chunk boundary, which is exactly what
        :meth:`_load_resume_state` does for ``resume=True``.  Doing it blindly
        would let the integer entry indices and the reflection offsets silently
        stop matching, so a plain (non-resuming) open still fails loudly.
        """
        try:
            existing = int(
                self.conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
            )
        except sqlite3.Error:
            existing = 0
        if existing:
            raise RuntimeError(
                f"{self.path} already contains {existing} entries; "
                "refusing to append to a non-empty .xinDB "
                "(pass resume=True or build a fresh file)"
            )

    def _load_resume_state(self) -> None:
        """Adopt the state of an interrupted build so it can be continued.

        Three things make this safe where a naive append would corrupt the file:

        * **The two tables are flushed independently.**  ``entries`` rows are
          committed every ``_ENTRY_BATCH`` phases, ``peak_chunks`` once per
          ``chunk_size`` phases, and neither commit waits for the other.  A kill
          can therefore leave *entries ahead* (rows whose peak payload never
          landed -- the common case) or, in a narrow window right after a chunk
          commit, *chunks ahead* (a chunk whose tail entries were still
          buffered).  Both are wound back to the last boundary the two tables
          agree on, and only whole chunks are kept.
        * **Geometry cannot be read from ``meta``.**  ``meta`` is written by
          :meth:`close`, so an interrupted file has none at all.  The stored
          chunks are the only witness of how they were laid out, so
          ``chunk_size`` and ``top_k`` are verified against their byte lengths
          instead.  A resumed run that disagrees is refused rather than allowed
          to append rows the reader would mis-slice.
        * **``element_inv``.**  It is now rebuilt from ``entries`` at
          :meth:`close`, so the inverted index of the inherited phases comes
          back for free instead of being lost with the dead process.

        Raises:
            RuntimeError: The file is in a shape that cannot be continued
                (non-contiguous chunks, a mismatched ``chunk_size``/``top_k``,
                or an offset chain that no longer adds up).
        """
        try:
            entries_n = int(
                self.conn.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
            )
            rows = self.conn.execute(
                "SELECT chunk_id, n, length(top_d), length(all_d) "
                "FROM peak_chunks ORDER BY chunk_id"
            ).fetchall()
            closed = self.conn.execute(
                "SELECT 1 FROM meta WHERE key = 'schema_version'"
            ).fetchone()
        except sqlite3.Error as exc:
            raise RuntimeError(f"cannot read resume state from {self.path}: {exc}") from exc

        if closed is not None:
            # ``meta`` is written by close() alone, so its presence means the
            # previous run finished.  Appending to a complete library is not
            # what --resume is for, and its last chunk may be partial.
            raise RuntimeError(
                f"{self.path} is a complete library ({entries_n:,} phases, "
                "schema_version recorded); there is nothing to resume"
            )

        if not entries_n and not rows:
            return  # nothing was written before the interruption

        # -- geometry, verified against the payload instead of meta --------
        for expected, (chunk_id, n, top_bytes, all_bytes) in enumerate(rows):
            if int(chunk_id) != expected:
                raise RuntimeError(
                    f"{self.path}: peak_chunks are not contiguous (chunk_id "
                    f"{chunk_id} where {expected} was expected); refusing to resume"
                )
            if int(n) != self.chunk_size:
                raise RuntimeError(
                    f"{self.path}: chunk {chunk_id} holds {n} phases but the writer "
                    f"was opened with chunk_size={self.chunk_size}; resume with the "
                    "original geometry"
                )
            want_top = int(n) * self.top_k * 4
            if int(top_bytes) != want_top:
                raise RuntimeError(
                    f"{self.path}: chunk {chunk_id} stores {int(top_bytes) // 4} "
                    f"d-spacings for {n} phases ({want_top // 4} expected for "
                    f"top_k={self.top_k}); resume with the original geometry"
                )
            if int(all_bytes) < 4:
                raise RuntimeError(
                    f"{self.path}: chunk {chunk_id} has an empty reflection payload; "
                    "refusing to resume"
                )
        covered = len(rows) * self.chunk_size

        # -- wind both tables back to the last shared chunk boundary -------
        target = (min(entries_n, covered) // self.chunk_size) * self.chunk_size
        dropped_entries = entries_n - target
        dropped_chunks = len(rows) - target // self.chunk_size
        if dropped_entries or dropped_chunks:
            if dropped_entries:
                # The first discarded row is where the CIF walk has to restart,
                # so remember it before it is deleted.
                row = self.conn.execute(
                    "SELECT entry_id FROM entries WHERE rowid = ?", (target,)
                ).fetchone()
                self._resume_from_entry = str(row[0]) if row else None
                self.conn.execute("DELETE FROM entries WHERE rowid >= ?", (target,))
            if dropped_chunks:
                self.conn.execute(
                    "DELETE FROM peak_chunks WHERE chunk_id >= ?",
                    (target // self.chunk_size,),
                )
            self.conn.commit()
            logger.warning(
                "resume: wound %s back to %d phases (dropped %d trailing entries "
                "and %d partial chunk(s))",
                self.path.name,
                target,
                dropped_entries,
                dropped_chunks,
            )

        # -- the offset chain must still add up ----------------------------
        peaks_entries = int(
            self.conn.execute("SELECT COALESCE(SUM(peak_n), 0) FROM entries").fetchone()[0]
        )
        peaks_chunks = int(
            self.conn.execute(
                "SELECT COALESCE(SUM(length(all_d)), 0) FROM peak_chunks"
            ).fetchone()[0]
        ) // 4
        if peaks_entries != peaks_chunks:
            raise RuntimeError(
                f"{self.path}: entries account for {peaks_entries:,} reflections but "
                f"the chunks hold {peaks_chunks:,}; the offset chain is broken and "
                "cannot be resumed"
            )

        # A *closed* library records the geometry it was built with; when the
        # resumed run uses a different one the old chunks no longer mean the
        # same thing.
        for key, mine in (("top_k", self.top_k), ("chunk_size", self.chunk_size)):
            stored = self.conn.execute(
                "SELECT value FROM meta WHERE key = ?", (key,)
            ).fetchone()
            if stored is not None and str(stored[0]) != str(mine):
                raise RuntimeError(
                    f"{self.path}: stored {key}={stored[0]!r} but the writer was "
                    f"opened with {key}={mine!r}; resume with the original geometry"
                )

        self._global_peak = peaks_entries
        self._entry_count = target
        self._next_index = target
        self._resumed_from = target

        # ``mask_truncated`` counts the phases whose element list does not fit
        # into the 63-bit ``elem_mask`` column.  The criterion is "carries an
        # element whose bit is >= 63", which is exactly the unmapped set, so the
        # inherited population can be recounted from the stored element strings
        # instead of being lost with the dead process.
        if self._unmapped_elements:
            unmapped = set(self._unmapped_elements)
            self._mask_truncated = sum(
                1
                for (elems,) in self.conn.execute("SELECT elements FROM entries")
                if elems and not unmapped.isdisjoint(str(elems).split(","))
            )

        self._chunk_id = target // self.chunk_size
        self._chunk_count = 0

        # Duplicate detection has to see the inherited ids too, otherwise a CIF
        # that was already stored would be added a second time.
        self._seen_ids = {
            str(r[0]) for r in self.conn.execute("SELECT entry_id FROM entries")
        }

        # ``element_inv`` may still hold the previous run's rows when resuming a
        # file that had been closed; the close-time rebuild replaces all of it.
        with contextlib.suppress(sqlite3.Error):
            self.conn.execute("DELETE FROM element_inv")
            self.conn.commit()

    # ── public API ───────────────────────────────────────────────
    @property
    def resumed_from(self) -> int:
        """Phases inherited from an earlier, interrupted run (0 for a fresh build)."""
        return self._resumed_from

    @property
    def resume_from_entry_id(self) -> str | None:
        """First inherited phase that must be recomputed on resume.

        ``None`` means the file ends exactly on a chunk boundary with every
        stored phase intact, so the build can simply continue after the last
        one.  A non-``None`` value is the entry id of the earliest discarded
        row: the caller has to restart its CIF walk *at* that phase, not after
        it, otherwise the dropped phases are lost for good.
        """
        return self._resume_from_entry

    def add(
        self,
        meta: PhaseMeta,
        d: np.ndarray,
        i: np.ndarray,
    ) -> None:
        """Append one phase to the database.

        Parameters
        ----------
        meta:
            Phase metadata (formula, cell, space group, elements, ...).
            Any object exposing those attributes works; the canonical
            shape is :class:`xinmatch.core.db.cif_meta.PhaseMeta`.
        d:
            d-spacings in Å, float64.  Stored in the order supplied; the
            peak engine emits descending d, which is what readers assume
            for display.  Scoring does not depend on this order (see
            :mod:`xinmatch.core.search.correlate`).
        i:
            Intensities parallel to *d*, float64.  Rescaled to the
            1000-scale inside this method, so raw values are accepted.

        Raises
        ------
        ValueError
            When *d* and *i* have different lengths, or when the phase
            carries no usable ``entry_id``.
        RuntimeError
            When the writer has already been closed.
        """
        if self._closed:
            raise RuntimeError("writer is closed")
        d_arr: np.ndarray = np.asarray(d, dtype=np.float64).ravel()
        i_arr: np.ndarray = np.asarray(i, dtype=np.float64).ravel()
        if d_arr.size != i_arr.size:
            raise ValueError(
                f"d and i must have the same length ({d_arr.size} != {i_arr.size})"
            )

        entry_id = str(getattr(meta, "entry_id", "") or "").strip()
        if not entry_id:
            raise ValueError("phase metadata carries no entry_id")
        if self._is_duplicate(entry_id):
            # A second row for the same phase would consume an entry
            # index that no ``entries`` row owns and would desynchronise
            # the rowid <-> index mapping the inverted index relies on.
            self._skipped += 1
            logger.debug("skipping duplicate entry_id %s", entry_id)
            return

        # The file declares ``peak_convention = "norm-1000"`` in its
        # meta table, so the writer guarantees that scale rather than
        # trusting the caller.  The upstream peak calculator already
        # normalises, which makes this a no-op there (the factor is
        # exactly 1.0); for a caller that hands over raw intensities it
        # is what keeps the library on one shared scale.
        if i_arr.size:
            peak = float(np.max(i_arr))
            if peak > 0.0:
                i_arr = i_arr * (INTENSITY_SCALE / peak)

        # The peak table keeps the caller's order (for this project: the peak
        # engine's descending-d order, matching the v2 library).  The prescreen
        # matrix however needs the *strongest* reflections, so its row is filled
        # from an explicit intensity ranking instead of assuming the input is
        # already intensity-sorted -- ``search.correlate`` sorts internally, so
        # the stored order has no effect on scoring.
        top_idx = np.argsort(-i_arr, kind="stable")[: self.top_k] if i_arr.size else None

        n_peaks = int(d_arr.size)
        index = self._next_index
        self._next_index += 1

        # The element list is normalised once and shared by the entries
        # row (bitmask + display string) and the inverted index, so the
        # two can never disagree about a phase's composition.  The
        # inverted index itself is rebuilt from these rows at close time.
        elements = _normalize_elements(
            getattr(meta, "elements", None), self._element_bits
        )
        self._append_entry_row(meta, entry_id, elements, n_peaks, index)
        self._append_chunk_data(d_arr, i_arr, top_idx)

        self._entry_count += 1
        self._global_peak += n_peaks

    def close(self) -> dict[str, int]:
        """Finalise the database and return build statistics.

        Flushes the last (partial) chunk, writes the element inverted
        index, records the ``meta`` table and closes the connection.

        Returns
        -------
        dict[str, int]
            Statistics with the keys ``entries``, ``chunks``, ``peaks``,
            ``elements``, ``skipped``, ``mask_truncated`` and ``bytes``.

        Raises
        ------
        RuntimeError
            When a write error occurs; the connection is closed either
            way.
        """
        if self._closed:
            return dict(self._stats)
        n_elements = 0
        try:
            self._flush_entries()
            self._flush_chunk()
            n_elements = self._write_element_inv()
            self._write_meta()
            self.conn.commit()
        finally:
            self._closed = True
            with contextlib.suppress(sqlite3.Error):
                self.conn.close()

        try:
            size = int(self.path.stat().st_size)
        except OSError:  # pragma: no cover - file always exists here
            size = 0
        self._stats = {
            "entries": self._entry_count,
            "chunks": self._chunk_id,
            "peaks": self._global_peak,
            "elements": n_elements,
            "skipped": self._skipped,
            "mask_truncated": self._mask_truncated,
            "resumed_from": self._resumed_from,
            "bytes": size,
        }
        if self._mask_truncated:
            logger.warning(
                "%d phases carry an element whose bit does not fit into "
                "entries.elem_mask (%d bits); their element list is kept in "
                "entries.elements and element_inv, but the bitmask cannot "
                "express it",
                self._mask_truncated,
                _STORABLE_MASK_BITS,
            )
        logger.info(
            "wrote %d phases / %d reflections into %d chunks (%.1f MB) -> %s",
            self._stats["entries"],
            self._stats["peaks"],
            self._stats["chunks"],
            size / 1e6,
            self.path,
        )
        return dict(self._stats)

    def __enter__(self) -> XinDBWriter:
        """Enter the context manager."""
        return self

    def __exit__(self, exc_type: object, exc: object, tb: object) -> None:
        """Close the writer when leaving the context manager.

        When the block raised, the database is closed *without* writing
        the ``meta`` table: an aborted build must not look like a
        complete one, and a missing ``schema_version`` is exactly what
        the reader's ``probe()`` rejects.
        """
        if self._closed:
            return
        if exc_type is not None:
            self._closed = True
            with contextlib.suppress(sqlite3.Error):
                self.conn.close()
            return
        self.close()

    # ── entries ──────────────────────────────────────────────────
    def _is_duplicate(self, entry_id: str) -> bool:
        """Return ``True`` when *entry_id* was already written.

        The set holds every id written during this run.  Keeping the ids
        in memory rather than querying ``entries`` keeps ``add()`` free
        of a per-phase SQL round trip, which measurably dominates the
        write cost of a full-COD build.  The set is bounded by the
        number of phases (a few tens of MB at COD scale) and the file
        itself is guarded against being appended to, so the set is an
        exact mirror of what the ``entries`` table holds.
        """
        return entry_id in self._seen_ids

    def _append_entry_row(
        self,
        meta: Any,
        entry_id: str,
        elements: list[str],
        n_peaks: int,
        index: int,
    ) -> None:
        """Build and buffer the ``entries`` row of one phase.

        The phase's ``rowid`` is written explicitly and equals the
        integer index used by the ``element_inv`` inverted index, so
        ``SELECT rowid FROM entries WHERE entry_id = ?`` returns exactly
        the number stored in the inverted index.

        The element mask is truncated to the 63 bits that fit into the
        ``INTEGER`` column (see the module docstring).  A phase with an
        element outside that window keeps the element in
        ``entries.elements`` and in ``element_inv``; only the fast
        bitmask path cannot express it, and the count of such phases is
        reported in the build statistics.
        """
        low, high = _split_mask(self._bits.mask_fn(elements))
        mask = low & _MAX_STORABLE_MASK
        if high or low > _MAX_STORABLE_MASK:
            # The phase carries at least one element whose bit does not
            # fit into the column.  The element itself is not lost (it
            # is in ``entries.elements`` and in ``element_inv``), but the
            # bitmask cannot express it, so the phase is counted.
            self._mask_truncated += 1

        a, b, c, alpha, beta, gamma = _cell_values(meta)
        volume = getattr(meta, "volume", None)
        if volume is None:
            volume = _cell_volume(a, b, c, alpha, beta, gamma)
        density = getattr(meta, "density", None)
        year = getattr(meta, "year", None)
        tier = getattr(meta, "tier", None)
        if tier is None:
            tier = self.tier

        row = (
            index,
            entry_id,
            getattr(meta, "formula", None),
            getattr(meta, "space_group", None),
            getattr(meta, "crystal_system", None),
            getattr(meta, "mineral_name", None),
            ",".join(elements),
            mask,
            a, b, c, alpha, beta, gamma,
            None if volume is None else float(volume),
            None if density is None else float(density),
            None if year is None else int(year),
            getattr(meta, "author", None),
            getattr(meta, "journal", None),
            getattr(meta, "quality", None),
            int(n_peaks),
            self._global_peak,
            int(n_peaks),
            int(tier),
        )
        self._entry_batch.append(row)
        self._seen_ids.add(entry_id)
        if len(self._entry_batch) >= _ENTRY_BATCH:
            self._flush_entries()

    def _flush_entries(self) -> None:
        """Write the buffered ``entries`` rows to SQLite."""
        if not self._entry_batch:
            return
        self.conn.executemany(
            """
            INSERT INTO entries (
                rowid, entry_id, formula, space_group, crystal_system,
                mineral_name, elements, elem_mask,
                a, b, c, alpha, beta, gamma, volume, density,
                year, author, journal, quality,
                n_peaks, peak_off, peak_n, tier
            ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
            """,
            self._entry_batch,
        )
        self.conn.commit()
        self._entry_batch.clear()

    # ── peak chunks ──────────────────────────────────────────────
    def _append_chunk_data(
        self,
        d_arr: np.ndarray,
        i_arr: np.ndarray,
        top_idx: np.ndarray | None = None,
    ) -> None:
        """Append one phase's reflections to the chunk being assembled.

        *top_idx* selects the reflections that go into the prescreen matrix (the
        strongest ones); the full CSR arrays keep the caller's order.
        """
        row = self._chunk_count
        # The buffer is reused for every chunk, so the row has to be
        # cleared first: a phase with fewer than top_k reflections would
        # otherwise inherit the reflections of the phase that occupied
        # the same row in the previous chunk, and the stale d-spacings
        # would create phantom prescreen hits.
        self._top_d_buf[row, :] = 0.0
        self._top_i_buf[row, :] = 0.0
        if top_idx is not None:
            take = int(min(self.top_k, top_idx.size))
            sel = top_idx[:take]
            self._top_d_buf[row, :take] = d_arr[sel].astype("<f4")
            self._top_i_buf[row, :take] = i_arr[sel].astype("<f4")
        # Phases with fewer than top_k reflections keep the 0.0 padding
        # written above: a padded d of 0.0 can never fall into a
        # measurement window (real d-spacings start at 0.77 A) and a
        # padded intensity of 0.0 carries no weight, so the padding
        # cannot create a false hit or skew the ranking.
        self._chunk_all_d.append(d_arr.astype("<f4"))
        self._chunk_all_i.append(i_arr.astype("<f4"))
        self._chunk_count += 1
        if self._chunk_count >= self.chunk_size:
            self._flush_chunk()

    def _flush_chunk(self) -> None:
        """Write the assembled chunk row to ``peak_chunks``."""
        if self._chunk_count == 0:
            return
        n = self._chunk_count
        all_d = (
            np.concatenate(self._chunk_all_d)
            if self._chunk_all_d
            else np.zeros(0, dtype="<f4")
        )
        all_i = (
            np.concatenate(self._chunk_all_i)
            if self._chunk_all_i
            else np.zeros(0, dtype="<f4")
        )
        self.conn.execute(
            """
            INSERT INTO peak_chunks (chunk_id, n, top_d, top_i, all_d, all_i)
            VALUES (?,?,?,?,?,?)
            """,
            (
                self._chunk_id,
                n,
                self._top_d_buf[:n].tobytes(),
                self._top_i_buf[:n].tobytes(),
                all_d.tobytes(),
                all_i.tobytes(),
            ),
        )
        self.conn.commit()
        self._chunk_id += 1
        self._chunk_count = 0
        self._chunk_all_d.clear()
        self._chunk_all_i.clear()

    # ── element inverted index ───────────────────────────────────
    def _build_element_index(self) -> dict[str, array]:
        """Rebuild the inverted index from ``entries.elements``.

        Reading the table instead of accumulating in memory is what makes
        ``resume`` possible: a continued build inherits the phases of the dead
        process, and their elements have to come back into the index.  It also
        removes the only piece of per-phase state that had to survive the whole
        run, so 400 k+ phases no longer pin a Python list per element.

        ``entries.elements`` is written as ``",".join(symbols)`` by
        :meth:`_append_entry_row`; an empty string means the phase carries no
        element and contributes nothing.

        Returns:
            Element symbol -> ascending array of entry indices.
        """
        buckets: dict[str, array] = {}
        cursor = self.conn.execute("SELECT rowid, elements FROM entries ORDER BY rowid")
        for rowid, elements in cursor:
            text = str(elements or "")
            if not text:
                continue
            for symbol in text.split(","):
                symbol = symbol.strip()
                if not symbol:
                    continue
                bucket = buckets.get(symbol)
                if bucket is None:
                    bucket = array("I")
                    buckets[symbol] = bucket
                bucket.append(int(rowid))
        return buckets

    def _write_element_inv(self) -> int:
        """Write the ``element_inv`` table from the ``entries`` table.

        Returns:
            Number of distinct elements written.
        """
        buckets = self._build_element_index()
        # A resumed file may still hold the previous run's rows; DELETE first so
        # an element that exists only in the inherited data cannot keep a stale
        # (shorter) list.
        self.conn.execute("DELETE FROM element_inv")
        if buckets:
            self.conn.executemany(
                "INSERT OR REPLACE INTO element_inv (element, eids) VALUES (?,?)",
                [
                    (element, pack_uint32(bucket))
                    for element, bucket in sorted(buckets.items())
                ],
            )
        self.conn.commit()
        return len(buckets)

    # ── meta ─────────────────────────────────────────────────────
    def _write_meta(self) -> None:
        """Write the standard ``meta`` set plus the writer-specific keys."""
        write_standard_meta(
            self.conn,
            entry_count=self._entry_count,
            cif_source=self.cif_source,
            build_time=self.build_time,
            library_kind=self.library_kind,
            build_tool=self.build_tool,
            data_date=self.data_date,
            top_k=self.top_k,
            chunk_size=self.chunk_size,
        )
        set_meta(self.conn, "tier_default", self.tier)
        # The element bit layout is *baked into* ``entries.elem_mask``,
        # so the table that produced it is recorded in the file itself.
        # A reader can therefore decode the stored masks correctly even
        # if its own table ever differs, and it can detect a mismatch
        # instead of silently returning wrong filter results.
        set_meta(self.conn, "element_mask_scheme", self._bits.scheme)
        set_meta(self.conn, "element_mask_bits", _STORABLE_MASK_BITS)
        set_meta(self.conn, "element_mask_full_bits", self._bits.full_bits)
        set_meta(self.conn, "element_bits", self._encode_element_bits())
        set_meta(
            self.conn,
            "element_mask_unmapped",
            self._encode_unmapped_elements(),
        )
        for key, value in self.extra_meta.items():
            set_meta(self.conn, key, value)
        self.conn.commit()

    def _encode_element_bits(self) -> str:
        """Serialise the element bit table into a compact meta string.

        The format is ``"Symbol:bit;Symbol:bit;..."`` with the symbols in
        bit order, so it is both human-readable and trivially parseable
        by a reader that wants to validate its own table.
        """
        ordered = sorted(self._element_bits.items(), key=lambda kv: kv[1])
        return ";".join(f"{symbol}:{bit}" for symbol, bit in ordered)

    def _encode_unmapped_elements(self) -> str:
        """Serialise the elements that have no bit in ``elem_mask``.

        These are the elements whose bit number does not fit into the
        63-bit window the ``entries.elem_mask`` column offers.  They are
        listed so a reader (or a human auditing the file) can see that
        the bitmask does not cover them — the exact element data is in
        ``entries.elements`` and in the ``element_inv`` index.
        """
        return ";".join(
            f"{symbol}:{self._element_bits[symbol]}"
            for symbol in self._unmapped_elements
        )
