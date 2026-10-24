"""Read-only acceptance checker for a built ``.xinDB`` reference library.

The build tools report what they *think* they wrote.  This module verifies what
is actually on disk, from the outside, using nothing but the public schema -- so
it can be run on a library built days ago, on a different machine, or by a
different tool (``build_cod_xindb`` / ``build_cod_inorg_xindb`` / a future
import path).

Checks are grouped and each carries an id:

``S`` -- schema
    The tables a v3 reader requires exist.
``M`` -- metadata provenance
    Every ``META_KEYS`` entry is present, ``peak_convention`` equals the frozen
    constant from :mod:`peak_compute` (not a copy), and ``library_kind`` /
    ``data_date`` agree with the file name.
``C`` -- chunk accounting
    ``SUM(peak_chunks.n)`` and ``SUM(entries.n_peaks)`` agree with the entry
    count, the chunk count matches ``ceil(entry_count / chunk_size)``, and the
    global ``peak_off`` ladders are contiguous -- i.e. no phase can point at a
    reflection outside its chunk.
``E`` -- entry sanity
    Lattice parameters and volumes are physical, ``peak_n == n_peaks``, no
    entry carries more reflections than ``top_k``, and the *reflection count*
    itself is sane:

    ``E5.boundary_peak``
        The pre-2026-10 peak calculator had two defects that survive into a
        built library, and both are caught here **exhaustively** (every phase,
        not a sample):

        * **180 deg boundary artifact** -- the strongest reflection sits at the
          d-spacing cutoff (``DEFAULT_DMIN`` = 0.7703 A = lambda/2) and the rest
          of the pattern has collapsed towards zero, because the old
          Lorentz-polarisation factor diverged at 2theta = 180 deg and the
          normalized pattern became that one line at 1000.  The collapse
          condition is what keeps the check specific -- any *real* reflection
          near the cutoff also becomes the strongest line for the same LP
          reason, so "strongest d at the cutoff" alone would flag 3-5 % of a
          healthy library.  Detected from ``peak_chunks.top_d`` column 0 (the
          strongest reflection, rows are intensity-ordered); also catches the
          collapsed two-line variant (e.g. ``96-254-30089``: bogus 180 deg line
          at 1000 plus one real line at 0.013).  **FAIL**;
        * **placeholder-cell artifact** -- a lone reflection at exactly 1.0 A,
          the placeholder unit cell gemmi invents for a CIF without
          ``_cell_*`` tags.  **FAIL**;
        * any other lone reflection is physically possible (a large cell with
          one in-window line), so it is reported as **WARN**;
        * the **semi-collapsed** residual -- strongest reflection at the cutoff
          but a still-populated pattern -- is reported as **WARN** too (up to 10
          examples + total, smallest remaining intensity first, e.g. PCOD
          ``96-244-04910``).  It shares the root cause but is not the build
          defect, so it does not block acceptance.

        Measured on the shipped libraries: 8 phases in ``COD-Inorg``
        (5 boundary + 3 placeholder), 7 in ``PCOD-Inorg`` (6 boundary singles +
        ``96-254-30089``), 0 in ``TCOD-Inorg``.  Semi-collapsed WARN: 2,251 /
        5,607 / 1,670 (COD / PCOD / TCOD).

    ``E6.peak_count``
        A phase with ``n_peaks <= 0`` (or ``NULL``) cannot be scored and the
        writer is not supposed to emit one; the acceptance pass still has to
        catch it.  **FAIL**.
``P`` -- peak payload
    Sampled chunks (or every chunk with ``--all-chunks``) decode to the
    documented shapes; every stored reflection is finite and positive;
    intensities stay on the ``norm-1000`` scale; the zero-padding after
    ``n_peaks`` is actually zero; and the prescreen ``top`` block is a subset of
    the phase's full reflection list.
``X`` -- element inverted index
    Sampled phases find themselves in ``element_inv`` under every symbol listed
    in ``entries.elements``, and every ``eid`` is a valid entry index.

A ``FAIL`` makes the exit code 1; ``WARN`` and ``INFO`` do not.

Usage::

    python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB
    python -m tools.validate_xindb data/COD-Full-20261020.xinDB --json out.json
    python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --all-chunks
"""

from __future__ import annotations

import argparse
import json
import random
import sqlite3
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from xinmatch.core.db.peak_compute import (  # noqa: E402
    DEFAULT_DMIN,
    DEFAULT_MAX_PEAKS,
    INTENSITY_SCALE,
    PEAK_CONVENTION,
)
from xinmatch.core.db.schema_v3 import META_KEYS, REQUIRED_TABLES  # noqa: E402

#: Intensity scale the ``norm-1000`` convention promises.  A phase's strongest
#: reflection must land exactly here, otherwise the prescreen scores are not
#: comparable across libraries.
NORM_SCALE = float(INTENSITY_SCALE)

#: Slack for the normalisation check: the writer rescales with a float32
#: division, so the maximum can be a few ULP off 1000.
NORM_TOL = 0.05

#: Lower slack of the ``E5`` boundary window.  ``0.7703`` is stored as float32
#: ``0.77029997``, so the artifact never lands exactly on ``DEFAULT_DMIN``.
DMIN_LOW_SLACK = 1e-6

#: Upper slack of the ``E5`` boundary window.  The bogus line is an ulp *above*
#: ``lambda/2`` (``0.7703000000000001`` before float32 rounding); a few 1e-3 A
#: is generous but still far below the ~0.01 A matching window, so a genuine
#: reflection cannot be swallowed.
DMIN_HIGH_SLACK = 1e-3

#: A *collapsed* 180 deg artifact puts essentially the whole normalised
#: intensity into the bogus boundary line: everything else sums to less than
#: this on the ``norm-1000`` scale (0.5 %).  This second condition is what makes
#: the strongest-reflection scan usable -- every *real* reflection that happens
#: to fall near 2theta = 180 deg also becomes the strongest one, because the
#: Lorentz-polarisation factor diverges there, so "strongest d is at the cutoff"
#: alone would flag 3-5 % of a healthy library (measured: 2,256 COD / 5,614 PCOD
#: / 1,670 TCOD phases).  Only a pattern whose intensity has actually collapsed
#: onto that line is the pre-2026-10 engine bug.
BOUNDARY_COLLAPSE_MAX = 5.0

#: d-spacing of the placeholder unit cell gemmi invents for a CIF with no
#: ``_cell_*`` tags (a = b = c = 1.0 A).  The old engine emitted its single
#: (1, 1, 1) line at exactly this value.
PLACEHOLDER_D = 1.0

#: Tolerance of the placeholder-cell comparison (float32 rounding of 1.0 is
#: exact, so this is only there to absorb a sloppy legacy writer).
PLACEHOLDER_TOL = 1e-3

#: How many offending entry ids a check detail line quotes (the acceptance
#: report has to be able to drive a targeted rebuild, so a small set of bad
#: phases is listed in full; only a systemic failure is truncated).
MAX_QUOTED_IDS = 500

#: How many of the *semi-collapsed* boundary phases (strongest reflection at the
#: cutoff, but a still-populated pattern) a WARN line quotes, plus the total.
MAX_BOUNDARY_WARN_QUOTED = 10

#: Lattice parameters outside these bounds are not physical for a crystal
#: structure; COD does contain a handful of garbage cells and the build is
#: supposed to have dropped them.
MIN_LATTICE = 0.5
MAX_LATTICE = 2000.0


@dataclass
class Report:
    """Collected check results."""

    rows: list[tuple[str, str, str]] = field(default_factory=list)
    fails: int = 0
    warns: int = 0

    def add(self, code: str, status: str, detail: str) -> None:
        self.rows.append((code, status, detail))
        if status == "FAIL":
            self.fails += 1
        elif status == "WARN":
            self.warns += 1

    def ok(self, code: str, detail: str = "ok") -> None:
        self.add(code, "PASS", detail)

    def fail(self, code: str, detail: str) -> None:
        self.add(code, "FAIL", detail)

    def warn(self, code: str, detail: str) -> None:
        self.add(code, "WARN", detail)

    def info(self, code: str, detail: str) -> None:
        self.add(code, "INFO", detail)

    def render(self) -> str:
        width = max((len(r[0]) for r in self.rows), default=4)
        lines = []
        for code, status, detail in self.rows:
            lines.append(f"  [{status}] {code:<{width}}  {detail}")
        return "\n".join(lines)


def _decode_matrix(blob: object, n: int, k: int) -> np.ndarray:
    data = b"" if blob is None else bytes(blob)
    expected = n * k * 4
    if len(data) < expected:
        data = data + b"\x00" * (expected - len(data))
    return np.frombuffer(data[:expected], dtype="<f4").reshape(n, k)


def _decode_flat(blob: object) -> np.ndarray:
    data = b"" if blob is None else bytes(blob)
    count = len(data) // 4
    return np.frombuffer(data[: count * 4], dtype="<f4")


def _meta(conn: sqlite3.Connection) -> dict[str, str]:
    try:
        return {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
    except sqlite3.Error:
        return {}


def _check_schema(conn: sqlite3.Connection, rep: Report) -> bool:
    names = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    missing = [t for t in REQUIRED_TABLES if t not in names]
    if missing:
        rep.fail("S1.tables", f"missing required tables: {', '.join(missing)}")
        return False
    rep.ok("S1.tables", f"{len(REQUIRED_TABLES)} required tables present")

    indexes = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='index'")}
    wanted = {
        "idx_entries_space_group",
        "idx_entries_crystal_system",
        "idx_entries_formula",
        "idx_entries_tier",
    }
    lost = sorted(wanted - indexes)
    if lost:
        rep.warn("S2.indexes", f"missing lookup indexes: {', '.join(lost)}")
    else:
        rep.ok("S2.indexes", f"{len(wanted)} lookup indexes present")
    return True


def _check_meta(conn: sqlite3.Connection, rep: Report, db: Path) -> dict[str, str]:
    meta = _meta(conn)
    blanks = [k for k in META_KEYS if k not in meta or not str(meta[k]).strip()]
    # ``data_date`` is validated separately: it is optional for libraries built
    # before the snapshot-date convention existed.
    blanks = [k for k in blanks if k != "data_date"]
    if blanks:
        rep.fail("M1.keys", f"empty or missing meta keys: {', '.join(blanks)}")
    else:
        rep.ok("M1.keys", f"{len(META_KEYS) - 1} required keys populated")

    stored = meta.get("peak_convention", "")
    if stored == PEAK_CONVENTION:
        rep.ok("M2.convention", stored)
    else:
        rep.fail(
            "M2.convention",
            f"stored {stored!r} != engine {PEAK_CONVENTION!r} (peaks were computed "
            "under a different cutoff than the file claims)",
        )

    # The file name carries the same identity, and ``multi_db`` picks the newest
    # snapshot by parsing it.  If the two disagree, the library will be loaded
    # under the wrong date.
    stem = db.stem
    kind = meta.get("library_kind", "")
    date = meta.get("data_date", "")
    head, _, tail = stem.rpartition("-")
    if head and head != kind:
        rep.warn("M3.identity", f"file name prefix {head!r} != meta.library_kind {kind!r}")
    elif date and tail != date:
        rep.fail(
            "M3.identity",
            f"file name says {tail!r} but meta.data_date says {date!r}",
        )
    elif date:
        rep.ok("M3.identity", f"kind={kind} date={date} (name and meta agree)")
    else:
        rep.warn("M3.identity", f"kind={kind} has no meta.data_date")

    rep.info("M4.provenance", f"build_time={meta.get('build_time', '?')}")
    rep.info("M5.provenance", f"cif_source={meta.get('cif_source', '?')}")
    rep.info("M6.provenance", f"build_tool={meta.get('build_tool', '?')}")
    return meta


def _check_counts(conn: sqlite3.Connection, rep: Report, meta: dict[str, str]) -> tuple[int, int, int]:
    entries = int(conn.execute("SELECT count(*) FROM entries").fetchone()[0])
    declared = meta.get("entry_count", "")
    try:
        declared_n = int(declared)
    except (TypeError, ValueError):
        declared_n = -1
    if declared_n == entries:
        rep.ok("C1.entry_count", f"meta.entry_count = {entries:,}")
    else:
        rep.fail(
            "C1.entry_count",
            f"meta.entry_count = {declared_n:,} but entries table has {entries:,}",
        )

    chunks, chunk_n = conn.execute(
        "SELECT count(*), coalesce(sum(n), 0) FROM peak_chunks"
    ).fetchone()
    if int(chunk_n) == entries:
        rep.ok("C2.chunk_rows", f"{int(chunks):,} chunks cover {entries:,} entries")
    else:
        rep.fail(
            "C2.chunk_rows",
            f"SUM(peak_chunks.n) = {int(chunk_n):,} but entries = {entries:,}",
        )

    try:
        chunk_size = int(meta.get("chunk_size", "0"))
    except ValueError:
        chunk_size = 0
    if chunk_size > 0:
        expected_chunks = -(-entries // chunk_size)
        if int(chunks) == expected_chunks:
            rep.ok("C3.chunk_size", f"ceil({entries:,}/{chunk_size}) = {int(chunks):,}")
        else:
            rep.warn(
                "C3.chunk_size",
                f"{int(chunks):,} chunks, expected {expected_chunks:,} "
                f"for chunk_size={chunk_size}",
            )
    else:
        rep.warn("C3.chunk_size", f"meta.chunk_size unreadable: {meta.get('chunk_size')!r}")

    row = conn.execute(
        "SELECT coalesce(sum(n_peaks), 0), coalesce(max(n_peaks), 0), "
        "coalesce(min(n_peaks), 0), "
        "coalesce(sum(CASE WHEN peak_n <> n_peaks THEN 1 ELSE 0 END), 0) "
        "FROM entries"
    ).fetchone()
    total_peaks, max_peaks, min_peaks, mismatched = (int(x) for x in row)

    total_all = 0
    for (blob,) in conn.execute("SELECT all_d FROM peak_chunks"):
        total_all += len(blob or b"") // 4
    if total_all == total_peaks:
        rep.ok("C4.peak_total", f"{total_peaks:,} reflections (chunks agree with entries)")
    else:
        rep.fail(
            "C4.peak_total",
            f"entries declare {total_peaks:,} reflections, chunk payload holds {total_all:,}",
        )

    if mismatched == 0:
        rep.ok("C5.peak_n", "peak_n == n_peaks for every entry")
    else:
        rep.fail("C5.peak_n", f"{mismatched:,} entries have peak_n != n_peaks")

    # ``peak_off`` is a global reflection index; the ladder must close exactly,
    # otherwise a phase's slice lands in the wrong chunk.
    tail = conn.execute(
        "SELECT peak_off + peak_n FROM entries ORDER BY rowid DESC LIMIT 1"
    ).fetchone()
    if tail is not None and int(tail[0]) == total_peaks:
        rep.ok("C6.peak_off", f"last phase ends at reflection {total_peaks:,}")
    else:
        rep.fail(
            "C6.peak_off",
            f"last phase ends at {int(tail[0]) if tail else -1:,}, expected {total_peaks:,}",
        )

    rep.info(
        "C7.peaks_per_phase",
        f"min={min_peaks} max={max_peaks} mean={total_peaks / entries:.1f}" if entries else "n/a",
    )
    return entries, total_peaks, max_peaks


def _check_entries(conn: sqlite3.Connection, rep: Report, meta: dict[str, str], max_peaks: int) -> None:
    # ``a <= 0`` is NULL-safe only if the NULL is tested first: in SQLite a
    # comparison against NULL yields NULL, which ``CASE`` treats as false, so an
    # entry with an entirely NULL cell would silently pass a plain ``a <= 0``
    # test.  Every predicate therefore spells out ``IS NULL OR``.
    bounds = (MAX_LATTICE,) * 3
    row = conn.execute(
        "SELECT "
        " coalesce(sum(CASE WHEN a IS NULL OR a <= 0 OR a > ?  THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN b IS NULL OR b <= 0 OR b > ?  THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN c IS NULL OR c <= 0 OR c > ?  THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN volume IS NULL OR volume <= 0 THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN alpha IS NULL OR alpha < 1 OR alpha > 179 THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN beta  IS NULL OR beta  < 1 OR beta  > 179 THEN 1 ELSE 0 END), 0),"
        " coalesce(sum(CASE WHEN gamma IS NULL OR gamma < 1 OR gamma > 179 THEN 1 ELSE 0 END), 0)"
        " FROM entries",
        bounds,
    ).fetchone()
    labels = ("a", "b", "c", "volume", "alpha", "beta", "gamma")
    bad = {lbl: int(v) for lbl, v in zip(labels, row, strict=True) if int(v)}
    total = int(conn.execute("SELECT count(*) FROM entries").fetchone()[0] or 0)
    if bad:
        detail = ", ".join(f"{k}={v:,}" for k, v in bad.items())
        samples = [
            r[0]
            for r in conn.execute(
                "SELECT entry_id FROM entries WHERE a IS NULL OR a <= 0 OR volume IS NULL "
                "OR volume <= 0 OR alpha IS NULL OR alpha < 1 LIMIT 5"
            )
        ]
        worst = max(bad.values())
        share = worst / total if total else 1.0
        msg = f"non-physical cell values: {detail}"
        if samples:
            msg += f" (e.g. {', '.join(samples)})"
        # A handful of degenerate CIFs is a known COD data blemish, not a build
        # bug: the build keeps them because they still carry a usable cell.  It
        # only becomes a FAIL when it looks systemic.
        if share > 0.001:
            rep.fail("E1.lattice", msg)
        else:
            rep.warn("E1.lattice", msg)
    else:
        rep.ok("E1.lattice", "all lattice parameters and volumes physical")

    empty = int(
        conn.execute(
            "SELECT count(*) FROM entries WHERE elements IS NULL OR trim(elements) = ''"
        ).fetchone()[0]
    )
    if empty:
        samples = [
            r[0]
            for r in conn.execute(
                "SELECT entry_id FROM entries WHERE elements IS NULL OR trim(elements) = '' LIMIT 5"
            )
        ]
        rep.warn(
            "E2.elements",
            f"{empty:,} entries carry no element list ({empty * 100.0 / total:.2f}%) "
            f"-- e.g. {', '.join(samples)}",
        )
    else:
        rep.ok("E2.elements", "every entry lists at least one element")

    no_sg = int(
        conn.execute(
            "SELECT count(*) FROM entries WHERE space_group IS NULL OR trim(space_group) = ''"
        ).fetchone()[0]
    )
    if no_sg:
        rep.info("E3.space_group", f"{no_sg:,} entries without a space group")

    # ``top_k`` is the *prescreen* width (32 reflections per phase), not the
    # per-phase storage cap: the peak table keeps up to DEFAULT_MAX_PEAKS lines
    # so a GUI can plot the full stick pattern.  Conflating the two would flag
    # every library as broken.
    if max_peaks > DEFAULT_MAX_PEAKS:
        rep.fail(
            "E4.max_peaks",
            f"max n_peaks = {max_peaks} exceeds DEFAULT_MAX_PEAKS = {DEFAULT_MAX_PEAKS}",
        )
    else:
        rep.ok("E4.max_peaks", f"n_peaks <= {DEFAULT_MAX_PEAKS} (max {max_peaks})")


def _decode_lone_reflections(
    conn: sqlite3.Connection,
    rows: list[tuple[int, str, int]],
    chunk_size: int,
) -> list[tuple[str, float | None]]:
    """Return ``(entry_id, d)`` for every ``n_peaks == 1`` entry in *rows*.

    Only the handful of single-reflection phases are decoded, but every
    ``peak_chunks`` row they live in is fetched at most once: on the 1.7 GiB
    PCOD library this costs six small block reads instead of a full payload
    scan, which is what makes ``E5`` affordably *exhaustive*.

    ``d`` is ``None`` when the stored slice cannot be located (a symptom the
    ``C``/``P`` checks report separately; ``E5`` must not crash on it).
    """
    cache: dict[int, tuple[int, bytes]] = {}
    decoded: list[tuple[str, float | None]] = []
    for rowid, entry_id, peak_off in rows:
        chunk_id = int(rowid) // chunk_size
        cached = cache.get(chunk_id)
        if cached is None:
            row = conn.execute(
                "SELECT n, all_d FROM peak_chunks WHERE chunk_id = ?", (chunk_id,)
            ).fetchone()
            if row is None:
                cached = (-1, b"")
            else:
                base = conn.execute(
                    "SELECT peak_off FROM entries WHERE rowid = ?",
                    (chunk_id * chunk_size,),
                ).fetchone()
                cached = (
                    int(base[0]) if base is not None else -1,
                    bytes(row[1] or b""),
                )
            cache[chunk_id] = cached
        base, payload = cached
        lo = int(peak_off) - base
        value: float | None = None
        if base >= 0 and lo >= 0 and (lo + 1) * 4 <= len(payload):
            value = float(np.frombuffer(payload[lo * 4 : (lo + 1) * 4], dtype="<f4")[0])
        decoded.append((str(entry_id), value))
    return decoded


def _format_d(value: float | None) -> str:
    return "?" if value is None else f"{value:.4f}"


def _scan_boundary_strongest(
    conn: sqlite3.Connection,
    rep: Report,
    chunk_size: int,
    top_k: int,
) -> tuple[
    list[tuple[str, float, int]],
    list[tuple[str, float, int, float]],
]:
    """Return the two boundary-phase classes keyed on the strongest reflection.

    The prescreen block is exactly what is needed: ``top_d`` rows are ordered by
    intensity, so column 0 is the phase's *strongest* reflection, and column 0
    of ``top_i`` is ~1000 by the ``norm-1000`` convention.  A phase is bucketed
    by the strongest reflection's d and by the sum of every *other* stored
    intensity:

    * ``collapsed`` -- strongest d in the boundary window **and** the rest of
      the pattern summed below :data:`BOUNDARY_COLLAPSE_MAX` (the actual
      pre-2026-10 normalization bug: the exploded 180 deg line became the whole
      pattern).  ``(entry_id, d, n_peaks)``.
    * ``semi_collapsed`` -- strongest d in the boundary window but the pattern
      is still populated.  ``(entry_id, d, n_peaks, rest_sum)``; these are
      **WARN**, not FAIL: a real reflection near 2theta = 180 deg legitimately
      becomes the strongest line because the Lorentz-polarisation factor
      diverges there, and thousands of healthy phases fall in this bucket.
      Quoting the smallest ``rest_sum`` first surfaces the borderline cases
      that look most like the collapse artifact.

    The scan is exhaustive (every chunk, every phase) and independent of
    ``--sample-chunks``; it costs one ``top_d``/``top_i`` decode per chunk, i.e.
    ~136 MB on the 1.06 M-phase PCOD library.
    """
    collapsed: list[tuple[str, float, int]] = []
    semi_collapsed: list[tuple[str, float, int, float]] = []
    for chunk_id, n, top_d_blob, top_i_blob in conn.execute(
        "SELECT chunk_id, n, top_d, top_i FROM peak_chunks ORDER BY chunk_id"
    ):
        n = int(n)
        top_d = _decode_matrix(top_d_blob, n, top_k)
        top_i = _decode_matrix(top_i_blob, n, top_k)
        strongest = top_d[:, 0]
        rest_sum = top_i[:, 1:].sum(axis=1)
        window = (strongest >= DEFAULT_DMIN - DMIN_LOW_SLACK) & (
            strongest <= DEFAULT_DMIN + DMIN_HIGH_SLACK
        )
        if not window.any():
            continue
        collapsed_mask = window & (rest_sum < BOUNDARY_COLLAPSE_MAX)
        semi_mask = window & ~collapsed_mask
        entries = conn.execute(
            "SELECT entry_id, n_peaks FROM entries WHERE rowid >= ? AND rowid < ? "
            "ORDER BY rowid",
            (int(chunk_id) * chunk_size, int(chunk_id) * chunk_size + n),
        ).fetchall()
        if len(entries) != n:
            rep.warn(
                "E5.boundary_peak",
                f"chunk {int(chunk_id)}: {len(entries)} entries at that offset, "
                f"chunk says {n} (see C2/C5)",
            )
            continue
        for index in np.flatnonzero(collapsed_mask):
            entry_id, n_peaks = entries[int(index)]
            collapsed.append((str(entry_id), float(strongest[int(index)]), int(n_peaks)))
        for index in np.flatnonzero(semi_mask):
            entry_id, n_peaks = entries[int(index)]
            semi_collapsed.append(
                (
                    str(entry_id),
                    float(strongest[int(index)]),
                    int(n_peaks),
                    float(rest_sum[int(index)]),
                )
            )
    return collapsed, semi_collapsed


def _check_peak_inventory(
    conn: sqlite3.Connection,
    rep: Report,
    meta: dict[str, str],
) -> None:
    """``E5.boundary_peak`` / ``E6.peak_count``: the two known bad-phase classes.

    ``E5`` is anchored on the **strongest reflection** of every phase (encoded
    in ``peak_chunks.top_d`` column 0), not just on single-reflection phases, so
    it also catches a pattern that collapsed to two lines (the real structure
    keeps one tiny line that survived the bogus normalization).

    Three signatures, exhaustive over every ``entries`` row:

    * **180 deg boundary artifact** -- strongest reflection at the
      ``DEFAULT_DMIN`` cutoff (``lambda/2``) *and* the rest of the pattern
      summed far below the norm scale (see :data:`BOUNDARY_COLLAPSE_MAX`) ->
      **FAIL**;
    * **placeholder-cell artifact** -- a lone reflection at 1.0 A, which is the
      ``a = b = c = 1.0 A`` cell gemmi invents for a CIF without ``_cell_*``
      tags -> **FAIL**;
    * any other lone reflection is physically possible (a large cell with one
      in-window line) -> **WARN**.

    The first class used to be detected only through ``n_peaks == 1``; that
    missed e.g. PCOD ``96-254-30089`` (n_peaks = 2: the bogus 180 deg line at
    1000 plus one real line at 0.013).  Anchoring on the strongest reflection
    generalises the check without changing its cost class.
    """
    try:
        chunk_size = int(meta.get("chunk_size", "0"))
        top_k = int(meta.get("top_k", "0"))
    except ValueError:
        chunk_size = top_k = 0
    if chunk_size <= 0 or top_k <= 0:
        rep.warn(
            "E5.boundary_peak",
            f"cannot enumerate bad phases: meta.chunk_size={meta.get('chunk_size')!r} "
            f"meta.top_k={meta.get('top_k')!r}",
        )
        return

    # -- E6: a phase that stores no reflection at all -----------------------
    empty = [
        str(row[0])
        for row in conn.execute(
            "SELECT entry_id FROM entries WHERE n_peaks IS NULL OR n_peaks <= 0 ORDER BY rowid"
        )
    ]
    if empty:
        rep.fail(
            "E6.peak_count",
            f"{len(empty):,} phases store no reflection at all "
            f"(n_peaks <= 0; e.g. {', '.join(empty[:5])})",
        )
    else:
        rep.ok("E6.peak_count", "every phase stores at least one reflection")

    # -- E5a: strongest reflection at the 2theta=180 deg cutoff -------------
    boundary, semi_collapsed = _scan_boundary_strongest(conn, rep, chunk_size, top_k)
    boundary_ids = {entry_id for entry_id, _d, _n in boundary}

    # -- E5b/E5c: lone reflections (1.0 A placeholder vs plausible) ---------
    singles = [
        (int(rowid), str(entry_id), int(peak_off))
        for rowid, entry_id, peak_off in conn.execute(
            "SELECT rowid, entry_id, peak_off FROM entries WHERE n_peaks = 1 ORDER BY rowid"
        )
    ]
    placeholder: list[tuple[str, float]] = []
    normal: list[tuple[str, float]] = []
    undecodable: list[str] = []
    for entry_id, value in _decode_lone_reflections(conn, singles, chunk_size):
        if value is None:
            undecodable.append(entry_id)
        elif entry_id in boundary_ids:
            continue  # already reported by the strongest-reflection scan
        elif abs(value - PLACEHOLDER_D) <= PLACEHOLDER_TOL:
            placeholder.append((entry_id, value))
        else:
            normal.append((entry_id, value))

    if boundary or placeholder:
        parts: list[str] = []
        if boundary:
            parts.append(
                f"{len(boundary)} with the strongest reflection at the "
                f"{DEFAULT_DMIN} A cutoff (2theta=180 deg) and the rest of the "
                f"pattern collapsed"
            )
        if placeholder:
            parts.append(
                f"{len(placeholder)} holding a lone {PLACEHOLDER_D} A "
                f"placeholder-cell line"
            )
        samples = [
            f"{entry_id}(d={_format_d(d)})" for entry_id, d, _n in boundary[:MAX_QUOTED_IDS]
        ]
        samples += [
            f"{entry_id}(d={_format_d(d)})" for entry_id, d in placeholder[:MAX_QUOTED_IDS]
        ]
        quoted = samples[:MAX_QUOTED_IDS]
        total = len(boundary) + len(placeholder)
        more = "" if total <= MAX_QUOTED_IDS else f" ... {total - MAX_QUOTED_IDS} more"
        rep.fail(
            "E5.boundary_peak",
            f"{total} bogus phases: " + "; ".join(parts) + f"; e.g. {', '.join(quoted)}{more}",
        )
    if normal:
        samples = ", ".join(
            f"{entry_id}(d={_format_d(value)})"
            for entry_id, value in normal[:MAX_QUOTED_IDS]
        )
        rep.warn(
            "E5.boundary_peak",
            f"{len(normal)} phases hold exactly one reflection at a plausible "
            f"d-spacing (not a known engine artifact): {samples}",
        )
    if semi_collapsed:
        # Same root cause (a reflection computed at lambda/2), but the pattern
        # is still populated, so a real reflection near 180 deg is the strongest
        # line.  Not a build defect: WARN, most-collapse-like cases first.
        semi_collapsed.sort(key=lambda item: item[3])
        samples = ", ".join(
            f"{entry_id}(d={_format_d(d)}, n={n_peaks}, rest_sum={rest_sum:.3f})"
            for entry_id, d, n_peaks, rest_sum in semi_collapsed[:MAX_BOUNDARY_WARN_QUOTED]
        )
        rep.warn(
            "E5.boundary_peak",
            f"{len(semi_collapsed)} phases have their strongest reflection in the "
            f"{DEFAULT_DMIN} A cutoff window but a still-populated pattern, so they are "
            f"NOT counted as bad phases (a real reflection near 2theta=180 deg is "
            f"legitimately the strongest line there); closest to the collapse boundary: "
            f"{samples}",
        )
    if undecodable:
        rep.warn(
            "E5.boundary_peak",
            f"{len(undecodable)} lone reflections could not be decoded "
            f"(see C/P checks): {', '.join(undecodable[:5])}",
        )
    if not boundary and not placeholder and not normal:
        rep.ok("E5.boundary_peak", "no bogus boundary/placeholder phase")


def _check_payload(
    conn: sqlite3.Connection,
    rep: Report,
    meta: dict[str, str],
    sample: int,
    seed: int,
) -> None:
    try:
        top_k = int(meta.get("top_k", "0"))
        chunk_size = int(meta.get("chunk_size", "0"))
    except ValueError:
        rep.fail("P0.config", "meta.top_k / meta.chunk_size unreadable")
        return
    if top_k <= 0 or chunk_size <= 0:
        rep.fail("P0.config", f"top_k={top_k} chunk_size={chunk_size}")
        return

    ids = [int(r[0]) for r in conn.execute("SELECT chunk_id FROM peak_chunks ORDER BY chunk_id")]
    if not ids:
        rep.fail("P1.sample", "no chunks to sample")
        return
    rng = random.Random(seed)
    # ``sample < 0`` is the ``--all-chunks`` sentinel: decode every block.
    if sample < 0:
        picked = sorted(ids)
        rep.info("P1.all_chunks", f"decoded all {len(picked):,} chunks")
    else:
        picked = sorted(rng.sample(ids, min(sample, len(ids))))

    problems: list[str] = []
    checked_entries = 0
    for cid in picked:
        row = conn.execute(
            "SELECT n, top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?",
            (cid,),
        ).fetchone()
        n = int(row[0])
        top_d = _decode_matrix(row[1], n, top_k)
        top_i = _decode_matrix(row[2], n, top_k)
        all_d = _decode_flat(row[3])
        all_i = _decode_flat(row[4])

        if len(all_d) != len(all_i):
            problems.append(f"chunk {cid}: all_d/all_i length {len(all_d)} vs {len(all_i)}")
            continue

        rows = conn.execute(
            "SELECT n_peaks, peak_off FROM entries WHERE rowid >= ? AND rowid < ? ORDER BY rowid",
            (cid * chunk_size, cid * chunk_size + n),
        ).fetchall()
        if len(rows) != n:
            problems.append(f"chunk {cid}: {len(rows)} entries at that offset, chunk says {n}")
            continue

        # ``peak_off`` is a *global* reflection index, but each chunk's ``all_d``
        # holds only its own reflections -- so the slice origin is the peak_off
        # of the chunk's first phase, not its entry index.  Subtracting the entry
        # index instead shifts every slice by the number of reflections already
        # stored and makes a perfectly good library look corrupt.
        chunk_peak_base = int(rows[0][1])

        if not np.all(np.isfinite(top_d)) or not np.all(np.isfinite(top_i)):
            problems.append(f"chunk {cid}: non-finite values in prescreen block")
        if not np.all(np.isfinite(all_d)) or not np.all(np.isfinite(all_i)):
            problems.append(f"chunk {cid}: non-finite values in full peak block")
        if top_d.size and float(top_d.min()) < 0:
            problems.append(f"chunk {cid}: negative d-spacing in prescreen block")
        if all_d.size and float(all_d.min()) <= 0:
            problems.append(f"chunk {cid}: non-positive d-spacing in full block")

        base = cid * chunk_size
        for i, (n_peaks, peak_off) in enumerate(rows):
            n_peaks = int(n_peaks)
            lo = int(peak_off) - chunk_peak_base
            checked_entries += 1
            if n_peaks <= 0:
                continue
            if lo < 0 or lo + n_peaks > len(all_d):
                problems.append(
                    f"chunk {cid} entry {base + i}: slice [{lo}, {lo + n_peaks}) "
                    f"outside payload of {len(all_d)}"
                )
                continue
            d_slice = all_d[lo : lo + n_peaks]
            i_slice = all_i[lo : lo + n_peaks]
            if float(i_slice.max()) > NORM_SCALE + NORM_TOL:
                problems.append(
                    f"chunk {cid} entry {base + i}: intensity {float(i_slice.max()):.1f} "
                    f"above norm scale {NORM_SCALE:.0f}"
                )
            if abs(float(i_slice.max()) - NORM_SCALE) > NORM_TOL:
                problems.append(
                    f"chunk {cid} entry {base + i}: strongest reflection "
                    f"{float(i_slice.max()):.3f} != {NORM_SCALE:.0f}"
                )
            # The prescreen row must be a subset of the full list: same values,
            # same count -- a mismatch means the two blocks came from different
            # computations.
            take = min(top_k, n_peaks)
            stored = top_d[i, :take]
            if take and not np.all(np.isin(np.round(stored, 5), np.round(d_slice, 5))):
                problems.append(
                    f"chunk {cid} entry {base + i}: prescreen d values absent from full list"
                )
            if take and not np.all(np.isin(np.round(top_i[i, :take], 3), np.round(i_slice, 3))):
                problems.append(
                    f"chunk {cid} entry {base + i}: prescreen intensities absent from full list"
                )
            pad_d = top_d[i, take:]
            pad_i = top_i[i, take:]
            if pad_d.size and (np.any(pad_d != 0.0) or np.any(pad_i != 0.0)):
                problems.append(f"chunk {cid} entry {base + i}: prescreen padding is not zero")
            if problems and len(problems) >= 20:
                break
        if problems and len(problems) >= 20:
            break

    if problems:
        for msg in problems[:10]:
            rep.add("P1.payload", "FAIL", msg)
        if len(problems) > 10:
            rep.add("P1.payload", "FAIL", f"... {len(problems) - 10} more")
    else:
        rep.ok(
            "P1.payload",
            f"{len(picked)} chunks / {checked_entries:,} entries decode cleanly, norm-1000 held",
        )


def _check_element_index(
    conn: sqlite3.Connection,
    rep: Report,
    entries: int,
    sample: int,
    seed: int,
) -> None:
    try:
        total = int(conn.execute("SELECT count(*) FROM element_inv").fetchone()[0])
    except sqlite3.Error:
        rep.fail("X1.present", "element_inv table missing or unreadable")
        return
    if total == 0:
        rep.fail("X1.present", "element_inv is empty -- element filtering cannot work")
        return
    rep.ok("X1.present", f"{total} element entries in the inverted index")

    bad_range = 0
    bad_order = 0
    for element, blob in conn.execute("SELECT element, eids FROM element_inv"):
        arr = _decode_eids(blob)
        if arr.size == 0:
            continue
        if int(arr.max()) >= entries:
            bad_range += 1
            rep.add(
                "X2.range",
                "FAIL",
                f"element {element}: eid {int(arr.max()):,} >= entry count {entries:,}",
            )
        if arr.size > 1 and np.any(np.diff(arr.astype(np.int64)) <= 0):
            bad_order += 1
    if bad_range == 0:
        rep.ok("X2.range", "every eid is a valid entry index")
    if bad_order == 0:
        rep.ok("X3.order", "every eid list is strictly ascending and unique")
    else:
        rep.warn("X3.order", f"{bad_order} element lists are not strictly ascending")

    rng = random.Random(seed + 1)
    offsets = sorted(rng.sample(range(entries), min(sample, entries)))
    if not offsets:
        return
    # ``entries`` is keyed by an explicit ``rowid`` equal to the insertion index
    # (0-based) -- the same index ``element_inv`` stores.  Reading it back as
    # ``rowid - 1`` shifts every lookup by one phase and reports the whole index
    # as broken.
    placeholders = ",".join("?" * len(offsets))
    picked = conn.execute(
        f"SELECT rowid, elements FROM entries WHERE rowid IN ({placeholders}) ORDER BY rowid",
        offsets,
    ).fetchall()
    if not picked:
        rep.warn("X4.coverage", "no sampled entries could be read back")
        return

    symbols = {s for _rid, txt in picked for s in _split_elements(txt)}
    wanted: dict[str, set[int]] = {}
    for sym in symbols:
        row = conn.execute("SELECT eids FROM element_inv WHERE element = ?", (sym,)).fetchone()
        wanted[sym] = set(_decode_eids(row[0]).tolist()) if row else set()

    misses = 0
    for rid, txt in picked:
        idx = int(rid)
        for sym in _split_elements(txt):
            if idx not in wanted.get(sym, set()):
                misses += 1
                if misses <= 5:
                    rep.add(
                        "X4.coverage",
                        "FAIL",
                        f"entry {idx:,} lists {sym} but is absent from element_inv[{sym}]",
                    )
    if misses == 0:
        rep.ok(
            "X4.coverage",
            f"{len(picked)} sampled phases found under all {len(symbols)} of their elements",
        )
    else:
        rep.add("X4.coverage", "FAIL", f"{misses} element/phase lookups failed")


def _decode_eids(blob: object) -> np.ndarray:
    data = b"" if blob is None else bytes(blob)
    count = len(data) // 4
    return np.frombuffer(data[: count * 4], dtype=np.uint32)


def _split_elements(raw: object) -> list[str]:
    text = str(raw or "")
    for sep in (";", ","):
        text = text.replace(sep, " ")
    return [tok.strip().title() for tok in text.split() if tok.strip()]


def validate(db: Path, sample_chunks: int, sample_entries: int, seed: int = 0) -> Report:
    """Run every check against *db* and return the populated report.

    ``sample_chunks`` is the number of ``peak_chunks`` blocks decoded by the
    payload check; a **negative** value means "all chunks" (the
    ``--all-chunks`` mode).  The ``E5``/``E6`` reflection-count checks are
    always exhaustive and do not depend on this parameter.
    """
    rep = Report()
    if not db.exists():
        rep.fail("F0.exists", f"no such file: {db}")
        return rep
    size_mb = db.stat().st_size / (1024.0 * 1024.0)
    rep.info("F1.file", f"{db.name}  {size_mb:,.1f} MB")

    uri = f"file:{db.as_posix()}?mode=ro"
    try:
        conn = sqlite3.connect(uri, uri=True, timeout=30.0)
    except sqlite3.Error as exc:
        rep.fail("F0.exists", f"cannot open: {exc}")
        return rep

    t0 = time.time()
    try:
        if not _check_schema(conn, rep):
            return rep
        meta = _check_meta(conn, rep, db)
        entries, _total_peaks, max_peaks = _check_counts(conn, rep, meta)
        _check_entries(conn, rep, meta, max_peaks)
        _check_peak_inventory(conn, rep, meta)
        _check_payload(conn, rep, meta, sample_chunks, seed)
        _check_element_index(conn, rep, entries, sample_entries, seed)
    except sqlite3.Error as exc:
        rep.fail("F2.sqlite", f"{type(exc).__name__}: {exc}")
    finally:
        conn.close()

    rep.info("F3.elapsed", f"{time.time() - t0:.1f}s")
    return rep


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("db", type=Path, help=".xinDB file to validate (opened read-only)")
    ap.add_argument(
        "--sample-chunks",
        type=int,
        default=8,
        help="chunks decoded and checked (default 8; see --all-chunks)",
    )
    ap.add_argument(
        "--all-chunks",
        action="store_true",
        help="decode every chunk instead of a sample (wins over --sample-chunks); "
        "the report states how many chunks were actually decoded",
    )
    ap.add_argument("--sample-entries", type=int, default=300, help="entries checked against element_inv")
    ap.add_argument("--seed", type=int, default=0, help="sampling seed")
    ap.add_argument("--json", type=Path, default=None, help="write the report as JSON here")
    args = ap.parse_args(argv)

    sample_chunks = -1 if args.all_chunks else args.sample_chunks
    rep = validate(args.db, sample_chunks, args.sample_entries, args.seed)
    print(f"validate_xindb: {args.db}")
    print(rep.render())
    verdict = "FAILED" if rep.fails else ("PASSED WITH WARNINGS" if rep.warns else "PASSED")
    print(f"\n  {verdict}: {rep.fails} fail, {rep.warns} warn, {len(rep.rows)} checks")

    if args.json:
        args.json.write_text(
            json.dumps(
                {
                    "db": str(args.db),
                    "verdict": verdict,
                    "fails": rep.fails,
                    "warns": rep.warns,
                    "checks": [
                        {"code": c, "status": s, "detail": d} for c, s, d in rep.rows
                    ],
                },
                indent=2,
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        print(f"  json -> {args.json}")
    return 1 if rep.fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
