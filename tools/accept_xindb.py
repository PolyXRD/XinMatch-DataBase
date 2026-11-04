"""One-command acceptance report for a built ``.xinDB`` reference library.

:mod:`tools.validate_xindb` answers "does this file satisfy the v3 format that
any reader assumes?" -- schema, metadata, chunk accounting, payload shapes,
element inverted index.  It is genuinely library-agnostic: it does not know (or
care) *which* of the six shipped products it is looking at.

This tool answers the product-level question:

    **Is this particular library correctly built from the archive it claims to
    come from, with every id owned by the right source and no known bad phase
    hiding in it?**

It covers the six products ``COD / PCOD / TCOD`` x ``Inorg / Full`` and derives
everything it needs from ``meta`` and the file name:

======================  ============  ===================================
product                 id offset      default source root
======================  ============  ===================================
``COD-Inorg``           0              ``<root>/inorg_cif/COD``
``COD-Full``            0              ``<root>/COD``
``PCOD-Inorg``          20,000,000     ``<root>/inorg_cif/PCOD``
``PCOD-Full``           20,000,000     ``<root>/PCOD``
``TCOD-Inorg``          100,000,000    ``<root>/inorg_cif/TCOD``
``TCOD-Full``           100,000,000    ``<root>/TCOD``
======================  ============  ===================================

``<root>`` is ``$XRD_DATA_ROOT`` (default ``D:/Project/XRD/sync-COD``); a
``meta.cif_source`` that still exists wins over it, so a library built from a
different archive is checked against that archive.  ``--manifest`` points at
the shared ``inorg_manifest.sqlite`` (the scan of the three full archives) and
is used as a fallback ownership witness, plus as a cross-check of the
inorganic membership flag.

Checks
------
``F`` -- file / schema
    The file exists, opens read-only, and carries the four v3 tables.
``M`` -- metadata
    ``schema_version`` is 3, ``top_k`` / ``chunk_size`` / ``entry_count`` are
    readable and consistent, and the source / dataset can be derived.
``C`` -- counts
    ``meta.entry_count`` == ``entries`` rows, ``SUM(peak_chunks.n)`` covers
    them, and ``SUM(entries.n_peaks)`` equals the stored reflection payload.
``R`` -- id ownership (full scan)
    Every ``entry_id`` is the canonical ``96-XXX-XXXX`` form of
    ``stored id = archive id + offset`` and lies inside this source's id band
    (COD ``[0, 20M)``, PCOD ``[20M, 100M)``, TCOD ``[100M, 200M)``).  A wrong
    offset or a cross-source id collision is caught here, not by the reader.
``A`` -- archive ownership (sample)
    ``--sample-ids`` (40) sampled ids are mapped back to ``<archive id>.cif``
    in the source root (or looked up in the manifest) so "every phase really
    came from this archive" is evidenced, not assumed.
``E`` -- element index
    The ``element_inv`` symbol set equals the set used by ``entries.elements``
    (full scan), every ``eid`` is a valid entry index, and sampled memberships
    hold in both directions.
``P`` -- peak payload
    Every chunk (default) or ``--sample-chunks`` blocks decode to the
    documented shapes; reflections are finite, positive and on the declared
    ``norm-1000`` scale; the prescreen block is a subset of the full list and
    its padding is zero.
``B`` -- bad phases
    Exhaustive scan anchored on the **strongest reflection** of every phase
    (``peak_chunks.top_d`` column 0): a phase whose strongest line sits at the
    d-spacing cutoff (2theta = 180 deg) *and* whose remaining intensity has
    collapsed onto it (the old Lorentz-polarisation/normalisation defect) fails
    the library -- this catches both the fully collapsed single-line pattern and
    the partially collapsed two-line case (PCOD ``96-254-30089``, n_peaks = 2).
    A lone 1.0 A reflection (placeholder-cell artifact) fails too, as does
    ``n_peaks <= 0``.  Phases whose strongest reflection is in the cutoff
    window but whose pattern is still populated are **WARN** only (up to 10
    examples + total, reported as ``B1.boundary_warn`` / ``stats.boundary_warn``
    -- e.g. PCOD ``96-244-04910``; they share the root cause but are not the
    build defect).

The tool imports **only the standard library and NumPy** -- it must run on a
machine that has neither XinMatch nor its peak engine installed.

Exit code is 1 when any check fails.

Usage::

    python -m tools.accept_xindb data/COD-Inorg-20261020.xinDB --json out.json
    python -m tools.accept_xindb data/PCOD-Inorg-20261020.xinDB --sample-chunks 8
    python -m tools.accept_xindb data/TCOD-Inorg-20261020.xinDB \\
        --manifest D:/Project/XRD/sync-COD/inorg_manifest.sqlite
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import sqlite3
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

#: Shared source root, same convention as the build tools (``XRD_DATA_ROOT``).
DEFAULT_DATA_ROOT = Path(os.environ.get("XRD_DATA_ROOT") or r"D:\Project\XRD\sync-COD")

#: Materials Project CIF root (``MP_CIF_ROOT`` overrides it).
DEFAULT_MP_ROOT = Path(os.environ.get("MP_CIF_ROOT") or r"D:\Project\MP-summary\cif")

#: Source id offsets.  ``MP`` uses native ``mp-<n>`` ids (no offset).
SOURCES: dict[str, tuple[int, str]] = {
    "COD": (0, "COD"),
    "PCOD": (20_000_000, "PCOD"),
    "TCOD": (100_000_000, "TCOD"),
    "MP": (0, ""),
}

#: The stored-id band each source owns.  A COD phase that lands in the PCOD
#: band means the build applied the wrong offset; that is a FAIL, not a WARN.
BANDS: dict[str, tuple[int, int]] = {
    "COD": (0, 20_000_000),
    "PCOD": (20_000_000, 100_000_000),
    "TCOD": (100_000_000, 200_000_000),
    "MP": (0, 20_000_000),
}

#: Intensity scale the ``norm-1000`` convention promises (fallback when the
#: stored ``peak_convention`` does not spell it out).
NORM_SCALE = 1000.0

#: Slack for the normalisation check (float32 division can be a few ULP off).
NORM_TOL = 0.05

#: Fallback d-spacing cutoff (A) when ``meta.peak_convention`` is unreadable.
DEFAULT_DMIN = 0.7703

#: Boundary window: the bogus 180 deg line sits an ulp above lambda/2; the
#: float32 rounding of 0.7703 is up to 3e-8 low, so [d_min-1e-6, d_min+1e-3].
DMIN_LOW_SLACK = 1e-6
DMIN_HIGH_SLACK = 1e-3

#: A collapsed 180 deg artifact puts essentially the whole normalised intensity
#: into the bogus boundary line: everything else sums to less than this on the
#: ``norm-1000`` scale (0.5 %).  This second condition is what keeps the
#: strongest-reflection test specific: the Lorentz-polarisation factor diverges
#: near 2theta = 180 deg, so a *real* reflection there also becomes the
#: strongest line (thousands of healthy phases would otherwise be flagged).
BOUNDARY_COLLAPSE_MAX = 5.0

#: The placeholder unit cell gemmi invents for a CIF without ``_cell_*`` tags.
PLACEHOLDER_D = 1.0
PLACEHOLDER_TOL = 1e-3

#: Cap on how many bad-phase ids are quoted in a check detail line.
MAX_QUOTED_IDS = 12

#: Cap on how many *semi-collapsed* boundary phases (strongest reflection at the
#: cutoff, pattern still populated) a WARN line quotes, plus the total.
MAX_BOUNDARY_WARN_QUOTED = 10

#: Cap on how many bad-phase ids are carried in the JSON stats block.
MAX_STATS_IDS = 500

#: Maximum number of accumulated payload problems before the scan stops.
MAX_PAYLOAD_PROBLEMS = 50

#: Tables a v3 library must carry.  ``element_inv`` is required here (unlike
#: the generic validator) because element filtering is an acceptance criterion.
REQUIRED_TABLES: tuple[str, ...] = ("meta", "entries", "peak_chunks", "element_inv")


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
@dataclass
class Report:
    """Collected check results plus the acceptance statistics."""

    db: Path
    kind: str = "?"
    dataset: str = "?"
    rows: list[tuple[str, str, str]] = field(default_factory=list)
    fails: int = 0
    warns: int = 0
    stats: dict[str, object] = field(default_factory=dict)

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

    @property
    def verdict(self) -> str:
        """``PASS`` / ``PASS WITH WARNINGS`` / ``FAIL``."""
        if self.fails:
            return "FAIL"
        return "PASS WITH WARNINGS" if self.warns else "PASS"

    def render(self) -> str:
        width = max((len(r[0]) for r in self.rows), default=4)
        return "\n".join(
            f"  [{status}] {code:<{width}}  {detail}" for code, status, detail in self.rows
        )

    def as_dict(self) -> dict[str, object]:
        """JSON-friendly view (``verdict/fails/warns/checks/stats``)."""
        return {
            "db": str(self.db),
            "kind": self.kind,
            "dataset": self.dataset,
            "verdict": self.verdict,
            "fails": self.fails,
            "warns": self.warns,
            "checks": [
                {"code": code, "status": status, "detail": detail}
                for code, status, detail in self.rows
            ],
            "stats": dict(self.stats),
        }


# ---------------------------------------------------------------------------
# Small decoding helpers (kept local on purpose: this tool must stay importable
# without the application package)
# ---------------------------------------------------------------------------
def _read_meta(conn: sqlite3.Connection) -> dict[str, str]:
    try:
        return {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
    except sqlite3.Error:
        return {}


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


def _decode_eids(blob: object) -> np.ndarray:
    data = b"" if blob is None else bytes(blob)
    count = len(data) // 4
    return np.frombuffer(data[: count * 4], dtype="<u4")


def _split_elements(raw: object) -> list[str]:
    text = str(raw or "").replace(";", ",")
    return [tok.strip() for tok in text.split(",") if tok.strip()]


def _entry_number(entry_id: object) -> int | None:
    """Recover the *stored* numeric id from a ``96-XXX-XXXX`` display id."""
    parts = str(entry_id or "").split("-")
    if len(parts) < 2:
        return None
    digits = "".join(parts[1:])
    if not digits.isascii() or not digits.isdigit():
        return None
    return int(digits)


def _entry_number_for(entry_id: object, kind: str) -> int | None:
    """Numeric id of an MP (``mp-<n>``) or COD-family (``96-...``) entry."""
    if kind == "MP":
        text = str(entry_id or "")
        if text.startswith("mp-") and text[3:].isascii() and text[3:].isdigit():
            return int(text[3:])
        return None
    return _entry_number(entry_id)


def _display_id(stored: int) -> str:
    """Canonical display id of a *stored* (offset-applied) numeric id.

    Must mirror ``cif_scan.entry_id_for_cod`` byte for byte: seven-digit
    zero padding, split after the first three digits.
    """
    text = str(int(stored)).zfill(7)
    return f"96-{text[:3]}-{text[3:]}"


def _source_cif_candidates(archive_id: int) -> list[str]:
    """Archive-relative CIF paths to try for *archive_id* (source id, no offset).

    The three archives fan out differently but all slice the zero-padded id:

    * COD -- three levels, two digits each: ``1/00/00/1000005.cif``;
    * TCOD (8-digit ids) -- three levels, two digits each:
      ``10/00/00/10000001.cif``;
    * PCOD -- two levels where the second level is the **first three digits**
      (verified on the real archive: ``1/100/1000000.cif``, ``3/310/3100000.cif``).

    Returning all shapes keeps the ownership check *data-driven*: it stats the
    paths instead of hard-coding which archive this library came from.
    """
    text = str(int(archive_id))
    if len(text) < 7:
        # No archive contains ids shorter than the 7-digit COD layout.
        return []
    candidates = [f"{text[:-6]}/{text[-6:-4]}/{text[-4:-2]}/{text}.cif"]
    if len(text) >= 4:
        candidates.append(f"{text[0]}/{text[:3]}/{text}.cif")
        candidates.append(f"{text[0]}/{text[1:4]}/{text}.cif")
    return candidates


def _cutoff_from_meta(meta: dict[str, str]) -> float:
    """d-spacing cutoff recorded in ``meta.peak_convention`` (fallback 0.7703)."""
    match = re.search(r"-d([0-9]+\.[0-9]+)-", meta.get("peak_convention", ""))
    return float(match.group(1)) if match else DEFAULT_DMIN


def _norm_scale_from_meta(meta: dict[str, str]) -> float:
    """Intensity scale recorded in ``meta.peak_convention`` (fallback 1000)."""
    match = re.search(r"norm([0-9]+)", meta.get("peak_convention", ""))
    return float(match.group(1)) if match else NORM_SCALE


# ---------------------------------------------------------------------------
# Derivation of kind / dataset / source root
# ---------------------------------------------------------------------------
def _derive_kind(db: Path, meta: dict[str, str], explicit: str) -> str | None:
    """Source archive (``COD``/``PCOD``/``TCOD``/``MP``) of the library."""
    if explicit and explicit != "auto":
        return explicit.upper()
    for text in (meta.get("library_kind", ""), db.stem):
        match = re.match(r"(MP|PCOD|TCOD|COD)(?:-|$)", str(text), re.IGNORECASE)
        if match:
            return match.group(1).upper()
    return None


def _derive_dataset(meta: dict[str, str]) -> str:
    """``inorg`` / ``full`` (``unknown`` when neither meta nor path says)."""
    library_kind = str(meta.get("library_kind", ""))
    if library_kind.endswith("-Inorg"):
        return "inorg"
    if library_kind.endswith("-Full") or library_kind.upper() == "MP":
        return "full"
    source = str(meta.get("cif_source", "")).replace("\\", "/").lower()
    if "inorg_cif" in source:
        return "inorg"
    if source:
        return "full"
    return "unknown"


def _source_root(
    meta: dict[str, str], kind: str, dataset: str, data_root: Path
) -> tuple[Path | None, str]:
    """Best existing CIF archive root for this library, plus how it was found."""
    candidates: list[tuple[Path, str]] = []
    stored = str(meta.get("cif_source", "")).strip()
    if stored:
        candidates.append((Path(stored), "meta.cif_source"))
    if dataset == "inorg":
        candidates.append((data_root / "inorg_cif" / kind, "XRD_DATA_ROOT/inorg_cif"))
        candidates.append((data_root / kind, "XRD_DATA_ROOT/full"))
    else:
        candidates.append((data_root / kind, "XRD_DATA_ROOT"))
        candidates.append((data_root / "inorg_cif" / kind, "XRD_DATA_ROOT/inorg_cif"))
    for path, origin in candidates:
        if path.is_dir():
            return path, origin
    return None, "not found"


# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------
def _check_schema(conn: sqlite3.Connection, rep: Report) -> bool:
    try:
        names = {
            str(r[0])
            for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")
        }
    except sqlite3.Error as exc:
        rep.fail("F2.sqlite", f"{type(exc).__name__}: {exc}")
        return False
    missing = [table for table in REQUIRED_TABLES if table not in names]
    if missing:
        rep.fail("S1.tables", f"missing required tables: {', '.join(missing)}")
        return False
    rep.ok("S1.tables", f"{len(REQUIRED_TABLES)} required tables present")
    return True


def _check_meta(
    conn: sqlite3.Connection,
    rep: Report,
    meta: dict[str, str],
) -> tuple[int, int, int]:
    """``M1``/``M2``: schema version and build geometry."""
    version = meta.get("schema_version", "")
    if version not in {"", "3"}:
        rep.fail("M1.schema", f"schema_version = {version!r}, expected 3")
    else:
        rep.ok("M1.schema", f"schema_version = {version or '(unset)'}")

    try:
        top_k = int(meta.get("top_k", "0"))
        chunk_size = int(meta.get("chunk_size", "0"))
    except ValueError:
        top_k = chunk_size = 0
    if top_k <= 0 or chunk_size <= 0:
        rep.fail("M2.geometry", f"top_k={top_k!r} chunk_size={chunk_size!r} unreadable")
        return -1, 0, 0
    rep.ok("M2.geometry", f"top_k={top_k} chunk_size={chunk_size}")

    try:
        declared_n = int(meta.get("entry_count", ""))
    except ValueError:
        declared_n = -1
    return declared_n, top_k, chunk_size


def _check_counts(
    conn: sqlite3.Connection, rep: Report, meta: dict[str, str], declared_n: int
) -> int:
    """``C``: entries / chunks / reflections must agree three ways."""
    entries_n = int(conn.execute("SELECT count(*) FROM entries").fetchone()[0])
    chunks_n = int(conn.execute("SELECT count(*) FROM peak_chunks").fetchone()[0])
    try:
        chunk_size = int(meta.get("chunk_size", "0"))
    except ValueError:
        chunk_size = 0

    if declared_n == entries_n:
        rep.ok("C1.entry_count", f"meta.entry_count = {entries_n:,}")
    else:
        rep.fail(
            "C1.entry_count",
            f"meta.entry_count = {declared_n:,} but entries has {entries_n:,}",
        )

    chunk_rows = int(
        conn.execute("SELECT coalesce(sum(n), 0) FROM peak_chunks").fetchone()[0]
    )
    if chunk_rows == entries_n:
        rep.ok("C2.chunk_rows", f"{chunks_n:,} chunks cover {entries_n:,} entries")
    else:
        rep.fail(
            "C2.chunk_rows",
            f"SUM(peak_chunks.n) = {chunk_rows:,} but entries = {entries_n:,}",
        )

    peaks_entries = int(
        conn.execute("SELECT coalesce(sum(n_peaks), 0) FROM entries").fetchone()[0]
    )
    peaks_payload = (
        int(
            conn.execute(
                "SELECT coalesce(sum(length(all_d)), 0) FROM peak_chunks"
            ).fetchone()[0]
        )
        // 4
    )
    if peaks_entries == peaks_payload:
        rep.ok("C3.reflections", f"{peaks_entries:,} reflections in both tables")
    else:
        rep.fail(
            "C3.reflections",
            f"entries declare {peaks_entries:,} reflections, payload holds {peaks_payload:,}",
        )

    if chunk_size > 0:
        expected = -(-entries_n // chunk_size)
        if chunks_n == expected:
            rep.ok("C4.chunk_size", f"ceil({entries_n:,}/{chunk_size}) = {chunks_n:,}")
        else:
            rep.warn(
                "C4.chunk_size",
                f"{chunks_n:,} chunks, expected {expected:,} for chunk_size={chunk_size}",
            )
    return entries_n


def _check_id_band(
    conn: sqlite3.Connection,
    rep: Report,
    kind: str,
    offset: int,
    entries_n: int,
) -> tuple[int | None, int | None]:
    """``R1``: exhaustive id <-> numeric ownership and source-band check.

    Returns ``(archive_id_min, archive_id_max)`` so the report can state the
    id range it actually saw.
    """
    band_lo, band_hi = BANDS[kind]
    bad = 0
    bad_samples: list[str] = []
    raw_min: int | None = None
    raw_max: int | None = None
    for (entry_id,) in conn.execute("SELECT entry_id FROM entries ORDER BY rowid"):
        stored = _entry_number_for(entry_id, kind)
        problem: str | None = None
        if stored is None:
            expected = "mp-<number>" if kind == "MP" else "96-XXX-XXXX"
            problem = f"{entry_id!r} is not a {expected} id"
        elif kind == "MP" and str(entry_id) != f"mp-{stored}":
            problem = f"{entry_id!r} is not the canonical display form of mp-{stored}"
        elif kind != "MP" and _display_id(stored) != str(entry_id):
            problem = f"{entry_id!r} is not the canonical display form of {stored:,}"
        elif not band_lo <= stored < band_hi:
            problem = (
                f"{entry_id!r} (stored id {stored:,}) is outside the {kind} band "
                f"[{band_lo:,}, {band_hi:,})"
            )
        if problem is not None:
            bad += 1
            if len(bad_samples) < 5:
                bad_samples.append(problem)
            continue
        assert stored is not None
        raw = stored - offset
        raw_min = raw if raw_min is None else min(raw_min, raw)
        raw_max = raw if raw_max is None else max(raw_max, raw)

    rep.stats["id_checked"] = entries_n
    rep.stats["id_raw_min"] = raw_min
    rep.stats["id_raw_max"] = raw_max
    if bad:
        rep.fail(
            "R1.id_band",
            f"{bad:,}/{entries_n:,} ids violate the {kind} ownership rule "
            f"(offset {offset:,}); e.g. {'; '.join(bad_samples)}",
        )
    elif entries_n:
        rep.ok(
            "R1.id_band",
            f"all {entries_n:,} ids map back to {kind} archive ids "
            f"[{raw_min:,}, {raw_max:,}] with offset {offset:,}",
        )
    return raw_min, raw_max


def _open_manifest(rep: Report, manifest: Path | None) -> sqlite3.Connection | None:
    """Open the manifest read-only, warning (not failing) when unusable."""
    if manifest is None:
        return None
    try:
        conn = sqlite3.connect(f"file:{manifest.as_posix()}?mode=ro", uri=True, timeout=10.0)
        conn.execute("SELECT 1 FROM files LIMIT 1")
    except sqlite3.Error as exc:
        rep.warn("A0.manifest", f"cannot use manifest {manifest}: {exc}")
        return None
    return conn


def _check_mp_ownership(
    conn: sqlite3.Connection,
    rep: Report,
    source_root: Path | None,
    entries_n: int,
) -> None:
    """``A1`` for MP: the entry-id set must equal the archive's mp-id set.

    The MP archive is formula-named, so there is no ``id -> path`` layout to
    stat.  Instead the tool re-scans the archive headers (165k small reads) and
    compares the *sets* exactly -- stronger than the sampled check used for the
    COD family.
    """
    if source_root is None:
        rep.warn(
            "A1.ownership",
            "cannot verify MP ids: no archive root (set MP_CIF_ROOT, pass --mp-root, "
            "or keep meta.cif_source alive)",
        )
        return
    try:
        from xinmatch.core.db.mp_scan import scan_mp_archive
    except ImportError as exc:  # pragma: no cover - repo checkout always has it
        rep.warn("A1.ownership", f"cannot import the MP scanner: {exc}")
        return
    try:
        archive = scan_mp_archive(source_root)
    except Exception as exc:
        rep.warn("A1.ownership", f"cannot scan MP archive {source_root}: {exc}")
        return

    archive_ids = {mp_id for mp_id, _path in archive.entries}
    db_ids: set[int] = set()
    bad: list[str] = []
    for (entry_id,) in conn.execute("SELECT entry_id FROM entries ORDER BY rowid"):
        stored = _entry_number_for(entry_id, "MP")
        if stored is None:
            bad.append(str(entry_id))
        else:
            db_ids.add(stored)
    missing = archive_ids - db_ids
    extra = db_ids - archive_ids
    rep.stats["mp_archive_ids"] = len(archive_ids)
    rep.stats["mp_db_ids"] = len(db_ids)
    rep.stats["mp_missing"] = len(missing)
    rep.stats["mp_extra"] = len(extra)
    rep.stats["mp_no_id"] = len(archive.no_id)

    duplicate_n = entries_n - len(db_ids) - len(bad) if entries_n else 0
    if bad or missing or extra or duplicate_n:
        parts = []
        if bad:
            parts.append(f"{len(bad)} unparseable ids (e.g. {bad[:3]})")
        if missing:
            sample = sorted(missing)[:5]
            parts.append(f"{len(missing)} archive ids missing (e.g. {sample})")
        if extra:
            sample = sorted(extra)[:5]
            parts.append(f"{len(extra)} ids absent from the archive (e.g. {sample})")
        if duplicate_n:
            parts.append(f"{duplicate_n} duplicate ids")
        rep.fail(
            "A1.ownership",
            f"MP id set does not match {source_root}: " + "; ".join(parts),
        )
    else:
        rep.ok(
            "A1.ownership",
            f"{len(db_ids):,} mp ids match the archive exactly "
            f"({archive.total:,} CIF; {len(archive.no_id)} without an mp-id header skipped)",
        )


def _check_archive_ownership(
    conn: sqlite3.Connection,
    rep: Report,
    kind: str,
    offset: int,
    entries_n: int,
    dataset: str,
    source_root: Path | None,
    source_origin: str,
    manifest_conn: sqlite3.Connection | None,
    sample_ids: int,
    seed: int,
) -> None:
    """``A1``: sampled ids must resolve to real CIFs (or manifest rows)."""
    if entries_n == 0:
        rep.fail("A1.ownership", "library holds no entries")
        return
    if kind == "MP":
        _check_mp_ownership(conn, rep, source_root, entries_n)
        return
    if source_root is None and manifest_conn is None:
        rep.warn(
            "A1.ownership",
            "cannot verify id ownership: no source root found "
            "(set XRD_DATA_ROOT, pass --manifest, or keep meta.cif_source alive)",
        )
        return

    rng = random.Random(seed)
    offsets = sorted(rng.sample(range(entries_n), min(sample_ids, entries_n)))
    placeholders = ",".join("?" * len(offsets))
    picked = conn.execute(
        f"SELECT rowid, entry_id FROM entries WHERE rowid IN ({placeholders}) ORDER BY rowid",
        offsets,
    ).fetchall()

    checked = 0
    located = 0
    manifest_hits = 0
    manifest_flags = 0
    misses: list[str] = []
    for _rowid, entry_id in picked:
        stored = _entry_number(entry_id)
        if stored is None:
            misses.append(f"{entry_id!r} (unparseable id)")
            continue
        checked += 1
        archive_id = stored - offset
        located_here = False
        if source_root is not None:
            for rel in _source_cif_candidates(archive_id):
                if (source_root / rel).is_file():
                    located_here = True
                    break
        manifest_rel: str | None = None
        if manifest_conn is not None:
            row = manifest_conn.execute(
                "SELECT rel, is_inorg FROM files WHERE cod_id = ? AND lib = ?",
                (stored, kind),
            ).fetchone()
            if row is not None:
                manifest_hits += 1
                manifest_rel = str(row[0] or "")
                if dataset == "inorg" and not int(row[1]):
                    manifest_flags += 1
                if source_root is None:
                    # No CIF tree to stat; a manifest row is the ownership
                    # evidence (the manifest is a scan of the real archives).
                    located_here = True
                elif manifest_rel:
                    rel = manifest_rel.replace("\\", "/")
                    if (source_root / rel).is_file():
                        located_here = True
                    elif not located_here:
                        # The manifest points at a path the archive lacks:
                        # ownership is still evidenced by the row, but the
                        # disagreement is worth a warning.
                        manifest_flags += 1
        if located_here:
            located += 1
        else:
            misses.append(
                f"{entry_id}({kind} id {archive_id} -> "
                f"{' / '.join(_source_cif_candidates(archive_id)) or 'no layout'})"
            )

    rep.stats["id_sampled"] = checked
    rep.stats["id_located"] = located
    rep.stats["manifest_hits"] = manifest_hits
    if manifest_flags:
        rep.warn(
            "A2.manifest_flags",
            f"{manifest_flags} sampled ids disagree with the manifest "
            f"(is_inorg=0 for an {dataset} library, or a manifest path the source tree "
            f"does not contain)",
        )
    if misses:
        rep.fail(
            "A1.ownership",
            f"{len(misses)}/{checked or len(picked)} sampled ids do not resolve to a source "
            f"CIF (root: {source_root or 'none'}; origin: {source_origin}); "
            f"e.g. {', '.join(misses[:5])}",
        )
    else:
        extra = ""
        if manifest_hits:
            extra = f"; {manifest_hits} confirmed by the manifest"
        where = str(source_root) if source_root is not None else "the manifest"
        rep.ok(
            "A1.ownership",
            f"{located}/{checked} sampled ids resolve to a source CIF "
            f"({where}; origin: {source_origin}){extra}",
        )


def _check_element_index(
    conn: sqlite3.Connection,
    rep: Report,
    entries_n: int,
    sample_elements: int,
    seed: int,
) -> None:
    """``E``: full symbol-set equality, eid range/order, sampled membership."""
    used: set[str] = set()
    declared_pairs = 0
    for (raw,) in conn.execute("SELECT elements FROM entries"):
        symbols_here = _split_elements(raw)
        declared_pairs += len(symbols_here)
        for symbol in symbols_here:
            used.add(symbol.title())
    index_row = conn.execute("SELECT element, eids FROM element_inv ORDER BY element").fetchall()
    indexed = {str(symbol) for symbol, _blob in index_row}
    rep.stats["elements"] = len(indexed)

    missing = sorted(used - indexed)
    extra = sorted(indexed - used)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"{len(missing)} used but not indexed ({', '.join(missing[:8])})")
        if extra:
            parts.append(f"{len(extra)} indexed but not used ({', '.join(extra[:8])})")
        rep.fail(
            "E1.element_sets",
            "element_inv disagrees with entries.elements: " + "; ".join(parts),
        )
    else:
        rep.ok(
            "E1.element_sets",
            f"element_inv symbol set equals entries.elements set ({len(indexed)} elements)",
        )

    indexed_pairs = sum(int(_decode_eids(blob).size) for _symbol, blob in index_row)
    if indexed_pairs == declared_pairs:
        rep.ok(
            "E5.pairs",
            f"element_inv holds exactly the {declared_pairs:,} declared (element, phase) pairs",
        )
    else:
        rep.fail(
            "E5.pairs",
            f"element_inv holds {indexed_pairs:,} (element, phase) pairs but "
            f"entries.elements declares {declared_pairs:,}",
        )

    bad_range = 0
    bad_order = 0
    for symbol, blob in index_row:
        arr = _decode_eids(blob)
        if arr.size == 0:
            rep.warn("E2.eids", f"element {symbol} has an empty eid list")
            continue
        if int(arr.max()) >= entries_n:
            bad_range += 1
            if bad_range <= 5:
                rep.add(
                    "E2.eids",
                    "FAIL",
                    f"element {symbol}: eid {int(arr.max()):,} >= entry count {entries_n:,}",
                )
        if arr.size > 1 and np.any(np.diff(arr.astype(np.int64)) <= 0):
            bad_order += 1
    if bad_range == 0:
        rep.ok("E2.eids", "every eid is a valid entry index")
    if bad_order:
        rep.warn("E3.eids_order", f"{bad_order} element lists are not strictly ascending")

    if not indexed or entries_n == 0:
        return
    rng = random.Random(seed + 1)
    symbols = sorted(indexed)
    chosen = rng.sample(symbols, min(sample_elements, len(symbols)))
    misses = 0
    checked = 0
    for symbol in chosen:
        row = conn.execute("SELECT eids FROM element_inv WHERE element = ?", (symbol,)).fetchone()
        arr = _decode_eids(row[0]) if row else np.zeros(0, dtype="<u4")
        if arr.size == 0:
            continue
        all_eids = [int(value) for value in arr.tolist()]
        eids = rng.sample(all_eids, min(8, len(all_eids)))
        statement = ",".join("?" * len(eids))
        rows = conn.execute(
            f"SELECT rowid, elements FROM entries WHERE rowid IN ({statement})", eids
        ).fetchall()
        checked += len(rows)
        for rowid, elements in rows:
            symbols_here = {token.title() for token in _split_elements(elements)}
            if symbol not in symbols_here:
                misses += 1
                if misses <= 5:
                    rep.add(
                        "E4.membership",
                        "FAIL",
                        f"element_inv[{symbol}] lists entry {rowid:,} but entries.elements "
                        f"is {elements!r}",
                    )
    rep.stats["element_memberships_checked"] = checked
    if misses:
        rep.add("E4.membership", "FAIL", f"{misses} sampled element/phase memberships do not hold")
    else:
        rep.ok(
            "E4.membership",
            f"{checked} sampled (element, phase) memberships hold across {len(chosen)} elements",
        )


def _check_payload(
    conn: sqlite3.Connection,
    rep: Report,
    meta: dict[str, str],
    top_k: int,
    chunk_size: int,
    sample_chunks: int | None,
    seed: int,
) -> None:
    """``P``: chunk payload shapes, values, normalisation and subset rule."""
    norm_scale = _norm_scale_from_meta(meta)
    ids = [int(r[0]) for r in conn.execute("SELECT chunk_id FROM peak_chunks ORDER BY chunk_id")]
    if not ids:
        rep.fail("P1.payload", "no chunks to check")
        return
    rng = random.Random(seed)
    if sample_chunks is None or sample_chunks <= 0 or sample_chunks >= len(ids):
        picked = sorted(ids)
        rep.info("P0.scope", f"decoded all {len(picked):,} chunks")
    else:
        picked = sorted(rng.sample(ids, sample_chunks))
        rep.info("P0.scope", f"decoded {len(picked):,} of {len(ids):,} chunks")

    problems: list[str] = []
    entries_checked = 0
    total = len(picked)
    for step, cid in enumerate(picked, 1):
        if total > 50 and step % 25 == 0:
            print(f"  payload: {step}/{total} chunks", file=sys.stderr, flush=True)
        row = conn.execute(
            "SELECT n, top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?",
            (cid,),
        ).fetchone()
        if row is None:
            problems.append(f"chunk {cid}: row missing")
            continue
        n = int(row[0])
        top_d_raw = b"" if row[1] is None else bytes(row[1])
        top_i_raw = b"" if row[2] is None else bytes(row[2])
        all_d_raw = b"" if row[3] is None else bytes(row[3])
        all_i_raw = b"" if row[4] is None else bytes(row[4])
        expected_top = n * top_k * 4
        if len(top_d_raw) != expected_top or len(top_i_raw) != expected_top:
            problems.append(
                f"chunk {cid}: prescreen block is {len(top_d_raw)}/{len(top_i_raw)} bytes, "
                f"expected {expected_top}"
            )
            continue
        if len(all_d_raw) != len(all_i_raw):
            problems.append(
                f"chunk {cid}: all_d/all_i are {len(all_d_raw)} vs {len(all_i_raw)} bytes"
            )
            continue
        if len(all_d_raw) % 4:
            problems.append(f"chunk {cid}: all_d payload is not a whole number of float32")
            continue

        top_d = _decode_matrix(top_d_raw, n, top_k)
        top_i = _decode_matrix(top_i_raw, n, top_k)
        all_d = _decode_flat(all_d_raw)
        all_i = _decode_flat(all_i_raw)
        if not np.all(np.isfinite(top_d)) or not np.all(np.isfinite(top_i)):
            problems.append(f"chunk {cid}: non-finite value in the prescreen block")
        if not np.all(np.isfinite(all_d)) or not np.all(np.isfinite(all_i)):
            problems.append(f"chunk {cid}: non-finite value in the full peak block")
        if top_d.size and float(top_d.min()) < 0:
            problems.append(f"chunk {cid}: negative d-spacing in the prescreen block")
        if all_d.size and float(all_d.min()) <= 0:
            problems.append(f"chunk {cid}: non-positive d-spacing in the full block")
        if all_i.size and float(all_i.min()) < 0:
            problems.append(f"chunk {cid}: negative intensity in the full block")

        entries = conn.execute(
            "SELECT entry_id, n_peaks, peak_off FROM entries "
            "WHERE rowid >= ? AND rowid < ? ORDER BY rowid",
            (cid * chunk_size, cid * chunk_size + n),
        ).fetchall()
        if len(entries) != n:
            problems.append(
                f"chunk {cid}: {len(entries)} entries at that offset, chunk says {n}"
            )
            continue
        chunk_peak_base = int(entries[0][2])

        for index, (entry_id, n_peaks, peak_off) in enumerate(entries):
            n_peaks = int(n_peaks)
            entries_checked += 1
            if n_peaks <= 0:
                continue
            lo = int(peak_off) - chunk_peak_base
            if lo < 0 or lo + n_peaks > len(all_d):
                problems.append(
                    f"{entry_id}: slice [{lo}, {lo + n_peaks}) outside payload of {len(all_d)}"
                )
                continue
            d_slice = all_d[lo : lo + n_peaks]
            i_slice = all_i[lo : lo + n_peaks]
            if float(d_slice.min()) <= 0:
                problems.append(f"{entry_id}: non-positive d-spacing in its slice")
            if float(i_slice.min()) < 0:
                problems.append(f"{entry_id}: negative intensity in its slice")
            strongest = float(i_slice.max())
            if strongest > norm_scale + NORM_TOL:
                problems.append(
                    f"{entry_id}: strongest reflection {strongest:.1f} above norm scale "
                    f"{norm_scale:.0f}"
                )
            elif abs(strongest - norm_scale) > NORM_TOL:
                problems.append(
                    f"{entry_id}: strongest reflection {strongest:.3f} != {norm_scale:.0f}"
                )
            take = min(top_k, n_peaks)
            if take:
                if not np.all(
                    np.isin(np.round(top_d[index, :take], 5), np.round(d_slice, 5))
                ):
                    problems.append(
                        f"{entry_id}: prescreen d values absent from the full reflection list"
                    )
                if not np.all(
                    np.isin(np.round(top_i[index, :take], 3), np.round(i_slice, 3))
                ):
                    problems.append(
                        f"{entry_id}: prescreen intensities absent from the full list"
                    )
            pad_d = top_d[index, take:]
            pad_i = top_i[index, take:]
            if pad_d.size and (np.any(pad_d != 0.0) or np.any(pad_i != 0.0)):
                problems.append(f"{entry_id}: prescreen padding is not zero")
            if len(problems) >= MAX_PAYLOAD_PROBLEMS:
                break
        if len(problems) >= MAX_PAYLOAD_PROBLEMS:
            break

    rep.stats["chunks_checked"] = len(picked)
    rep.stats["entries_payload_checked"] = entries_checked
    if problems:
        for message in problems[:10]:
            rep.add("P1.payload", "FAIL", message)
        if len(problems) > 10:
            rep.add("P1.payload", "FAIL", f"... {len(problems) - 10} more payload problems")
    else:
        rep.ok(
            "P1.payload",
            f"{len(picked):,} chunks / {entries_checked:,} entries decode cleanly on the "
            f"norm-{norm_scale:.0f} scale",
        )


def _check_bad_phases(
    conn: sqlite3.Connection,
    rep: Report,
    meta: dict[str, str],
    chunk_size: int,
) -> None:
    """``B1``/``B2``: exhaustive bad-phase scan.

    ``B1`` is anchored on the **strongest reflection** of every phase
    (``peak_chunks.top_d`` column 0; rows are intensity-ordered) so it also
    catches a pattern that only partially collapsed -- e.g. PCOD
    ``96-254-30089`` (n_peaks = 2: the bogus 180 deg line at 1000 plus one real
    line at 0.013), which an ``n_peaks == 1`` test misses.

    A 180 deg boundary artifact is only flagged when the rest of the pattern
    has collapsed too: the Lorentz-polarisation factor diverges near
    2theta = 180 deg, so any *real* reflection there also becomes the strongest
    line, and "strongest d at the cutoff" alone would flag thousands of healthy
    phases.  The placeholder-cell signature (a lone 1.0 A line) is unchanged.
    """
    cutoff = _cutoff_from_meta(meta)
    try:
        top_k = int(meta.get("top_k", "0"))
    except ValueError:
        top_k = 0
    stats = rep.stats

    # -- B2: a phase that stores no reflection at all -----------------------
    empty = [
        str(row[0])
        for row in conn.execute(
            "SELECT entry_id FROM entries WHERE n_peaks IS NULL OR n_peaks <= 0 ORDER BY rowid"
        )
    ]
    if empty:
        rep.fail(
            "B2.peak_count",
            f"{len(empty):,} phases store no reflection (n_peaks <= 0); "
            f"e.g. {', '.join(empty[:MAX_QUOTED_IDS])}",
        )
    else:
        rep.ok("B2.peak_count", "every phase stores at least one reflection")

    bad: list[dict[str, object]] = []
    semi_collapsed: list[dict[str, object]] = []
    # -- B1a: strongest reflection at the 2theta = 180 deg cutoff -----------
    if chunk_size > 0 and top_k > 0:
        for chunk_id, n, top_d_blob, top_i_blob in conn.execute(
            "SELECT chunk_id, n, top_d, top_i FROM peak_chunks ORDER BY chunk_id"
        ):
            n = int(n)
            top_d = _decode_matrix(top_d_blob, n, top_k)
            top_i = _decode_matrix(top_i_blob, n, top_k)
            strongest = top_d[:, 0]
            rest = top_i[:, 1:].sum(axis=1)
            window = (strongest >= cutoff - DMIN_LOW_SLACK) & (
                strongest <= cutoff + DMIN_HIGH_SLACK
            )
            if not window.any():
                continue
            collapsed_mask = window & (rest < BOUNDARY_COLLAPSE_MAX)
            semi_mask = window & ~collapsed_mask
            rows = conn.execute(
                "SELECT entry_id, n_peaks FROM entries WHERE rowid >= ? AND rowid < ? "
                "ORDER BY rowid",
                (int(chunk_id) * chunk_size, int(chunk_id) * chunk_size + n),
            ).fetchall()
            if len(rows) != n:
                rep.warn(
                    "B1.single_peak",
                    f"chunk {int(chunk_id)}: {len(rows)} entries at that offset, "
                    f"chunk says {n} (see C2)",
                )
                continue
            for index in np.flatnonzero(collapsed_mask):
                entry_id, n_peaks = rows[int(index)]
                bad.append(
                    {
                        "entry_id": str(entry_id),
                        "d": round(float(strongest[int(index)]), 5),
                        "n_peaks": int(n_peaks),
                        "signature": "boundary",
                    }
                )
            for index in np.flatnonzero(semi_mask):
                entry_id, n_peaks = rows[int(index)]
                semi_collapsed.append(
                    {
                        "entry_id": str(entry_id),
                        "d": round(float(strongest[int(index)]), 5),
                        "n_peaks": int(n_peaks),
                        "rest_sum": round(float(rest[int(index)]), 3),
                    }
                )
    boundary_ids = {item["entry_id"] for item in bad}

    # -- B1b/B1c: lone reflections (1.0 A placeholder vs plausible) ---------
    singles = [
        (int(rowid), str(entry_id), int(peak_off))
        for rowid, entry_id, peak_off in conn.execute(
            "SELECT rowid, entry_id, peak_off FROM entries WHERE n_peaks = 1 ORDER BY rowid"
        )
    ]
    normal: list[str] = []
    undecodable: list[str] = []
    if singles and chunk_size > 0:
        cache: dict[int, tuple[int, bytes]] = {}
        for rowid, entry_id, peak_off in singles:
            cid = rowid // chunk_size
            cached = cache.get(cid)
            if cached is None:
                chunk_row = conn.execute(
                    "SELECT n, all_d FROM peak_chunks WHERE chunk_id = ?", (cid,)
                ).fetchone()
                if chunk_row is None:
                    cached = (-1, b"")
                else:
                    base_row = conn.execute(
                        "SELECT peak_off FROM entries WHERE rowid = ?",
                        (cid * chunk_size,),
                    ).fetchone()
                    cached = (
                        int(base_row[0]) if base_row is not None else -1,
                        bytes(chunk_row[1] or b""),
                    )
                cache[cid] = cached
            base, payload = cached
            lo = int(peak_off) - base
            if base < 0 or lo < 0 or (lo + 1) * 4 > len(payload):
                undecodable.append(entry_id)
                continue
            value = float(np.frombuffer(payload[lo * 4 : (lo + 1) * 4], dtype="<f4")[0])
            if entry_id in boundary_ids:
                continue  # already reported by the strongest-reflection scan
            if abs(value - PLACEHOLDER_D) <= PLACEHOLDER_TOL:
                bad.append(
                    {
                        "entry_id": entry_id,
                        "d": round(value, 5),
                        "n_peaks": 1,
                        "signature": "placeholder",
                    }
                )
            else:
                normal.append(f"{entry_id}(d={value:.4f})")

    stats["single_peak_phases"] = len(singles)
    stats["bad_phase_ids"] = bad[:MAX_STATS_IDS]
    stats["boundary_phases"] = sum(1 for item in bad if item["signature"] == "boundary")
    stats["placeholder_phases"] = sum(1 for item in bad if item["signature"] == "placeholder")
    stats["normal_single_peak_phases"] = len(normal)
    # Semi-collapsed residue: strongest reflection at the cutoff but a populated
    # pattern.  Same root cause as the collapse artifact, but a real reflection
    # near 2theta = 180 deg is legitimately the strongest line there, so this is
    # informational only (WARN); the total is in stats, the closest-to-collapse
    # cases are quoted.
    semi_collapsed.sort(key=lambda item: float(item["rest_sum"]))
    stats["boundary_warn_phases"] = len(semi_collapsed)
    stats["boundary_warn_ids"] = semi_collapsed[:MAX_STATS_IDS]

    if bad:
        quoted = ", ".join(
            f"{item['entry_id']}(d={item['d']}, n={item['n_peaks']}, {item['signature']})"
            for item in bad[:MAX_QUOTED_IDS]
        )
        more = "" if len(bad) <= MAX_QUOTED_IDS else f" ... {len(bad) - MAX_QUOTED_IDS} more"
        rep.fail(
            "B1.single_peak",
            f"{len(bad)} bogus phases "
            f"({stats['boundary_phases']} with the strongest reflection at the {cutoff} A "
            f"cutoff and a collapsed pattern, {stats['placeholder_phases']} with a lone "
            f"{PLACEHOLDER_D} A placeholder line): {quoted}{more}",
        )
    else:
        rep.ok("B1.single_peak", "no bogus boundary/placeholder phase")
    if normal:
        rep.warn(
            "B1.single_peak",
            f"{len(normal)} phases hold exactly one reflection at a plausible d-spacing: "
            f"{', '.join(normal[:MAX_QUOTED_IDS])}",
        )
    if undecodable:
        rep.warn(
            "B1.single_peak",
            f"{len(undecodable)} lone reflections could not be decoded: "
            f"{', '.join(undecodable[:MAX_QUOTED_IDS])}",
        )
    if semi_collapsed:
        samples = ", ".join(
            f"{item['entry_id']}(d={item['d']}, n={item['n_peaks']}, "
            f"rest_sum={item['rest_sum']})"
            for item in semi_collapsed[:MAX_BOUNDARY_WARN_QUOTED]
        )
        rep.warn(
            "B1.boundary_warn",
            f"{len(semi_collapsed)} phases have their strongest reflection in the "
            f"{cutoff} A cutoff window but a still-populated pattern, so they are NOT "
            f"bad phases (a real reflection near 2theta=180 deg is legitimately the "
            f"strongest line there); closest to the collapse boundary: {samples}",
        )


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
def accept(
    db: Path,
    kind: str = "auto",
    manifest: Path | None = None,
    *,
    data_root: Path | None = None,
    mp_root: Path | None = None,
    sample_ids: int = 40,
    sample_elements: int = 20,
    sample_chunks: int | None = None,
    seed: int = 0,
) -> Report:
    """Run the acceptance checks against *db* and return the populated report."""
    db = Path(db)
    rep = Report(db=db)
    t0 = time.time()
    if not db.exists():
        rep.fail("F0.exists", f"no such file: {db}")
        rep.stats["elapsed_sec"] = 0.0
        return rep
    size_mb = db.stat().st_size / (1024.0 * 1024.0)
    rep.info("F1.file", f"{db.name}  {size_mb:,.1f} MB")

    uri = f"file:{db.as_posix()}?mode=ro"
    try:
        conn = sqlite3.connect(uri, uri=True, timeout=30.0)
    except sqlite3.Error as exc:
        rep.fail("F0.exists", f"cannot open: {exc}")
        rep.stats["elapsed_sec"] = round(time.time() - t0, 2)
        return rep

    manifest_conn: sqlite3.Connection | None = None
    try:
        if not _check_schema(conn, rep):
            return rep
        meta = _read_meta(conn)
        derived_kind = _derive_kind(db, meta, kind)
        if derived_kind is None:
            rep.fail(
                "M0.kind",
                "cannot derive the source archive from meta.library_kind / file name; "
                "pass --kind COD|PCOD|TCOD|MP",
            )
            return rep
        rep.kind = derived_kind
        dataset = _derive_dataset(meta)
        rep.dataset = dataset
        if rep.kind == "MP":
            offset = 0
            root = None
            origin = "not found"
            candidates = [
                (mp_root, "--mp-root" if mp_root is not None else ""),
                (Path(meta.get("cif_source", "")) if meta.get("cif_source") else None, "meta.cif_source"),
                (DEFAULT_MP_ROOT, "DEFAULT_MP_ROOT"),
            ]
            for candidate, candidate_origin in candidates:
                if candidate is not None and candidate.is_dir():
                    root, origin = candidate, candidate_origin
                    break
        else:
            offset = SOURCES[rep.kind][0]
            root, origin = _source_root(
                meta, rep.kind, dataset, data_root if data_root is not None else DEFAULT_DATA_ROOT
            )
        rep.stats["dataset"] = dataset
        rep.stats["id_offset"] = offset
        rep.stats["source_root"] = str(root) if root else None
        rep.stats["source_origin"] = origin

        declared_n, top_k, chunk_size = _check_meta(conn, rep, meta)
        if top_k <= 0 or chunk_size <= 0:
            return rep
        entries_n = _check_counts(conn, rep, meta, declared_n)
        rep.stats["entries"] = entries_n
        rep.stats["chunks"] = int(
            conn.execute("SELECT count(*) FROM peak_chunks").fetchone()[0]
        )

        _check_id_band(conn, rep, rep.kind, offset, entries_n)
        manifest_conn = _open_manifest(rep, Path(manifest) if manifest is not None else None)
        _check_archive_ownership(
            conn,
            rep,
            rep.kind,
            offset,
            entries_n,
            dataset,
            root,
            origin,
            manifest_conn,
            sample_ids,
            seed,
        )
        _check_element_index(conn, rep, entries_n, sample_elements, seed)
        _check_payload(conn, rep, meta, top_k, chunk_size, sample_chunks, seed)
        _check_bad_phases(conn, rep, meta, chunk_size)
    except sqlite3.Error as exc:
        rep.fail("F2.sqlite", f"{type(exc).__name__}: {exc}")
    finally:
        if manifest_conn is not None:
            manifest_conn.close()
        conn.close()

    rep.stats["elapsed_sec"] = round(time.time() - t0, 2)
    return rep


def main(argv: list[str] | None = None) -> int:
    """Command line entry point; returns the process exit code."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("db", type=Path, help=".xinDB file to accept (opened read-only)")
    parser.add_argument(
        "--kind",
        choices=["auto", "COD", "PCOD", "TCOD", "MP"],
        default="auto",
        help="source archive; auto reads meta.library_kind / the file name",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=None,
        help="inorg_manifest.sqlite (fallback ownership witness + inorganic flag check)",
    )
    parser.add_argument(
        "--sample-ids",
        type=int,
        default=40,
        help="entries sampled for archive ownership evidence (default 40)",
    )
    parser.add_argument(
        "--sample-elements",
        type=int,
        default=20,
        help="elements sampled for inverted-index membership (default 20)",
    )
    parser.add_argument(
        "--sample-chunks",
        type=int,
        default=None,
        help="peak chunks decoded (default: all; a positive value samples that many)",
    )
    parser.add_argument("--seed", type=int, default=0, help="sampling seed")
    parser.add_argument("--data-root", type=Path, default=None, help="CIF archive root override")
    parser.add_argument(
        "--mp-root",
        type=Path,
        default=None,
        help="Materials Project CIF directory (default: meta.cif_source / $MP_CIF_ROOT)",
    )
    parser.add_argument("--json", type=Path, default=None, help="write the report as JSON here")
    args = parser.parse_args(argv)

    rep = accept(
        args.db,
        kind=args.kind,
        manifest=args.manifest,
        data_root=args.data_root,
        mp_root=args.mp_root,
        sample_ids=args.sample_ids,
        sample_elements=args.sample_elements,
        sample_chunks=args.sample_chunks,
        seed=args.seed,
    )
    print(f"accept_xindb: {args.db}")
    print(f"  kind={rep.kind} dataset={rep.dataset}")
    print(rep.render())
    stats = rep.stats
    if stats:
        print(
            f"\n  entries={stats.get('entries', '?')} "
            f"chunks={stats.get('chunks', '?')} "
            f"id_range=[{stats.get('id_raw_min', '?')}, {stats.get('id_raw_max', '?')}] "
            f"elements={stats.get('elements', '?')}"
        )
    print(
        f"\n  {rep.verdict}: {rep.fails} fail, {rep.warns} warn, "
        f"{len(rep.rows)} checks ({stats.get('elapsed_sec', 0)}s)"
    )

    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(
            json.dumps(rep.as_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )
        print(f"  json -> {args.json}")
    return 1 if rep.fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
