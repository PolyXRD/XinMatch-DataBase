"""Tests for the upgraded ``tools/validate_xindb.py`` acceptance checker.

What is covered
---------------
* a real 12-phase build (from the ``tiny_archive`` fixture) passes every check;
* the CLI/JSON contract is unchanged: ``db/verdict/fails/warns/checks`` keys,
  exit code 0 on PASS and 1 on FAIL, ``--json`` is parseable;
* ``--all-chunks`` decodes every block (and says so) instead of sampling 8;
* ``E5.boundary_peak`` distinguishes the two known bad-phase signatures -- a
  lone reflection at the d-spacing cutoff (2theta = 180 deg artifact) and a
  lone 1.0 A reflection (placeholder-cell artifact) are FAILs, while a plausible
  lone reflection is only a WARN;
* ``E6.peak_count`` catches a phase that stores no reflection at all.

Every library used here is built under pytest's ``tmp_path``; the published
libraries in ``data/`` and the real CIF archives are never touched.
"""

from __future__ import annotations

import json
import sqlite3
import struct
from pathlib import Path

import pytest

from tools.validate_xindb import main as validate_main
from tools.validate_xindb import validate

#: Snapshot date used by every build, so file names satisfy ``M3.identity``.
DATA_DATE = "20261020"

#: Number of phases in the ``tiny_archive`` fixture.
N_PHASES = 12


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _build_small(
    cif_root: Path,
    out: Path,
    *,
    library_kind: str = "COD-Inorg",
    id_offset: int = 0,
) -> Path:
    """Build a 12-phase library from *cif_root* through the real builder."""
    from tools.build_cod_xindb import build

    build(
        cif_root,
        out,
        library_kind=library_kind,
        data_date=DATA_DATE,
        id_offset=id_offset,
        limit=N_PHASES,
        workers=1,
        log_every=1000,
    )
    assert out.is_file(), f"build did not produce {out}"
    return out


def _row(report, code: str) -> list[tuple[str, str]]:
    """Return ``(status, detail)`` for every row with *code*."""
    return [(status, detail) for row_code, status, detail in report.rows if row_code == code]


def _mutate_to_lone_peak(db: Path, d: float) -> str:
    """Rewrite the **last** phase of *db* into a single-reflection phase.

    Picking the last row keeps the global ``peak_off`` ladder valid: only its
    own reflections are dropped, so the ``C`` accounting checks still pass and
    the test isolates the ``E5``/``E6`` behaviour under test.
    """
    conn = sqlite3.connect(str(db))
    try:
        meta = dict(conn.execute("SELECT key, value FROM meta"))
        chunk_size = int(meta["chunk_size"])
        top_k = int(meta["top_k"])
        rowid = int(conn.execute("SELECT max(rowid) FROM entries").fetchone()[0])
        entry_id, peak_off = conn.execute(
            "SELECT entry_id, peak_off FROM entries WHERE rowid = ?", (rowid,)
        ).fetchone()
        chunk_id = rowid // chunk_size
        row_in_chunk = rowid - chunk_id * chunk_size
        top_d, top_i, all_d, all_i = conn.execute(
            "SELECT top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?",
            (chunk_id,),
        ).fetchone()
        chunk_base = int(
            conn.execute(
                "SELECT peak_off FROM entries WHERE rowid = ?", (chunk_id * chunk_size,)
            ).fetchone()[0]
        )
        offset = int(peak_off) - chunk_base
        new_all_d = bytes(all_d)[: offset * 4] + struct.pack("<f", d)
        new_all_i = bytes(all_i)[: offset * 4] + struct.pack("<f", 1000.0)
        row_bytes = top_k * 4
        new_top_d = bytearray(bytes(top_d))
        new_top_i = bytearray(bytes(top_i))
        new_top_d[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        new_top_i[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        new_top_d[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 4] = struct.pack("<f", d)
        new_top_i[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 4] = struct.pack(
            "<f", 1000.0
        )
        conn.execute("UPDATE entries SET n_peaks = 1, peak_n = 1 WHERE rowid = ?", (rowid,))
        conn.execute(
            "UPDATE peak_chunks SET top_d = ?, top_i = ?, all_d = ?, all_i = ? "
            "WHERE chunk_id = ?",
            (bytes(new_top_d), bytes(new_top_i), new_all_d, new_all_i, chunk_id),
        )
        conn.commit()
        return str(entry_id)
    finally:
        conn.close()


def _mutate_to_empty_phase(db: Path) -> str:
    """Rewrite the last phase of *db* so that it stores zero reflections."""
    conn = sqlite3.connect(str(db))
    try:
        meta = dict(conn.execute("SELECT key, value FROM meta"))
        chunk_size = int(meta["chunk_size"])
        top_k = int(meta["top_k"])
        rowid = int(conn.execute("SELECT max(rowid) FROM entries").fetchone()[0])
        entry_id, n_peaks, peak_off = conn.execute(
            "SELECT entry_id, n_peaks, peak_off FROM entries WHERE rowid = ?", (rowid,)
        ).fetchone()
        chunk_id = rowid // chunk_size
        row_in_chunk = rowid - chunk_id * chunk_size
        top_d, top_i, all_d, all_i = conn.execute(
            "SELECT top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?",
            (chunk_id,),
        ).fetchone()
        chunk_base = int(
            conn.execute(
                "SELECT peak_off FROM entries WHERE rowid = ?", (chunk_id * chunk_size,)
            ).fetchone()[0]
        )
        offset = int(peak_off) - chunk_base
        new_all_d = bytes(all_d)[: offset * 4]
        new_all_i = bytes(all_i)[: offset * 4]
        row_bytes = top_k * 4
        new_top_d = bytearray(bytes(top_d))
        new_top_i = bytearray(bytes(top_i))
        new_top_d[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        new_top_i[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        conn.execute("UPDATE entries SET n_peaks = 0, peak_n = 0 WHERE rowid = ?", (rowid,))
        conn.execute(
            "UPDATE peak_chunks SET top_d = ?, top_i = ?, all_d = ?, all_i = ? "
            "WHERE chunk_id = ?",
            (bytes(new_top_d), bytes(new_top_i), new_all_d, new_all_i, chunk_id),
        )
        conn.commit()
        return str(entry_id)
    finally:
        conn.close()


def _mutate_to_collapsed_double_peak(
    db: Path,
    *,
    boundary_d: float = 0.7703,
    real_d: float = 2.5,
    real_i: float = 0.0127,
) -> str:
    """Rewrite the last phase of *db* into a **two-line collapsed** pattern.

    This is the PCOD ``96-254-30089`` signature: the bogus 2theta = 180 deg line
    carries the whole normalised intensity (1000) while the one real reflection
    that survived the normalization is tiny (``real_i``).  An ``n_peaks == 1``
    check cannot see it, so ``E5`` must be anchored on the strongest reflection.
    """
    conn = sqlite3.connect(str(db))
    try:
        meta = dict(conn.execute("SELECT key, value FROM meta"))
        chunk_size = int(meta["chunk_size"])
        top_k = int(meta["top_k"])
        rowid = int(conn.execute("SELECT max(rowid) FROM entries").fetchone()[0])
        entry_id, peak_off = conn.execute(
            "SELECT entry_id, peak_off FROM entries WHERE rowid = ?", (rowid,)
        ).fetchone()
        chunk_id = rowid // chunk_size
        row_in_chunk = rowid - chunk_id * chunk_size
        top_d, top_i, all_d, all_i = conn.execute(
            "SELECT top_d, top_i, all_d, all_i FROM peak_chunks WHERE chunk_id = ?",
            (chunk_id,),
        ).fetchone()
        chunk_base = int(
            conn.execute(
                "SELECT peak_off FROM entries WHERE rowid = ?", (chunk_id * chunk_size,)
            ).fetchone()[0]
        )
        offset = int(peak_off) - chunk_base
        row_bytes = top_k * 4
        new_top_d = bytearray(bytes(top_d))
        new_top_i = bytearray(bytes(top_i))
        new_top_d[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        new_top_i[row_in_chunk * row_bytes : (row_in_chunk + 1) * row_bytes] = b"\x00" * row_bytes
        # Intensity order: strongest (boundary) first, the one real line second.
        new_top_d[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 8] = (
            struct.pack("<f", boundary_d) + struct.pack("<f", real_d)
        )
        new_top_i[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 8] = (
            struct.pack("<f", 1000.0) + struct.pack("<f", real_i)
        )
        new_all_d = (
            bytes(all_d)[: offset * 4] + struct.pack("<f", real_d) + struct.pack("<f", boundary_d)
        )
        new_all_i = (
            bytes(all_i)[: offset * 4] + struct.pack("<f", real_i) + struct.pack("<f", 1000.0)
        )
        conn.execute("UPDATE entries SET n_peaks = 2, peak_n = 2 WHERE rowid = ?", (rowid,))
        conn.execute(
            "UPDATE peak_chunks SET top_d = ?, top_i = ?, all_d = ?, all_i = ? "
            "WHERE chunk_id = ?",
            (bytes(new_top_d), bytes(new_top_i), new_all_d, new_all_i, chunk_id),
        )
        conn.commit()
        return str(entry_id)
    finally:
        conn.close()


@pytest.fixture
def good_db(tmp_path: Path, tiny_archive: Path) -> Path:
    """A complete 12-phase ``COD-Inorg`` slice built from ``tiny_archive``."""
    return _build_small(tiny_archive, tmp_path / f"COD-Inorg-{DATA_DATE}.xinDB")


# ---------------------------------------------------------------------------
# happy path and CLI contract
# ---------------------------------------------------------------------------
def test_good_library_passes(good_db: Path) -> None:
    report = validate(good_db, -1, N_PHASES)
    assert report.fails == 0, report.render()
    assert _row(report, "E5.boundary_peak")[0][0] == "PASS"
    assert _row(report, "E6.peak_count")[0][0] == "PASS"
    assert _row(report, "P1.payload")[0][0] == "PASS"
    # The payload check really decoded the only chunk of a 12-phase library.
    assert "12 entries" in _row(report, "P1.payload")[0][1]


def test_cli_json_contract_is_unchanged(good_db: Path, tmp_path: Path, capsys) -> None:
    out = tmp_path / "report.json"
    code = validate_main(
        [str(good_db), "--sample-chunks", "1", "--sample-entries", str(N_PHASES), "--json", str(out)]
    )
    assert code == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert set(payload) == {"db", "verdict", "fails", "warns", "checks"}
    assert payload["verdict"] == "PASSED"
    assert payload["fails"] == 0
    assert payload["warns"] == 0
    assert payload["db"] == str(good_db)
    assert all(set(check) == {"code", "status", "detail"} for check in payload["checks"])
    assert "PASSED" in capsys.readouterr().out


def test_all_chunks_decodes_every_chunk(good_db: Path, tmp_path: Path) -> None:
    out = tmp_path / "all_chunks.json"
    code = validate_main([str(good_db), "--all-chunks", "--json", str(out)])
    assert code == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    scope = [
        check for check in payload["checks"] if check["code"] == "P1.all_chunks"
    ]
    assert scope, "the --all-chunks report must state how many chunks it decoded"
    assert scope[0]["status"] == "INFO"
    assert "decoded all 1 chunks" in scope[0]["detail"]


def test_validate_keeps_keyword_argument_compatibility(good_db: Path) -> None:
    # Builder tests call it as ``validate(db, sample_chunks=..., sample_entries=...)``.
    report = validate(good_db, sample_chunks=1, sample_entries=N_PHASES)
    assert report.fails == 0


# ---------------------------------------------------------------------------
# E5: bogus single-reflection phases
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("d", "signature"),
    [
        (0.7703, "cutoff"),  # 2theta = 180 deg boundary artifact
        (1.0, "placeholder"),  # gemmi placeholder-cell artifact
    ],
)
def test_bogus_lone_peak_fails(good_db: Path, tmp_path: Path, d: float, signature: str) -> None:
    entry_id = _mutate_to_lone_peak(good_db, d)
    report = validate(good_db, 1, N_PHASES)
    statuses = _row(report, "E5.boundary_peak")
    assert any(status == "FAIL" for status, _detail in statuses), report.render()
    detail = next(detail for status, detail in statuses if status == "FAIL")
    assert entry_id in detail
    assert signature in detail
    # The accounting checks must still hold: the mutation is isolated to one phase.
    assert _row(report, "C4.peak_total")[0][0] == "PASS"
    assert _row(report, "C6.peak_off")[0][0] == "PASS"

    out = tmp_path / "bogus.json"
    assert validate_main([str(good_db), "--json", str(out)]) == 1
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "FAILED"
    assert payload["fails"] >= 1


def test_collapsed_double_line_peak_fails(good_db: Path, tmp_path: Path) -> None:
    """A two-line collapsed pattern is the same 180 deg artifact (PCOD 96-254-30089).

    ``n_peaks == 1`` cannot see it, so this pins the strongest-reflection anchor:
    the bogus boundary line at 1000 dominates the pattern while the surviving
    real line is 0.013.
    """
    entry_id = _mutate_to_collapsed_double_peak(good_db)
    report = validate(good_db, 1, N_PHASES)
    statuses = _row(report, "E5.boundary_peak")
    assert any(status == "FAIL" for status, _detail in statuses), report.render()
    detail = next(detail for status, detail in statuses if status == "FAIL")
    assert entry_id in detail
    assert "cutoff" in detail
    # The accounting checks must still hold: only the last phase was rewritten.
    assert _row(report, "C4.peak_total")[0][0] == "PASS"
    assert _row(report, "C6.peak_off")[0][0] == "PASS"

    out = tmp_path / "double.json"
    assert validate_main([str(good_db), "--json", str(out)]) == 1
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "FAILED"


def test_semi_collapsed_boundary_is_only_a_warning(good_db: Path) -> None:
    """Strongest reflection in the cutoff window + populated pattern -> WARN.

    This is the ``96-244-04910`` class: same 180 deg root cause, but the pattern
    still carries measurable intensity elsewhere, so it must not fail acceptance.
    """
    entry_id = _mutate_to_collapsed_double_peak(good_db, real_i=100.0)
    report = validate(good_db, 1, N_PHASES)
    statuses = _row(report, "E5.boundary_peak")
    assert statuses and statuses[0][0] == "WARN"
    detail = statuses[0][1]
    assert entry_id in detail
    assert "still-populated" in detail
    assert report.fails == 0, report.render()
    assert report.warns >= 1


def test_plausible_lone_peak_is_only_a_warning(good_db: Path) -> None:
    entry_id = _mutate_to_lone_peak(good_db, 2.5)
    report = validate(good_db, 1, N_PHASES)
    statuses = _row(report, "E5.boundary_peak")
    assert statuses and statuses[0][0] == "WARN"
    assert entry_id in statuses[0][1]
    assert report.fails == 0, report.render()
    assert report.warns >= 1


# ---------------------------------------------------------------------------
# E6: a phase with no reflection at all
# ---------------------------------------------------------------------------
def test_empty_phase_fails_peak_count(good_db: Path) -> None:
    entry_id = _mutate_to_empty_phase(good_db)
    report = validate(good_db, 1, N_PHASES)
    statuses = _row(report, "E6.peak_count")
    assert any(status == "FAIL" for status, _detail in statuses), report.render()
    detail = next(detail for status, detail in statuses if status == "FAIL")
    assert entry_id in detail
    # The reflection payload was shrunk with the count, so the accounting holds.
    assert _row(report, "C4.peak_total")[0][0] == "PASS"


# ---------------------------------------------------------------------------
# file-level failures
# ---------------------------------------------------------------------------
def test_missing_file_fails(tmp_path: Path) -> None:
    report = validate(tmp_path / "no-such.xinDB", 1, 1)
    assert report.fails == 1
    assert _row(report, "F0.exists")[0][0] == "FAIL"
    assert validate_main([str(tmp_path / "no-such.xinDB")]) == 1
