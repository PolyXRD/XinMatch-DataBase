"""Tests for the new ``tools/accept_xindb.py`` six-library acceptance tool.

What is covered
---------------
* a real 12-phase ``COD-Inorg`` build passes: metadata, counts, full id-band
  scan, sampled archive ownership, element-index equality and membership, peak
  payload, and the exhaustive bad-phase scan;
* the tool is genuinely source-agnostic: the same build slice is accepted as
  ``COD``, ``PCOD`` (``+20M`` offset) and ``TCOD`` (``+100M``), and an offset
  that contradicts the declared source fails the id-band check;
* the JSON contract (``verdict/fails/warns/checks/stats``) and the exit code
  (1 on FAIL) are stable;
* each bad-phase signature and each form of corruption is caught: lone
  boundary/placeholder reflections, an empty phase, a dropped element row, a
  truncated element list, a zeroed peak payload, a missing source archive and
  a broken ``meta.cif_source``;
* ``--manifest`` supplies ownership evidence when the CIF archive is absent.

Every library used here is built under pytest's ``tmp_path``; the published
libraries in ``data/`` and the real CIF archives are never touched.
"""

from __future__ import annotations

import json
import shutil
import sqlite3
import struct
from pathlib import Path

import pytest

from tools.accept_xindb import accept
from tools.accept_xindb import main as accept_main

#: Snapshot date used by every build, so file names stay stable.
DATA_DATE = "20261020"

#: Number of phases in the ``tiny_archive`` fixture.
N_PHASES = 12

#: Archive id of the first phase in the ``tiny_archive`` fixture.
FIRST_ID = 1_000_000


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


def _status(report, code: str) -> str | None:
    rows = _row(report, code)
    return rows[0][0] if rows else None


def _set_meta(db: Path, key: str, value: str) -> None:
    conn = sqlite3.connect(str(db))
    try:
        conn.execute("UPDATE meta SET value = ? WHERE key = ?", (value, key))
        conn.commit()
    finally:
        conn.close()


def _mutate_to_lone_peak(db: Path, d: float) -> str:
    """Rewrite the last phase of *db* into a single-reflection phase."""
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
        new_top_d[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 4] = struct.pack("<f", d)
        new_top_i[row_in_chunk * row_bytes : row_in_chunk * row_bytes + 4] = struct.pack(
            "<f", 1000.0
        )
        conn.execute("UPDATE entries SET n_peaks = 1, peak_n = 1 WHERE rowid = ?", (rowid,))
        conn.execute(
            "UPDATE peak_chunks SET top_d = ?, top_i = ?, all_d = ?, all_i = ? "
            "WHERE chunk_id = ?",
            (
                bytes(new_top_d),
                bytes(new_top_i),
                bytes(all_d)[: offset * 4] + struct.pack("<f", d),
                bytes(all_i)[: offset * 4] + struct.pack("<f", 1000.0),
                chunk_id,
            ),
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
        conn.execute("UPDATE entries SET n_peaks = 0, peak_n = 0 WHERE rowid = ?", (rowid,))
        conn.execute(
            "UPDATE peak_chunks SET top_d = ?, top_i = ?, all_d = ?, all_i = ? "
            "WHERE chunk_id = ?",
            (
                bytes(new_top_d),
                bytes(new_top_i),
                bytes(all_d)[: offset * 4],
                bytes(all_i)[: offset * 4],
                chunk_id,
            ),
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
    """Rewrite the last phase of *db* into a two-line collapsed 180 deg pattern.

    This is the PCOD ``96-254-30089`` signature: the bogus boundary line carries
    the whole normalised intensity (1000) while the one surviving real
    reflection is tiny.  The B1 check must be anchored on the strongest
    reflection, not on ``n_peaks == 1``, to catch it.
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


def _copy_library(source: Path, target: Path) -> Path:
    shutil.copy2(source, target)
    return target


@pytest.fixture
def good_db(tmp_path: Path, tiny_archive: Path) -> Path:
    """A complete 12-phase ``COD-Inorg`` slice built from ``tiny_archive``."""
    return _build_small(tiny_archive, tmp_path / f"COD-Inorg-{DATA_DATE}.xinDB")


# ---------------------------------------------------------------------------
# happy path and CLI contract
# ---------------------------------------------------------------------------
def test_good_library_passes(good_db: Path) -> None:
    report = accept(good_db)
    assert report.fails == 0, report.render()
    assert report.verdict == "PASS"
    assert report.kind == "COD"
    assert report.dataset == "inorg"
    for code in (
        "M1.schema",
        "C1.entry_count",
        "R1.id_band",
        "A1.ownership",
        "E1.element_sets",
        "E5.pairs",
        "E4.membership",
        "P1.payload",
        "B1.single_peak",
        "B2.peak_count",
    ):
        assert _status(report, code) == "PASS", f"{code} -> {_row(report, code)}"
    stats = report.stats
    assert stats["entries"] == N_PHASES
    assert stats["id_raw_min"] == FIRST_ID
    assert stats["id_raw_max"] == FIRST_ID + N_PHASES - 1
    assert stats["bad_phase_ids"] == []


def test_cli_json_contract(good_db: Path, tmp_path: Path, capsys) -> None:
    out = tmp_path / "accept.json"
    code = accept_main(
        [
            str(good_db),
            "--sample-ids",
            str(N_PHASES),
            "--sample-chunks",
            "1",
            "--json",
            str(out),
        ]
    )
    assert code == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert set(payload) == {"db", "kind", "dataset", "verdict", "fails", "warns", "checks", "stats"}
    assert payload["verdict"] == "PASS"
    assert payload["kind"] == "COD"
    assert payload["fails"] == 0
    assert payload["stats"]["entries"] == N_PHASES
    assert all(set(check) == {"code", "status", "detail"} for check in payload["checks"])
    assert "PASS" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# six-library generality: source / id-offset derivation
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("kind", "offset"),
    [("COD", 0), ("PCOD", 20_000_000), ("TCOD", 100_000_000)],
)
def test_accepts_each_source_with_its_offset(
    tmp_path: Path, tiny_archive: Path, kind: str, offset: int
) -> None:
    db = _build_small(
        tiny_archive,
        tmp_path / f"{kind}-Inorg-{DATA_DATE}.xinDB",
        library_kind=f"{kind}-Inorg",
        id_offset=offset,
    )
    report = accept(db)
    assert report.fails == 0, report.render()
    assert report.kind == kind
    assert report.stats["id_offset"] == offset
    assert _status(report, "R1.id_band") == "PASS"
    assert _status(report, "A1.ownership") == "PASS"
    assert report.stats["id_raw_min"] == FIRST_ID
    assert report.stats["id_raw_max"] == FIRST_ID + N_PHASES - 1


def test_accepts_full_dataset_library(tmp_path: Path, tiny_archive: Path) -> None:
    # Same builder, same archive slice, but declared as the Full product: the
    # dataset must be derived as "full" and the ownership check still resolve.
    db = _build_small(
        tiny_archive,
        tmp_path / f"COD-Full-{DATA_DATE}.xinDB",
        library_kind="COD-Full",
    )
    report = accept(db)
    assert report.fails == 0, report.render()
    assert report.dataset == "full"
    assert _status(report, "A1.ownership") == "PASS"


def test_wrong_offset_fails_id_band(tmp_path: Path, tiny_archive: Path) -> None:
    # A COD library whose ids were shifted into the PCOD band is not owned by
    # the COD archive and must not be accepted.
    db = _build_small(
        tiny_archive,
        tmp_path / f"COD-Inorg-{DATA_DATE}.xinDB",
        library_kind="COD-Inorg",
        id_offset=20_000_000,
    )
    report = accept(db, sample_chunks=1)
    assert report.fails >= 1
    assert _status(report, "R1.id_band") == "FAIL"
    assert accept_main([str(db), "--sample-chunks", "1"]) == 1


# ---------------------------------------------------------------------------
# bad phases and corruption are caught
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("d", "counter"),
    [(0.7703, "boundary_phases"), (1.0, "placeholder_phases")],
)
def test_bogus_lone_peak_fails(
    good_db: Path, tmp_path: Path, d: float, counter: str
) -> None:
    entry_id = _mutate_to_lone_peak(good_db, d)
    out = tmp_path / "bad.json"
    assert accept_main([str(good_db), "--json", str(out)]) == 1
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "FAIL"
    assert payload["stats"][counter] == 1
    assert entry_id in [item["entry_id"] for item in payload["stats"]["bad_phase_ids"]]
    checks = {check["code"]: check for check in payload["checks"]}
    assert checks["B1.single_peak"]["status"] == "FAIL"
    assert entry_id in checks["B1.single_peak"]["detail"]


def test_collapsed_double_line_peak_fails(good_db: Path, tmp_path: Path) -> None:
    """The PCOD 96-254-30089 signature: no n_peaks==1, still a bogus phase."""
    entry_id = _mutate_to_collapsed_double_peak(good_db)
    out = tmp_path / "double.json"
    assert accept_main([str(good_db), "--json", str(out)]) == 1
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "FAIL"
    assert payload["stats"]["boundary_phases"] == 1
    assert payload["stats"]["placeholder_phases"] == 0
    reported = {item["entry_id"]: item for item in payload["stats"]["bad_phase_ids"]}
    assert entry_id in reported
    assert reported[entry_id]["n_peaks"] == 2
    checks = {check["code"]: check for check in payload["checks"]}
    assert checks["B1.single_peak"]["status"] == "FAIL"
    assert "cutoff" in checks["B1.single_peak"]["detail"]


def test_semi_collapsed_boundary_is_warn_only(good_db: Path, tmp_path: Path) -> None:
    """The 96-244-04910 class is WARN, never FAIL (verdict PASS WITH WARNINGS)."""
    entry_id = _mutate_to_collapsed_double_peak(good_db, real_i=100.0)
    out = tmp_path / "semi.json"
    assert accept_main([str(good_db), "--json", str(out)]) == 0
    payload = json.loads(out.read_text(encoding="utf-8"))
    assert payload["verdict"] == "PASS WITH WARNINGS"
    assert payload["fails"] == 0
    stats = payload["stats"]
    assert stats["boundary_phases"] == 0
    assert stats["boundary_warn_phases"] == 1
    assert stats["bad_phase_ids"] == []
    item = stats["boundary_warn_ids"][0]
    assert item["entry_id"] == entry_id
    assert item["n_peaks"] == 2
    assert item["rest_sum"] >= 5.0
    checks = {check["code"]: check for check in payload["checks"]}
    assert checks["B1.boundary_warn"]["status"] == "WARN"
    assert entry_id in checks["B1.boundary_warn"]["detail"]


def test_empty_phase_fails_peak_count(good_db: Path) -> None:
    entry_id = _mutate_to_empty_phase(good_db)
    report = accept(good_db, sample_chunks=1)
    assert report.fails == 1, report.render()
    detail = _row(report, "B2.peak_count")[0][1]
    assert entry_id in detail


def test_missing_element_row_fails_set_check(good_db: Path) -> None:
    conn = sqlite3.connect(str(good_db))
    try:
        conn.execute("DELETE FROM element_inv WHERE element = 'Si'")
        conn.commit()
    finally:
        conn.close()
    report = accept(good_db)
    assert _status(report, "E1.element_sets") == "FAIL"
    assert report.fails >= 1
    assert accept_main([str(good_db)]) == 1


def test_truncated_element_list_fails_pair_count(good_db: Path) -> None:
    conn = sqlite3.connect(str(good_db))
    try:
        conn.execute(
            "UPDATE element_inv SET eids = substr(eids, 1, length(eids) - 4) WHERE element = 'O'"
        )
        conn.commit()
    finally:
        conn.close()
    report = accept(good_db)
    # The symbol set is unchanged, so the full pair count is what catches it.
    assert _status(report, "E1.element_sets") == "PASS"
    assert _status(report, "E5.pairs") == "FAIL"


def test_zeroed_peak_payload_fails(good_db: Path) -> None:
    conn = sqlite3.connect(str(good_db))
    try:
        conn.execute(
            "UPDATE peak_chunks SET all_d = zeroblob(length(all_d)) WHERE chunk_id = 0"
        )
        conn.commit()
    finally:
        conn.close()
    report = accept(good_db, sample_chunks=1)
    assert _status(report, "P1.payload") == "FAIL"
    assert report.fails >= 1


# ---------------------------------------------------------------------------
# archive ownership
# ---------------------------------------------------------------------------
def test_missing_source_files_fail_ownership(good_db: Path, tmp_path: Path) -> None:
    empty_archive = tmp_path / "empty_archive"
    empty_archive.mkdir()
    db = _copy_library(good_db, tmp_path / "moved.xinDB")
    _set_meta(db, "cif_source", str(empty_archive))
    report = accept(db, sample_chunks=1)
    assert _status(report, "A1.ownership") == "FAIL"
    assert report.fails == 1
    assert accept_main([str(db), "--sample-chunks", "1"]) == 1


def _write_manifest(path: Path, *, stored_ids: list[int], lib: str, inorg_all: bool) -> None:
    conn = sqlite3.connect(str(path))
    try:
        conn.executescript(
            """
            CREATE TABLE files (
                cod_id INTEGER PRIMARY KEY,
                lib TEXT NOT NULL,
                rel TEXT NOT NULL,
                formula TEXT,
                formula_red TEXT,
                elements TEXT,
                is_inorg INTEGER NOT NULL,
                size INTEGER,
                mtime INTEGER
            );
            """
        )
        rows = []
        for index, stored_id in enumerate(stored_ids):
            archive_id = stored_id if lib == "COD" else stored_id - (
                20_000_000 if lib == "PCOD" else 100_000_000
            )
            rel = f"1/00/00/{archive_id}.cif"
            is_inorg = 1 if inorg_all or index else 0
            rows.append((stored_id, lib, rel, None, None, None, is_inorg, None, None))
        conn.executemany(
            "INSERT INTO files (cod_id, lib, rel, formula, formula_red, elements, "
            "is_inorg, size, mtime) VALUES (?,?,?,?,?,?,?,?,?)",
            rows,
        )
        conn.commit()
    finally:
        conn.close()


def test_manifest_supplies_ownership_evidence(good_db: Path, tmp_path: Path) -> None:
    db = _copy_library(good_db, tmp_path / "moved.xinDB")
    _set_meta(db, "cif_source", str(tmp_path / "gone_archive"))
    manifest = tmp_path / "manifest.sqlite"
    _write_manifest(
        manifest,
        stored_ids=[FIRST_ID + index for index in range(N_PHASES)],
        lib="COD",
        inorg_all=True,
    )
    report = accept(
        db,
        manifest=manifest,
        data_root=tmp_path / "no-root",
        sample_ids=N_PHASES,
        sample_chunks=1,
    )
    assert report.fails == 0, report.render()
    assert _status(report, "A1.ownership") == "PASS"
    assert report.stats["manifest_hits"] == N_PHASES


def test_manifest_inorganic_flag_is_cross_checked(good_db: Path, tmp_path: Path) -> None:
    db = _copy_library(good_db, tmp_path / "moved.xinDB")
    _set_meta(db, "cif_source", str(tmp_path / "gone_archive"))
    manifest = tmp_path / "manifest.sqlite"
    _write_manifest(
        manifest,
        stored_ids=[FIRST_ID + index for index in range(N_PHASES)],
        lib="COD",
        inorg_all=False,
    )
    report = accept(
        db,
        manifest=manifest,
        data_root=tmp_path / "no-root",
        sample_ids=N_PHASES,
        sample_chunks=1,
    )
    # Ownership is still evidenced; the flag disagreement is a warning only.
    assert report.fails == 0, report.render()
    assert _status(report, "A1.ownership") == "PASS"
    assert _status(report, "A2.manifest_flags") == "WARN"
    assert report.warns >= 1


def test_missing_file_fails(tmp_path: Path) -> None:
    report = accept(tmp_path / "no-such.xinDB")
    assert report.fails == 1
    assert _status(report, "F0.exists") == "FAIL"
    assert accept_main([str(tmp_path / "no-such.xinDB")]) == 1


def test_source_cif_candidates_cover_all_three_archive_layouts() -> None:
    """The ownership check must know COD / TCOD / PCOD fan-out shapes.

    Real archives (verified on ``sync-COD``):
    COD ``1/00/00/1000005.cif``, TCOD ``10/00/00/10000001.cif``,
    PCOD ``1/100/1000000.cif`` / ``3/310/3100000.cif``.
    """
    from tools.accept_xindb import _source_cif_candidates

    assert "1/00/00/1000005.cif" in _source_cif_candidates(1_000_005)
    assert "10/00/00/10000001.cif" in _source_cif_candidates(10_000_001)
    assert "1/100/1000000.cif" in _source_cif_candidates(1_000_000)
    assert "3/310/3100000.cif" in _source_cif_candidates(3_100_000)
    assert _source_cif_candidates(123) == []
