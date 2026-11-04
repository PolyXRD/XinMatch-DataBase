"""Tests for the Materials Project external library (mp_scan + MP build/accept)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

from tests.conftest import render_cif
from tools.accept_xindb import accept
from tools.build_cod_xindb import build
from tools.validate_xindb import validate
from xinmatch.core.db.mp_scan import (
    crystal_system_from_symbol,
    entry_id_for_mp,
    parse_mp_header,
    scan_mp_archive,
)

DATA_DATE = "20261031"


def _mp_text(mp_id: int, *, space_group: str = "Fm-3m", cell: float = 5.431) -> str:
    """A minimal MP-style CIF: comment header + a valid P1 structure."""
    header = (
        "# Materials Project summary 2025-09-25\n"
        f"# mp-id: mp-{mp_id} | spacegroup: {space_group} | "
        "energy_above_hull: 0.000000 eV/atom\n"
    )
    return header + render_cif(str(mp_id), a=cell, b=cell, c=cell)


def test_scan_mp_archive_orders_ids_and_reports_missing_header(tmp_path: Path) -> None:
    root = tmp_path / "cif"
    root.mkdir()
    for mp_id in (300, 100, 200):
        (root / f"X{mp_id}.cif").write_text(_mp_text(mp_id), encoding="utf-8")
    (root / "no_header.cif").write_text(render_cif("5"), encoding="utf-8")

    archive = scan_mp_archive(root)

    assert [mp_id for mp_id, _path in archive.entries] == [100, 200, 300]
    assert [p.name for p in archive.no_id] == ["no_header.cif"]
    assert archive.total == 4

    header = parse_mp_header(root / "X100.cif")
    assert header is not None
    assert header.mp_id == 100
    assert header.space_group == "Fm-3m"
    assert header.summary_date == "2025-09-25"


def test_entry_id_and_crystal_system() -> None:
    assert entry_id_for_mp(862690) == "mp-862690"
    assert crystal_system_from_symbol("P6_3/mmc") == "hexagonal"
    assert crystal_system_from_symbol("Fm-3m") == "cubic"
    assert crystal_system_from_symbol(None) is None
    with pytest.raises(ValueError):
        entry_id_for_mp(-1)


def test_build_and_accept_mp_library(tmp_path: Path) -> None:
    root = tmp_path / "cif"
    root.mkdir()
    for mp_id in (300, 100, 200):
        (root / f"X{mp_id}.cif").write_text(_mp_text(mp_id), encoding="utf-8")

    out = tmp_path / f"MP-{DATA_DATE}.xinDB"
    stats = build(
        root,
        out,
        library_kind="MP",
        build_tool="tests.test_mp_build",
        data_date=DATA_DATE,
        tier=2,
        workers=1,
        id_mode="mp",
    )
    assert stats.entries == 3
    assert stats.entries_total == 3

    con = sqlite3.connect(f"file:{out.as_posix()}?mode=ro", uri=True)
    try:
        meta = dict(con.execute("SELECT key, value FROM meta"))
        assert meta["library_kind"] == "MP"
        assert meta["data_date"] == DATA_DATE
        assert meta["tier_default"] == "2"
        ids = [r[0] for r in con.execute("SELECT entry_id FROM entries ORDER BY rowid")]
        assert ids == ["mp-100", "mp-200", "mp-300"]
        row = con.execute(
            "SELECT space_group, crystal_system, tier, journal, author, quality "
            "FROM entries WHERE entry_id = 'mp-100'"
        ).fetchone()
        assert row == ("Fm-3m", "cubic", 2, "Materials Project", "Materials Project", "DFT")
    finally:
        con.close()

    report = validate(out, -1, 3)
    assert report.fails == 0, report.render()

    accepted = accept(out, kind="auto", mp_root=root, sample_chunks=None)
    assert accepted.kind == "MP"
    assert accepted.fails == 0, accepted.render()
    assert accepted.stats["entries"] == 3
    assert accepted.stats["id_raw_min"] == 100
    assert accepted.stats["id_raw_max"] == 300
