"""Build-entry tests for the six-library CLI (task-6).

What is covered
---------------
* the central :data:`~xinmatch.core.db.cif_scan.SOURCES` registry: the three
  archive families, both roots each and the 3COD id shifts;
* ``--source``/``--set`` -> ``library_kind`` / archive root / ``id_offset``
  derivation for all six combinations, the legacy COD-only ``--kind`` alias
  and its conflict errors;
* default output directory semantics (``<repo>/data`` unless
  ``XINMATCH_XINDB_DIR`` explicitly overrides it) and ``XRD_DATA_ROOT``
  isolation;
* real end-to-end builds of 12-phase slices from the ``tiny_archive`` fixture,
  including an equality check against the old two-call (``parse_meta`` +
  ``compute_pattern``) loop and a few ``validate_xindb`` checks.

Every build writes under pytest's ``tmp_path``.  The published libraries in
``data/`` and the real CIF archives (``sync-COD``) are never touched.
"""

from __future__ import annotations

import os
import sqlite3
import subprocess
import sys
from pathlib import Path

import pytest

from tools.build_cod_inorg_xindb import main as inorg_main
from tools.build_cod_xindb import (
    PROJECT_ROOT,
    default_out_dir,
    resolve_plan,
)
from tools.build_cod_xindb import (
    main as build_main,
)
from xinmatch.core.db.cif_scan import (
    DEFAULT_COD_ROOT,
    DEFAULT_INORG_COD_ROOT,
    DEFAULT_INORG_PCOD_ROOT,
    DEFAULT_INORG_TCOD_ROOT,
    DEFAULT_PCOD_ROOT,
    DEFAULT_TCOD_ROOT,
    FULL_ROOTS,
    INORG_ROOTS,
    SOURCES,
    SUBSETS,
    iter_cod_cifs,
)

#: Snapshot date used by every build in this module, so file names are stable.
DATA_DATE = "20261020"

#: Archive id of the first phase in the ``tiny_archive`` fixture.
FIRST_ID = 1_000_000

#: Number of phases in the ``tiny_archive`` fixture.
N_PHASES = 12


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _read_meta(db: Path) -> dict[str, str]:
    """Return the ``meta`` table of *db* as a dict."""
    conn = sqlite3.connect(f"file:{db.as_posix()}?mode=ro", uri=True)
    try:
        return {str(k): str(v) for k, v in conn.execute("SELECT key, value FROM meta")}
    finally:
        conn.close()


def _read_rows(db: Path, sql: str) -> list[tuple]:
    """Run one read-only query and return all rows."""
    conn = sqlite3.connect(f"file:{db.as_posix()}?mode=ro", uri=True)
    try:
        return [tuple(row) for row in conn.execute(sql)]
    finally:
        conn.close()


def _numeric_ids(db: Path) -> list[int]:
    """Numeric ids recovered from the stored ``96-XXX-XXXX`` entry ids."""
    rows = _read_rows(db, "SELECT entry_id FROM entries ORDER BY rowid")
    return [int("".join(str(row[0]).split("-")[1:])) for row in rows]


def _build_slice(
    out: Path,
    cif_root: Path,
    *args: str,
    wrapper=None,
) -> int:
    """Run a 12-phase build of *cif_root* into *out* and return the exit code."""
    entry = wrapper or build_main
    return entry(
        [
            "--cif-root", str(cif_root),
            "--limit", str(N_PHASES),
            "--workers", "1",
            "--data-date", DATA_DATE,
            "--out", str(out),
            *args,
        ]
    )


def _assert_slice(db: Path, *, library_kind: str, id_offset: int) -> dict[str, str]:
    """Assert the *db* is a complete 12-phase slice with the expected identity."""
    assert db.is_file(), f"build did not produce {db}"
    meta = _read_meta(db)
    assert meta["library_kind"] == library_kind
    assert meta["data_date"] == DATA_DATE
    assert meta["entry_count"] == str(N_PHASES)

    ids = _numeric_ids(db)
    assert ids == [FIRST_ID + i + id_offset for i in range(N_PHASES)]

    rows = _read_rows(db, "SELECT n_peaks FROM entries")
    assert len(rows) == N_PHASES
    assert all(int(row[0]) > 0 for row in rows), "every tiny_archive CIF has peaks"

    chunks = _read_rows(db, "SELECT chunk_id, n FROM peak_chunks ORDER BY chunk_id")
    assert [int(row[0]) for row in chunks] == list(range(len(chunks)))
    assert sum(int(row[1]) for row in chunks) == N_PHASES
    return meta


def _build_legacy_two_call(db: Path, cif_root: Path, *, library_kind: str) -> None:
    """Build *db* through the **old two-call loop**, as a parity oracle.

    Kept inside the test suite on purpose: the builder itself must never have a
    second, drifting code path, but the old ``parse_meta`` + ``compute_pattern``
    sequence stays the reference the single-parse pipeline has to match.
    """
    from xinmatch.core.db.cif_meta import parse_meta
    from xinmatch.core.db.cif_scan import entry_id_for_cod
    from xinmatch.core.db.peak_compute import (
        DEFAULT_DMIN,
        DEFAULT_MAX_PEAKS,
        compute_pattern,
    )
    from xinmatch.core.db.schema_v3 import CHUNK_SIZE, TOP_K
    from xinmatch.core.db.writer_v3 import XinDBWriter

    writer = XinDBWriter(
        db,
        top_k=TOP_K,
        chunk_size=CHUNK_SIZE,
        tier=1,
        cif_source=str(cif_root),
        library_kind=library_kind,
        build_tool="tools.build_cod_xindb",
        data_date=DATA_DATE,
        build_time="2026-10-20T00:00:00+08:00",
    )
    try:
        for cod_id, cif_path in iter_cod_cifs(cif_root):
            meta = parse_meta(cif_path, entry_id=entry_id_for_cod(cod_id))
            if meta is None:
                continue
            peaks = compute_pattern(
                cif_path,
                dmin=DEFAULT_DMIN,
                max_peaks=DEFAULT_MAX_PEAKS,
                normalize=True,
            )
            if peaks is None:
                continue
            writer.add(meta, peaks[0], peaks[1])
    finally:
        writer.close()


# ---------------------------------------------------------------------------
# source registry
# ---------------------------------------------------------------------------
def test_source_registry_reuses_the_default_constants() -> None:
    """``SOURCES`` must expose the existing roots, not second copies of them."""
    assert set(SOURCES) == {"COD", "PCOD", "TCOD"}
    assert SUBSETS == ("inorg", "full")

    assert SOURCES["COD"].full_root == DEFAULT_COD_ROOT
    assert SOURCES["COD"].inorg_root == DEFAULT_INORG_COD_ROOT
    assert SOURCES["PCOD"].full_root == DEFAULT_PCOD_ROOT
    assert SOURCES["PCOD"].inorg_root == DEFAULT_INORG_PCOD_ROOT
    assert SOURCES["TCOD"].full_root == DEFAULT_TCOD_ROOT
    assert SOURCES["TCOD"].inorg_root == DEFAULT_INORG_TCOD_ROOT

    assert [SOURCES[name].id_offset for name in ("COD", "PCOD", "TCOD")] == [
        0,
        20_000_000,
        100_000_000,
    ]
    assert {name: spec.full_root for name, spec in SOURCES.items()} == FULL_ROOTS
    assert {name: spec.inorg_root for name, spec in SOURCES.items()} == INORG_ROOTS


@pytest.mark.parametrize(
    ("source", "subset", "expected_kind", "offset", "root"),
    [
        ("COD", "inorg", "COD-Inorg", 0, DEFAULT_INORG_COD_ROOT),
        ("COD", "full", "COD-Full", 0, DEFAULT_COD_ROOT),
        ("PCOD", "inorg", "PCOD-Inorg", 20_000_000, DEFAULT_INORG_PCOD_ROOT),
        ("PCOD", "full", "PCOD-Full", 20_000_000, DEFAULT_PCOD_ROOT),
        ("TCOD", "inorg", "TCOD-Inorg", 100_000_000, DEFAULT_INORG_TCOD_ROOT),
        ("TCOD", "full", "TCOD-Full", 100_000_000, DEFAULT_TCOD_ROOT),
    ],
)
def test_six_combinations_derive_kind_root_and_offset(
    source: str, subset: str, expected_kind: str, offset: int, root: Path
) -> None:
    plan = resolve_plan(
        ["--source", source, "--set", subset, "--data-date", DATA_DATE]
    )
    assert plan.source == source
    assert plan.subset == subset
    assert plan.library_kind == expected_kind
    assert plan.id_offset == offset
    assert plan.cif_root == root
    assert plan.out.name == f"{expected_kind}-{DATA_DATE}.xinDB"
    assert plan.warnings == ()


def test_default_out_dir_is_repo_data_and_env_is_an_override(
    monkeypatch, tmp_path: Path
) -> None:
    monkeypatch.delenv("XINMATCH_XINDB_DIR", raising=False)
    assert default_out_dir() == PROJECT_ROOT / "data"

    override = tmp_path / "shared-xindb"
    monkeypatch.setenv("XINMATCH_XINDB_DIR", str(override))
    assert default_out_dir() == override


def test_legacy_kind_alias_maps_to_cod() -> None:
    full = resolve_plan(["--kind", "full", "--data-date", DATA_DATE])
    inorg = resolve_plan(["--kind", "inorg", "--data-date", DATA_DATE])
    assert (full.library_kind, full.id_offset, full.subset) == ("COD-Full", 0, "full")
    assert (inorg.library_kind, inorg.id_offset, inorg.subset) == (
        "COD-Inorg",
        0,
        "inorg",
    )


def test_kind_and_set_together_are_rejected() -> None:
    with pytest.raises(SystemExit) as excinfo:
        resolve_plan(["--kind", "full", "--set", "inorg"])
    assert excinfo.value.code == 2


def test_kind_is_rejected_for_non_cod_sources() -> None:
    with pytest.raises(SystemExit) as excinfo:
        resolve_plan(["--source", "PCOD", "--kind", "inorg"])
    assert excinfo.value.code == 2


def test_explicit_overrides_are_flagged(tmp_path: Path) -> None:
    plan = resolve_plan(
        [
            "--source", "PCOD", "--set", "full",
            "--cif-root", str(tmp_path),
            "--id-offset", "7",
            "--data-date", DATA_DATE,
        ]
    )
    # The overrides are applied ...
    assert plan.cif_root == tmp_path
    assert plan.id_offset == 7
    # ... but never silently: the library_kind still says PCOD-Full while the
    # root/offset no longer match the convention, so both are warned about.
    assert plan.library_kind == "PCOD-Full"
    assert len(plan.warnings) == 2
    assert "--cif-root" in plan.warnings[0]
    assert "--id-offset" in plan.warnings[1]


def test_xrd_data_root_env_is_honoured_on_import(tmp_path: Path) -> None:
    """``XRD_DATA_ROOT`` is read at import time; verify it in a fresh process."""
    data_root = tmp_path / "sync-COD"
    env = dict(os.environ, XRD_DATA_ROOT=str(data_root))
    code = (
        "import xinmatch.core.db.cif_scan as m;"
        "print(m.DEFAULT_COD_ROOT);"
        "print(m.SOURCES['PCOD'].full_root);"
        "print(m.INORG_ROOTS['TCOD'])"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        cwd=str(PROJECT_ROOT),
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.splitlines() == [
        str(data_root / "COD"),
        str(data_root / "PCOD"),
        str(data_root / "inorg_cif" / "TCOD"),
    ]


# ---------------------------------------------------------------------------
# end-to-end builds (12-phase slices of tiny_archive)
# ---------------------------------------------------------------------------
def test_build_cod_inorg_default_output_path(
    tmp_path: Path, tiny_archive: Path, monkeypatch, capsys
) -> None:
    """No ``--out``: the name must be derived into the env-overridden dir."""
    out_dir = tmp_path / "shared-xindb"
    monkeypatch.setenv("XINMATCH_XINDB_DIR", str(out_dir))
    code = build_main(
        [
            "--source", "COD", "--set", "inorg",
            "--cif-root", str(tiny_archive),
            "--limit", str(N_PHASES), "--workers", "1",
            "--data-date", DATA_DATE,
        ]
    )
    assert code == 0
    captured = capsys.readouterr()
    assert "[WARN] --cif-root" in captured.out

    db = out_dir / f"COD-Inorg-{DATA_DATE}.xinDB"
    meta = _assert_slice(db, library_kind="COD-Inorg", id_offset=0)
    assert meta["build_tool"] == "tools.build_cod_xindb"

    # The published libraries are not the only consumer of the format: the
    # validator shipped next to the builder must accept the slice as well.
    from tools.validate_xindb import validate

    report = validate(db, sample_chunks=1, sample_entries=N_PHASES)
    assert report.fails == 0, report.render()


def test_build_cod_full_slice(tmp_path: Path, tiny_archive: Path) -> None:
    db = tmp_path / "accept" / f"COD-Full-{DATA_DATE}.xinDB"
    code = _build_slice(db, tiny_archive, "--source", "COD", "--set", "full")
    assert code == 0
    meta = _assert_slice(db, library_kind="COD-Full", id_offset=0)
    assert meta["build_tool"] == "tools.build_cod_xindb"


def test_build_cod_inorg_legacy_kind_alias(tmp_path: Path, tiny_archive: Path) -> None:
    db = tmp_path / "accept" / f"COD-Full-{DATA_DATE}.xinDB"
    code = _build_slice(db, tiny_archive, "--kind", "full")
    assert code == 0
    _assert_slice(db, library_kind="COD-Full", id_offset=0)


def test_build_cod_inorg_wrapper_defaults_to_inorg(
    tmp_path: Path, tiny_archive: Path
) -> None:
    db = tmp_path / "accept" / f"COD-Inorg-{DATA_DATE}.xinDB"
    code = _build_slice(
        db, tiny_archive, "--source", "COD", wrapper=inorg_main
    )
    assert code == 0
    meta = _assert_slice(db, library_kind="COD-Inorg", id_offset=0)
    # Provenance still says which entry point built the published inorg libs.
    assert meta["build_tool"] == "tools.build_cod_inorg_xindb"


def test_build_pcod_inorg_slice_applies_id_shift(
    tmp_path: Path, tiny_archive: Path, capsys
) -> None:
    db = tmp_path / "accept" / f"PCOD-Inorg-{DATA_DATE}.xinDB"
    code = _build_slice(
        db, tiny_archive, "--source", "PCOD", "--set", "inorg"
    )
    assert code == 0
    captured = capsys.readouterr()
    assert "[WARN] --cif-root" in captured.out

    meta = _assert_slice(db, library_kind="PCOD-Inorg", id_offset=20_000_000)
    assert meta["library_kind"] == "PCOD-Inorg"
    # The display-id shape is preserved, only the numeric range moved.
    entry_ids = _read_rows(db, "SELECT entry_id FROM entries ORDER BY rowid")
    assert entry_ids[0][0] == "96-210-00000"


def test_start_and_limit_select_a_window(tmp_path: Path, tiny_archive: Path) -> None:
    db = tmp_path / "accept" / f"COD-Inorg-{DATA_DATE}.xinDB"
    code = build_main(
        [
            "--source", "COD", "--set", "inorg",
            "--cif-root", str(tiny_archive),
            "--start", str(FIRST_ID + 6),
            "--limit", "3",
            "--workers", "1",
            "--data-date", DATA_DATE,
            "--out", str(db),
        ]
    )
    assert code == 0
    assert _numeric_ids(db) == [FIRST_ID + 6, FIRST_ID + 7, FIRST_ID + 8]


@pytest.mark.parametrize("source", ["COD", "PCOD", "TCOD"])
@pytest.mark.parametrize("subset", ["inorg", "full"])
def test_all_six_combinations_build_with_every_flag(
    tmp_path: Path, tiny_archive: Path, source: str, subset: str
) -> None:
    """``--resume``/``--tier``/``--log-every``/-v must work in all six shapes."""
    kind = SOURCES[source].library_kind(subset)
    db = tmp_path / "accept" / f"{kind}-{DATA_DATE}.xinDB"
    code = build_main(
        [
            "--source", source, "--set", subset,
            "--cif-root", str(tiny_archive),
            "--limit", str(N_PHASES), "--workers", "1",
            "--resume", "--tier", "0", "--log-every", "2",
            "--data-date", DATA_DATE, "--out", str(db), "-v",
        ]
    )
    assert code == 0
    # ``--resume`` on a fresh file must behave like a fresh build.
    _assert_slice(db, library_kind=kind, id_offset=SOURCES[source].id_offset)
    tiers = _read_rows(db, "SELECT DISTINCT tier FROM entries")
    assert tiers == [(0,)]


def test_parallel_workers_produce_the_same_library(
    tmp_path: Path, tiny_archive: Path
) -> None:
    """``--workers 2`` must not change the output (archive order is kept)."""
    serial = tmp_path / "accept" / "serial.xinDB"
    parallel = tmp_path / "accept" / "parallel.xinDB"
    for out, workers in ((serial, "1"), (parallel, "2")):
        proc = subprocess.run(
            [
                sys.executable, "-m", "tools.build_cod_xindb",
                "--source", "COD", "--set", "inorg",
                "--cif-root", str(tiny_archive),
                "--limit", str(N_PHASES), "--workers", workers,
                "--data-date", DATA_DATE, "--out", str(out),
            ],
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            check=False,
        )
        assert proc.returncode == 0, proc.stderr

    for sql in (
        "SELECT * FROM entries ORDER BY rowid",
        "SELECT chunk_id, n, top_d, top_i, all_d, all_i "
        "FROM peak_chunks ORDER BY chunk_id",
        "SELECT element, eids FROM element_inv ORDER BY element",
    ):
        assert _read_rows(serial, sql) == _read_rows(parallel, sql)
    # Identical peaks + identical meta means the two files describe one library.
    serial_meta = _read_meta(serial)
    parallel_meta = _read_meta(parallel)
    serial_meta.pop("build_time", None)
    parallel_meta.pop("build_time", None)
    assert serial_meta == parallel_meta


def test_single_parse_build_matches_the_old_two_call_loop(
    tmp_path: Path, tiny_archive: Path
) -> None:
    """Parity: the CLI (``compute_phase``) == old ``parse_meta`` + ``compute_pattern``."""
    from_worker = tmp_path / "accept" / f"COD-Inorg-{DATA_DATE}.xinDB"
    legacy = tmp_path / "accept" / "legacy.xinDB"
    assert _build_slice(
        from_worker, tiny_archive, "--source", "COD", "--set", "inorg"
    ) == 0
    _build_legacy_two_call(legacy, tiny_archive, library_kind="COD-Inorg")

    # Entries: every scalar column, in archive order.
    entries_sql = "SELECT * FROM entries ORDER BY rowid"
    assert _read_rows(from_worker, entries_sql) == _read_rows(legacy, entries_sql)

    # Peaks: the columnar payloads are byte-identical (same engine, same order).
    chunks_sql = (
        "SELECT chunk_id, n, top_d, top_i, all_d, all_i "
        "FROM peak_chunks ORDER BY chunk_id"
    )
    assert _read_rows(from_worker, chunks_sql) == _read_rows(legacy, chunks_sql)

    # Element inverted index.
    elements_sql = "SELECT element, eids FROM element_inv ORDER BY element"
    assert _read_rows(from_worker, elements_sql) == _read_rows(legacy, elements_sql)

    # Meta: everything except the build timestamp the CLI stamps at run time.
    new_meta = _read_meta(from_worker)
    old_meta = _read_meta(legacy)
    new_meta.pop("build_time", None)
    old_meta.pop("build_time", None)
    assert new_meta == old_meta
