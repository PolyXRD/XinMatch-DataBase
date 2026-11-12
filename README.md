# XinMatch-DataBase

**Build code, format specification and release artifacts for the XinMatch reference phase database (`.xinDB` v3 format).**

This repository contains **only code and documentation**; the database itself is distributed as
compressed packages in [Releases](../../releases) (68 MB / 87 MB / 1.3 GB per database, not
committed to git history). The data sources are
[Crystallography Open Database (COD)](http://www.crystallography.net/cod/) and its
predicted/theoretical sub-databases PCOD / TCOD; this repository does not perform structure
refinement or rewrite the original CIFs, and only performs reflection calculation and index packaging.

> The sister repository [PolyXRD-DataBase](https://github.com/PolyXRD/PolyXRD-DataBase) publishes
> the **polyDB** format of the same source data (used by PolyXRD). **The two formats have different
> conventions and are not interchangeable** (see §2).

---

## 0. Build location convention (effective 2026-10-24, required reading)

| Item | Convention |
|---|---|
| **Where the database is built** | **In this repository, `XinMatch-DataBase`** (no longer scattered across other checkouts). |
| **Where the CIFs needed for the build are** | The parent directory's **`sync-COD`**: `D:\Project\XRD\sync-COD` (can be overridden with `XRD_DATA_ROOT`). This repository does not carry the original CIFs. |
| **Default artifact location** | `<this repository>/data/<Kind>-<YYYYMMDD>.xinDB`; `--out` has the highest priority, `XINMATCH_XINDB_DIR` overrides only when explicitly set. |
| **Supported scope** | **6 databases**: `COD` / `PCOD` / `TCOD` × `Inorg` (inorganic subset) / `Full` (full database). |
| **Complete operations manual** | [Build and acceptance handbook (Chinese)](docs/构建与验收手册.md) (build matrix, validation commands, release checklist, troubleshooting). |

**One chain**: `sync-COD (CIF) → ①scan_inorg (only needed for Inorg) → ②build (this repository's data/) → ③validate → ④accept/special review`.

---

## 1. Release artifacts (Release `20261020`)

> This repository now supports **6 databases** (`COD/PCOD/TCOD × Inorg/Full`). The `20261020`
> Release actually publishes **three Inorg databases** (COD/PCOD/TCOD-Inorg in the table below);
> the three `*-Full` databases can already be built with the same set of commands (see §3.3 and the
> build manual), but their size and build time are on the order of hours, so they are not included in
> this Release.

| Database | entry_id range | Phase count | Reflection count | Chunk count | Elements | Uncompressed | zip |
|---|---|---|---|---|---|---|---|
| **COD-Inorg** | `96-100-*` | 78,036 | 9,394,774 | 20 | 96 kinds | 118,710,272 B (113.2 MiB) | **87.2 MiB** (91,484,834 B) |
| **PCOD-Inorg** | `96-210-*` | 1,061,953 | 159,274,265 | 260 | 19 kinds | 1,776,353,280 B (1694.1 MiB) | **1.35 GiB** (1,447,251,707 B) |
| **TCOD-Inorg** | `96-110-*` | 46,846 | 5,874,790 | 12 | 91 kinds | 71,352,320 B (68.0 MiB) | **54.1 MiB** (56,694,762 B) |
| **Total** | — | **1,186,835** | **174,543,829** | 292 | — | 1.83 GiB | ≈1.52 GiB |

- The three databases are **mutually exclusive, install as needed** (XinMatch supports mounting
  multiple databases). The COD experimental database is the workhorse for everyday interactive
  search; PCOD (GRINSP predictions) and TCOD (theoretical) serve as a "fallback when the
  experimental database finds nothing".
- The `20261020` in the name is the **source data snapshot date** (`meta.data_date`), not the build
  date; build times are in `meta.build_time` (COD 2026-10-21 17:25, TCOD 2026-10-22 23:29,
  PCOD 2026-10-22 23:46).
- All distributed versions use `journal_mode=DELETE`, with no `-wal`/`-shm` sidecars, so they can
  be placed directly on read-only media (NAS read-only shares, disc images) without worrying about
  "read also requiring write permission".
- Checksum values are in the `SHA256SUMS.txt` attached to the Release.

```bash
unzip COD-Inorg-20261020.xinDB.zip                    # get COD-Inorg-20261020.xinDB
python -m tools.validate_xindb COD-Inorg-20261020.xinDB   # optional: 25-item self-check
```

> On this machine the three databases are already placed under the repository's `data/`
> (`data/COD-Inorg-20261020.xinDB`, etc., not committed to git), and XinMatch's library registry
> `~/.xinmatch/libraries.json` already points to that directory.

> ⚠️ The file name must be `<library_kind>-<data_date>.xinDB`, otherwise `validate_xindb`'s
> `M3.identity` will report FAIL (it reconciles the file name against `meta`). After renaming, just
> rebuild the database; there is no need to fix meta.

### 1.1 2026-11-01 rebuild (after the peak-engine fix, local artifacts in `data/`)

In 2026-10 two classes of "spurious-peak phases" were fixed: the 2θ=180° boundary (`d=λ/2`)
reflection that floating-point rounding put into the peak table + LP divergence collapsing the whole
pattern; and the `d=1.0` spurious peak computed from a default placeholder cell when gemmi cannot
read the cell. The three Inorg databases were rebuilt with the fixed engine,
**the source snapshot is still `20261020`** (`meta.data_date` unchanged). The pre-rebuild old
databases were once placed in `data/_prev/` for entry-by-entry comparison and
**were deleted on 2026-11-04 after acceptance was complete** (the comparison results are kept in the
acceptance report):

| Database | Phase count | Reflection count | Size (bytes) | SHA-256 (first 16) | Acceptance |
|---|---|---|---|---|---|
| **COD-Inorg** (rebuilt) | **78,035** | **9,395,681** | **118,640,640** (113.1 MiB) | `7C08B34DCACA9BAC…` | validate 0F/1W (30 items); accept 0F/1W; 7 bad phases fixed, 1 phase without a cell removed |
| **TCOD-Inorg** (rebuilt) | **46,846** | **5,874,790** | **71,352,320** (68.0 MiB) | `8F2A25D583B577A7…` | validate 0F/1W (29 items); accept 0F/1W; entry-by-entry `(d,i)` **exactly identical** to the old database |
| **PCOD-Inorg** (rebuilt) | **1,061,953** | **159,275,318** | **1,777,357,376** (1694.1 MiB) | `ABB7EC8C1E68CD6F…` | validate 0F/1W (29 items); accept 0F/1W; 7 bad phases fixed; 17/1,061,953 payloads changed |
| **TCOD-Full** (new database) | **361,701** | **53,099,238** | **623,173,632** (594.3 MiB) | `36838FBE1F18EDBF…` | validate 0F/2W (29 items); accept 0F/1W; 46,846/46,846 phases common with TCOD-Inorg byte-for-byte identical; self-recall 10/10 |
| **COD-Full** (new database) | **532,038** | **77,458,154** | **922,120,192** (879.4 MiB) | `E217F28FA03755A1…` | validate 0F/2W (30 items); accept 0F/2W (20 items); 78,035/78,035 phases common with COD-Inorg byte-for-byte identical; self-recall 10/10 |
| **MP** (external database) | **155,201** | **20,900,193** | **251,633,664** (240.0 MiB) | `26D7A9BED83A57C0…` | validate 0F/1W (29 items); accept 0F/1W (19 items); id set exactly equal to the archived 155,201 `mp-id`s; self-recall Top-1/Top-5 20/20 |

> **PCOD-Full is not built**: the 1,061,962 CIFs in the PCOD archive are **all inorganic**
> (manifest `is_inorg=1` 100%), essentially equivalent to `PCOD-Inorg` (the original archive differs
> by only 2 CIFs not ingested); it is skipped by decision.

- **MP external database**: `data/MP-20261031.xinDB` (Materials Project, 155,201 phases, ids are
  native `mp-<n>`, `tier=2` computed structures); XinMatch mounts it via "Add external database".
  For build/acceptance see [MP database build and acceptance (Chinese)](docs/MP-库构建与验收.md), published release
  [`20261031-mp`](https://github.com/PolyXRD/XinMatch-DataBase/releases/tag/20261031-mp)
  (zip + `SHA256SUMS.txt`).

- For acceptance details, the per-phase bad-phase list and before/after comparison see
  [`docs/ACCEPTANCE-20261101-rebuild.md`](docs/ACCEPTANCE-20261101-rebuild.md).
- **Release published**: [`20261020-rebuild`](https://github.com/PolyXRD/XinMatch-DataBase/releases/tag/20261020-rebuild)
  (2026-11-04) contains zips for 5 databases + `SHA256SUMS.txt`; release notes
  [`docs/release-notes-20261104.md`](docs/release-notes-20261104.md). Packaging command:
  `python -m tools.package_release --out _release\<date>` (manual §5).
- The pre-rebuild old databases in `data/_prev/` (COD 113 MB + TCOD 68 MB + PCOD 1.7 GB)
  **were deleted on 2026-11-04**; the old `20261020` Release (containing the three Inorg zips with
  bad phases) was also deleted, leaving only this rebuilt Release.

---

## 2. What `.xinDB` is

**A single-file SQLite** database that stores the "pre-screening matrix" and the "full reflection table" separately, with the reflection data stored only once:

```
COD-Inorg-20261020.xinDB
├── meta          self-describing key/value (peak convention, build parameters, traceability)
├── entries       one row per phase (chemical formula / space group / cell / density / elements + CSR offsets)
├── peak_chunks   columnar reflection chunks (one chunk per 4096 phases)
└── element_inv   element inverted index (element → 0-based phase index, uint32 array)
```

| Segment | Shape | Type | Purpose |
|---|---|---|---|
| `top_d`/`top_i` | `n × 32` | `<f4` | pre-screening: NumPy overlap with the measured-peak window, **prune only, no sorting** |
| `all_d`/`all_i` | `Σpeak_n` (CSR) | `<f4` | scoring: all reflections of the candidate phase |
| `element_inv.eids` | variable length | `<u4` | ascending, deduplicated, 0-based phase index |

**The iron rule for addressing**: `entries.rowid` is explicitly written starting from 0 = the
0-based phase index; `element_inv.eids` and the `peak_chunks` chunk number are both based on it, and
**must not** be derived backwards from `entry_id`.
For the complete field-level specification see [xinDB v3 format specification (Chinese)](docs/xinDB格式规范-v3.md).

### Peak convention (hard-coded into `meta.peak_convention`)

```
CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000
```

λ=1.5406 Å (Cu Kα1)｜d 0.7703–30.0 Å｜equivalent reflections merged by 2θ (including multiplicity)｜
**includes Debye-Waller** (`B = 8π²U_iso`, per site)｜intensity normalized to 1000｜at most 150 reflections per phase｜pre-screening width 32

| | xinDB (this repository) | polyDB (PolyXRD-DataBase) |
|---|---|---|
| Container | single-file SQLite | directory (8 segment files + mmap) |
| d range | 0.7703–30 Å | 1.0–26 Å (2θ 5–100.76°) |
| Temperature factor | **yes** (Debye-Waller) | no |
| Peaks per phase | 150 | 100 (top K by intensity) |
| Positioning | search + refinement integrated | GUI interactive search |

> The same CIF yields different absolute intensity systems in the two databases. **Mixing them
> produces systematic deviations**; do not use xinDB intensities to validate polyDB, or vice versa.

### Measured peak tables: `peak_origin = observed`

The convention above describes databases **whose peaks are computed from CIFs**
(COD / PCOD / TCOD / MP), whose `meta.peak_origin` defaults to `computed`. **Measured card
databases** such as ICDD PDF-2 are not computed; when ingested they deliberately do not follow this
convention:

| meta | Computed databases | Measured databases (PDF-2) |
|---|---|---|
| `peak_origin` | `computed` (default) | `observed` |
| `peak_convention` | `CuKa1-1.5406-…-norm1000` | `observed-ICDD-PDF2-2004-norm1000` |
| Intensity scale | norm-1000 | **also norm-1000** (otherwise the scoring weights of the two kinds of databases would not be comparable) |
| NULL cell / peak count >150 | build defect | data fact |

The read side does not need a second reader: `reader_v3` only recognizes `schema_version == 3`.
Build tool `tools/build_pdf2_xindb.py`; **ICDD commercial data does not enter this repository or any
Release**, see [PDF2 database build and acceptance (Chinese)](docs/PDF2-库构建与验收.md) for details.

---

## 3. Building from scratch

### 3.1 Environment

| Dependency | Requirement | Notes |
|---|---|---|
| Python | 3.10+ | tested on 3.10.11 |
| numpy | required | reflection calculation and columnar packaging |
| gemmi | **required** | CIF parsing (hard dependency of `peak_compute`) |
| Disk | ≈ 2× the database size | PCOD database 1.7 GB, build peaks at about 4 GB |
| Memory | 8 GB class | chunked writes, does not fully consume the reflection table |

```bash
pip install numpy gemmi
```

### 3.2 CIF archives and the inorganic subset

The code does not download data. The chain is: **original archive → inorganic subset → `.xinDB`**.

**The source data is not in this repository**; as of 2026-10-23 it is kept in one place at
`D:/Project/XRD/sync-COD` (shared with
[PolyXRD-DataBase](https://github.com/PolyXRD/PolyXRD-DataBase); for data-side tools and layout notes
see [PolyXRD/sync-COD](https://github.com/PolyXRD/sync-COD)):

```
D:/Project/XRD/sync-COD/       ← data root, can be overridden with the XRD_DATA_ROOT environment variable
├── COD/        535,166 CIFs   # full-database original archive (level 3)
├── PCOD/     1,061,964 CIFs   # GRINSP predictions (level 2)
├── TCOD/       362,294 CIFs   # theoretical (level 3, 8-digit id)
├── inorg_cif/COD|PCOD|TCOD/    # ← step ① output: CIF copies judged inorganic (79,198 / 1,061,962 / 46,943)
└── inorg_manifest.sqlite       # ← step ① output: one-line-per-file manifest (1,959,421 rows, id_clash=0)
```

Step ① (inorganic classification = "no C, or containing C but no H" + copy + manifest):

```bash
python tools/scan_inorg.py --workers 6          # by default reads/writes the data root above in place
XRD_DATA_ROOT=/other/root python tools/scan_inorg.py --workers 6   # override when switching disks
```

- The criteria and id offsets (PCOD +2e7 / TCOD +1e8) are **entry-by-entry identical** to the three
  databases already published by this repository (the script is reused verbatim).
- The `inorg_cif/` used by the published databases is currently at
  `D:\Project\XRD\sync-COD\inorg_cif\` (1,188,103 CIFs, shared with polyDB). This repository does
  not carry this batch of data and does not need to — read-only reference is enough.
- During a build, `--cif-root` defaults to the `DEFAULT_INORG_*_ROOT` in
  `xinmatch/core/db/cif_scan.py`, which is derived from `XRD_DATA_ROOT` (default is
  `D:\Project\XRD\sync-COD\inorg_cif\*`). When switching machines, **point to it via the environment
  variable or `--cif-root`**; no need to change the code.

### 3.3 Building the databases (six databases, unified entry point)

`python -m tools.build_cod_xindb` is the only build entry point: `--source {COD,PCOD,TCOD}` ×
`--set {inorg,full}`, with the CIF root / `library_kind` / id offset derived automatically;
`tools/build_cod_inorg_xindb.py` is a thin compatibility wrapper around it.

```bash
# Three inorganic databases (default databases for interactive search, small and quick to build)
python -m tools.build_cod_xindb --source COD  --set inorg --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source PCOD --set inorg --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source TCOD --set inorg --data-date 20261020 --workers 14 -v

# Three full databases (exhaustive reference; takes hours, see manual §3.6)
python -m tools.build_cod_xindb --source COD  --set full --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source PCOD --set full --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source TCOD --set full --data-date 20261020 --workers 14 -v
```

- Output defaults to `<repo>/data/<kind>-<data_date>.xinDB`; `--out` specifies any path,
  `XINMATCH_XINDB_DIR` overrides only when explicitly set.
- `--limit N` for a trial build first (**always run a slice and get validate fully green before
  running the full build**); `--resume` resumes from a breakpoint (== equivalent to a single full
  run, verified); `--start` starts from a specified number range; `--workers` is unrelated to the
  output (1 vs 2 is byte-for-byte identical, with regression tests).
- **Use `--workers` 14 (optimal measured on this machine's 8C/16T on real full databases)**:
  TCOD 46,846 phases in 15.2 min (the old implementation with 14 workers took 16.9 min);
  COD with 8 workers ran 61.6 min, slightly slower than the historical 58 min with 14 workers —
  the 500-phase slice benchmark underestimates the benefit of more workers, so do not lower workers
  to 8 on that basis.
- The build hot path uses `build_pipeline.compute_phase()`, **parsing each CIF only once** (the old
  implementation parsed 3 times + built the structure twice); structure-factor chunking adapts to the
  atom count (256–8192 rows), and both optimizations **leave the output bit-for-bit unchanged**,
  raising whole-database throughput by about +10–15% (see manual §3.7 for details).
- The source defaults are `max_peaks=150`, `dmin=0.7703` (= λ/2); if you change the peak convention,
  change `peak_compute.PEAK_CONVENTION` at the same time, otherwise `meta` will not match reality and
  the validator will block it.
- `--kind full|inorg` is only a backward-compatibility alias for COD (giving it together with
  `--set` raises an error).

### 3.4 Measured scale (after the 2026-11-01 rebuild)

| Database | Phase count | Reflection count | Avg peaks/phase | File | Build time (14 workers benchmark) |
|---|---|---|---|---|---|
| COD-Inorg | **78,035** | **9,395,681** | 120.4 (min 2 / max 150) | **113.1 MiB** | 61.6 min (`--workers 8`; 14 workers should be faster) |
| PCOD-Inorg | **1,061,953** | **159,275,318** | 150.0 (min 20 / max 150) | **1694.1 MiB** | **13.3 h** (old implementation 12.7 h; external tasks contended during the run) |
| TCOD-Inorg | **46,846** | **5,874,790** | 125.4 (min 3 / max 150) | **68.0 MiB** | **15.2 min** (old implementation 16.9 min) |
| TCOD-Full | **361,701** | **53,099,238** | 146.8 (min 3 / max 150) | **594.3 MiB** | **12.1 h** (499.3 phases/min, 14 workers) |
| COD-Full | **532,038** | **77,458,154** | 145.6 (min 1 / max 150) | **879.4 MiB** | **36.4 h** (243.7 phases/min, 14 workers) |
| MP (external database) | **155,201** | **20,900,193** | 134.7 (min 8 / max 150) | **240.0 MiB** | **9.0 min** (17,318 phases/min, 14 workers) |

> The old implementation did about 1,590 phases/min single-process (COD slice); after the 2026-10
> single-pass parse + adaptive chunking optimization, whole-database throughput rose by about
> +10–15%; see [Build and acceptance handbook (Chinese)](docs/构建与验收手册.md) §3.7 for details.

---

## 4. Validation

```bash
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --json tmp/validate.json
python -m tools.validate_xindb data/PCOD-Inorg-20261020.xinDB --json tmp/validate-pcod.json
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --all-chunks   # full-chunk decode before release
```

Read-only, using only the public schema; criteria groups: `S` tables/indexes complete, `M` meta
complete and consistent with the file name, `C` chunk accounting (`Σpeak_n` = `Σn_peaks` =
`entry_count`, chunk count = `ceil(N/4096)`, `peak_off` contiguous), `E` cell physical quantities /
elements non-empty / `n_peaks ≤ 150` / **`E5` single-peak bad phases (180° boundary, placeholder
cell)** / **`E6` empty phases**, `P` sampled or full-chunk decoding + norm-1000 + zero padding +
`top⊆all`, `X` element inverted index value range / ascending order / self-consistency.

Before release also run **`tools/accept_xindb.py`** (source/payload/bad-phase acceptance common to
all 6 databases) and two special reviews:

```bash
python -m tools.accept_xindb data/COD-Inorg-20261020.xinDB \
    --manifest D:/Project/XRD/sync-COD/inorg_manifest.sqlite --json tmp/accept.json
python -m tools._compare_xindb_peaks.py <old database> <new database>          # only metadata changed, physics untouched
python -m tools._verify_resume_equivalence                     # --resume is equivalent to a single full run
```

> ✅ **The three Inorg databases were rebuilt with the fixed engine on 2026-11-01/02**:
> all 8 COD bad phases cleared (7 recovered real peaks, 1 phase without a cell removed, phase count
> 78,036→78,035, 9,395,681 reflections);
> all 7 PCOD bad phases cleared (159,275,318 reflections, only 17/1,061,953 payloads changed);
> TCOD had no bad phases and its peak payloads are byte-for-byte unchanged. All three databases'
> `validate_xindb --all-chunks` and `accept_xindb` are **0 FAIL**.
> **TCOD-Full and COD-Full have both been built and accepted**: TCOD-Full 361,701 phases /
> 53,099,238 reflections; COD-Full 532,038 phases / 77,458,154 reflections (36.4 h); the phases
> common with the corresponding Inorg databases are **byte-for-byte identical**.
> The earlier `20261020` Release (COD 8 / PCOD 7 bad phases) **was deleted on 2026-11-04**, leaving
> only the rebuilt `20261020-rebuild`.
> If you still have old Inorg zips on hand, do not use them; switch to the rebuilt version.
> See [`docs/ACCEPTANCE-20261101-rebuild.md`](docs/ACCEPTANCE-20261101-rebuild.md) for details.

---

## 5. Directory structure

```
.
├── _atomic_scattering_params.json Cromer-Mann coefficients (fastxrd finds it via tools/../)
├── data/                          **local artifact directory** (not committed): after the rebuild, 5 databases (3 Inorg + 2 Full) + build logs
├── pytest.ini  tests/             toolchain regression tests (synthetic CIFs/archives, independent of the 1.8 GiB real databases)
├── tmp/                           local trial builds / bench / acceptance JSON (not committed, .gitignore)
├── xinmatch/                      **build subset** (full program in PolyXRD/XinMatch)
│   ├── __init__.py                empty shell, to avoid dragging scipy/GUI into the build chain
│   └── core/db/
│       ├── cif_meta.py            CIF metadata extraction (chemical formula / space group / cell / density / elements)
│       ├── cif_scan.py            CIF archive scan + six-database source registry (SOURCE/SUBSET/id offsets)
│       ├── peak_compute.py        reflection calculation + peak convention constants (with Debye-Waller / boundary and placeholder-cell fixes)
│       ├── build_pipeline.py      build hot path: compute_phase() single-pass parse (meta + peak table)
│       ├── schema_v3.py           v3 DDL / required tables / required meta keys
│       ├── writer_v3.py           chunked columnar writer + element_inv rebuild
│       └── data/atomic_scattering_params.json   Cromer-Mann coefficients (99 elements)
├── tools/                         four-step chain, index in tools/README.md
│   ├── scan_inorg.py              ① inorganic classification + copies + manifest (XRD_DATA_ROOT determines the root directory)
│   ├── cif_common.py              ① implementation of the criteria (parse_elements / is_inorganic / formula_mass)
│   ├── fastcif.py  fastxrd.py     ① CIF parsing / reflection calculation dependencies
│   ├── build_cod_xindb.py         ② **only build entry point**: --source × --set (6 databases) + resume
│   ├── build_cod_inorg_xindb.py   ② thin compatibility wrapper (default inorg)
│   ├── _bench_build.py            ② build throughput benchmark (old two-step vs single-pass parse, workers/chunksize)
│   ├── validate_xindb.py          ③ read-only structure validator (S/M/C/E/P/X, including E5/E6)
│   ├── accept_xindb.py            ④ acceptance common to all 6 databases (id range / source / elements / payload / bad phases)
│   ├── _compare_xindb_peaks.py    ④ chunk-by-chunk peak-payload comparison (proves only metadata changed, physics untouched)
│   └── _verify_resume_equivalence.py ④ proves --resume is equivalent to a single full run
└── docs/
    ├── 构建与验收手册.md           **build/acceptance/release checklist (authoritative entry point)**
    └── xinDB格式规范-v3.md        field-level specification (container / addressing / segments / meta / validation)
```

`xinmatch/` is a **trimmed subset**: it keeps only the modules used for building databases, and
`__init__.py` is replaced with an empty shell (the full version does `import scipy`, etc., which the
build does not need). The remaining search/refinement/GUI modules are in the full repository.

Not committed (`.gitignore`): `data/ *.xinDB *.db *.sqlite tmp/ _release/ __pycache__/`.
This machine's `data/` currently holds the rebuilt 5 databases (COD-Inorg 113.1 MB /
TCOD-Inorg 68.0 MB / PCOD-Inorg 1694.1 MB / TCOD-Full 594.3 MB / COD-Full 879.4 MB) and
`data/_build-*.log`; the old databases in `data/_prev/` were deleted on 2026-11-04.

---

## 6. Known gaps

| Item | Current state | Impact |
|---|---|---|
| **Bad phases in the pre-rebuild old databases** | old `20261020` Release: COD 8, PCOD 7, TCOD 0 | built with the old peak engine; **after the 2026-11-01/02 rebuild all three Inorg databases are at 0**. The old Release and its Inorg zips were deleted on 2026-11-04, leaving only the rebuilt version |
| Semi-collapsed boundary phases (WARN) | COD 2,251 / PCOD 5,598 (after rebuild) / TCOD 1,670 entries: the strongest line falls at d≈λ/2 (2θ→180°), other peaks exist (intensity sum ≥5/1000) | a **peak-convention edge phenomenon** caused by a real reflection + LP diverging at 180°, not a spurious peak, and **rebuilding will not eliminate it**; a real fix requires changing the peak convention (raising `DEFAULT_DMIN` to 0.78–0.79 Å or dropping 2θ>179°, a v4 decision). Currently only a WARN, does not block release |
| Degenerate cell | no longer present in the rebuilt COD database (source CIF of `96-210-4629` has no cell → removed during build) | — |
| No space group | 126 entries in the COD database (127 before the rebuild, `96-210-4629` removed; `E3` only gives a hint) | does not affect peak positions; such phases are excluded when filtering by space group |
| Full-database release | **TCOD-Full and COD-Full have both been built and accepted** (TCOD 361,701 phases / 53,099,238 reflections / 594.3 MiB; COD 532,038 phases / 77,458,154 reflections / 879.4 MiB; both validate/accept 0 FAIL, phases common with the corresponding Inorg databases byte-for-byte identical); `PCOD-Full` is equivalent to `PCOD-Inorg` (PCOD is all inorganic), so it is not built | once built, XinMatch auto-discovers them via the built-in `xindb_tcod_full` / `xindb_cod_full` rows; in TCOD-Full 48,732 phases and COD-Full 81,284 phases have elements exceeding the 63-bit `elem_mask` (`element_inv` is still complete) |
| Full-database release | **building already supports** 6 databases (`--source × --set`); the `20261020` release published only the three Inorg databases | organic/full-database search requires first building the `*-Full` databases (takes hours, see build manual §3.6) |
| Narrow PCOD element range | only 19 elements (O/Si/P/Al/V/S/B/F…) | the objective composition of the GRINSP prediction set, not data lost during the build |
| `top_i` is float32 | intensity uses `<f4` rather than quantized integers | size traded for precision; under norm-1000 the significant digits are sufficient |

---

## 7. Data sources and licensing

- **Code**: follows XinMatch's dual license (academic non-commercial / commercial by separate
  agreement), see `LICENSE`, `LICENSE-CN`.
- **Data**: CIF content copyright belongs to COD and the authors of the original literature. When
  using this repository's artifacts, please comply with
  [COD's terms of use](http://www.crystallography.net/cod/) and cite COD and the corresponding
  original literature; PCOD (GRINSP predictions) and TCOD (theoretical) entries have their own
  original sources. This repository distributes only **derived indexes**
  (reflection positions/intensities/metadata); it does not rewrite or replace the original CIFs.
- Intensities are **computed values**, including Debye-Waller but with no instrumental/sample
  broadening, so they still deviate from experimental cards (ICDD PDF); sufficient for phase
  identification; for quantitative analysis you must calibrate yourself.

## 8. Version history

| Date | Event |
|---|---|
| 2026-10-21 / 22 | `COD-Inorg`, `TCOD-Inorg`, `PCOD-Inorg` databases built (xinDB v3) |
| 2026-10-23 | this repository created and the `20261020` release published (code + three zips); COD database `journal_mode` normalized from WAL to DELETE, `integrity_check` re-checked and passed |
| 2026-10-23 | the three databases moved from `PolyXRD-DataBase/xindb/` to this repository's `data/`; added inorganic classification (①) and two special-review (④) tools, making the chain self-contained |
| 2026-10-23 | **source data moved out** → `D:/Project/XRD/sync-COD` (shared by both databases), `XRD_DATA_ROOT` default changed to that path; `build_dbs.py` stripped into the pure-function module `cif_common.py` (old 3COD sqlite packager deleted) |
| 2026-10-23 | data-side tools spun off into the public repository [**PolyXRD/sync-COD**](https://github.com/PolyXRD/sync-COD) (data still not committed) |
| 2026-10-24 | **build/acceptance tools reworked**: unified six-database build entry point (`--source × --set`, including PCOD/TCOD Full); added `build_pipeline.compute_phase()` single-pass parse; added `accept_xindb.py` generic acceptance and `E5/E6` bad-phase checks; established `tests/` regressions (synthetic CIFs, independent of real databases); added [Build and acceptance handbook (Chinese)](docs/构建与验收手册.md) and recorded the convention "databases are built in this repository, CIF sources are in `sync-COD`" |
| 2026-11-01 | **rebuilt the three Inorg databases with the fixed peak engine** (source snapshot still `20261020`): COD 78,036→**78,035** phases / 9,395,681 reflections (8 bad phases cleared); TCOD 46,846 phases, peak payloads byte-for-byte identical to the old database; PCOD 1,061,953 phases / 159,275,318 reflections (7 bad phases cleared); all three validate/accept **0 FAIL**. See [`docs/ACCEPTANCE-20261101-rebuild.md`](docs/ACCEPTANCE-20261101-rebuild.md) |
| 2026-11-02 | began building **TCOD-Full** (362,294 source CIFs) and **COD-Full** (535,166 source CIFs); `PCOD-Full` skipped because PCOD is entirely inorganic and equivalent to `PCOD-Inorg` |
| 2026-11-02 | **TCOD-Full built and accepted**: 361,701 phases / 53,099,238 reflections / 623,173,632 B / 12.1 h; validate 0F/2W, accept 0F/1W; 46,846/46,846 phases common with TCOD-Inorg byte-for-byte identical; self-recall 10/10. COD-Full build in progress |
| 2026-11-04 | **COD-Full built and accepted**: 532,038 phases / 77,458,154 reflections / 922,120,192 B / **36.4 h**; validate 0F/2W, accept 0F/2W (`bad_phase_ids=0`); 78,035/78,035 phases common with COD-Inorg byte-for-byte identical; self-recall 10/10. All 5 databases (3 Inorg + 2 Full) complete, PCOD-Full skipped as equivalent |
| 2026-11-04 | **MP-20261031 external database built and accepted** (Materials Project, 155,204 CIFs): 155,201 phases / 20,900,193 reflections / 251,633,664 B / 9.0 min; `mp-<n>` ids, validate 0F/1W, accept 0F/1W (id set exactly matches the archive); self-recall 20/20. See [MP database build and acceptance (Chinese)](docs/MP-库构建与验收.md) |
