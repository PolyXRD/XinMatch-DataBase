# XinMatch-DataBase

**XinMatch 参考物相库（`.xinDB` v3 格式）的构建代码、格式规范与发布产物。**

本仓库只放**代码与文档**；数据库本体以 [Releases](../../releases) 里的压缩包分发
（单库 68 MB / 87 MB / 1.3 GB，不入 git 历史）。数据源为
[Crystallography Open Database (COD)](http://www.crystallography.net/cod/) 及其预测/理论子库
PCOD / TCOD，本仓库不做结构精修、不改写原始 CIF，只做反射计算与索引打包。

> 姊妹仓库 [PolyXRD-DataBase](https://github.com/PolyXRD/PolyXRD-DataBase) 发布的是同源数据的
> **polyDB** 格式（PolyXRD 用）。**两种格式口径不同，不可互换**（见 §2）。

---

## 1. 发布产物（Release `20261020`）

| 库 | entry_id 段 | 相数 | 反射数 | 块数 | 元素 | 解压后 | zip |
|---|---|---|---|---|---|---|---|
| **COD-Inorg** | `96-100-*` | 78,036 | 9,394,774 | 20 | 96 种 | 118,710,272 B（113.2 MiB） | **87.2 MiB**（91,484,834 B） |
| **PCOD-Inorg** | `96-210-*` | 1,061,953 | 159,274,265 | 260 | 19 种 | 1,776,353,280 B（1694.1 MiB） | **1.35 GiB**（1,447,251,707 B） |
| **TCOD-Inorg** | `96-110-*` | 46,846 | 5,874,790 | 12 | 91 种 | 71,352,320 B（68.0 MiB） | **54.1 MiB**（56,694,762 B） |
| **合计** | — | **1,186,835** | **174,543,829** | 292 | — | 1.83 GiB | ≈1.52 GiB |

- 三个库**互不包含，按需装**（XinMatch 支持多库挂载）。COD 实验库是日常交互检索的主力；
  PCOD（GRINSP 预测）与 TCOD（理论）用于"实验库查不到时兜底"。
- 命名里的 `20261020` 是**源数据快照日**（`meta.data_date`），不是构建日；
  构建时刻见 `meta.build_time`（COD 2026-10-21 17:25、TCOD 2026-10-22 23:29、PCOD 2026-10-22 23:46）。
- 分发版全部 `journal_mode=DELETE`，无 `-wal`/`-shm` 侧车，可直接放只读介质
  （NAS 只读共享、光盘镜像）而不必担心"读也要写权限"。
- 校验值见 Release 附带的 `SHA256SUMS.txt`。

```bash
unzip COD-Inorg-20261020.xinDB.zip                    # 得到 COD-Inorg-20261020.xinDB
python -m tools.validate_xindb COD-Inorg-20261020.xinDB   # 可选：自检 25 项
```

> 本机已把三个库放在仓库的 `data/` 下（`data/COD-Inorg-20261020.xinDB` 等，不入 git），
> XinMatch 的库注册表 `~/.xinmatch/libraries.json` 已指向该目录。

> ⚠️ 文件名必须是 `<library_kind>-<data_date>.xinDB`，否则 `validate_xindb` 的
> `M3.identity` 会判 FAIL（它拿文件名和 `meta` 对账）。改名后重新建库即可，不必修 meta。

---

## 2. `.xinDB` 是什么

**单文件 SQLite**，把"预筛矩阵"和"全反射表"分开存，反射数据只存一份：

```
COD-Inorg-20261020.xinDB
├── meta          自描述 key/value（峰口径、建库参数、追溯）
├── entries       每相一行（化学式/空间群/晶胞/密度/元素 + CSR 偏移）
├── peak_chunks   列式反射块（每 4096 相一块）
└── element_inv   元素倒排（element → 0 起相序号，uint32 数组）
```

| 段 | 形状 | 类型 | 用途 |
|---|---|---|---|
| `top_d`/`top_i` | `n × 32` | `<f4` | 预筛：与测量峰窗口做 NumPy 交叠，**只剪枝不排序** |
| `all_d`/`all_i` | `Σpeak_n`（CSR） | `<f4` | 打分：候选相的全部反射 |
| `element_inv.eids` | 变长 | `<u4` | 升序去重、0 起相序号 |

**寻址铁律**：`entries.rowid` 显式从 0 写入 = 0 起相序号；`element_inv.eids` 与
`peak_chunks` 块号都基于它，**不能**用 `entry_id` 反推。
完整字段级规范见 [`docs/xinDB格式规范-v3.md`](docs/xinDB格式规范-v3.md)。

### 峰口径（写死进 `meta.peak_convention`）

```
CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000
```

λ=1.5406 Å（Cu Kα1）｜d 0.7703–30.0 Å｜等价反射按 2θ 合并（含多重性）｜
**含 Debye-Waller**（`B = 8π²U_iso`，逐位点）｜强度归一 1000｜每相最多 150 条反射｜预筛宽 32

| | xinDB（本仓库） | polyDB（PolyXRD-DataBase） |
|---|---|---|
| 容器 | 单文件 SQLite | 目录（8 个段文件 + mmap） |
| d 范围 | 0.7703–30 Å | 1.0–26 Å（2θ 5–100.76°） |
| 温度因子 | **有**（Debye-Waller） | 无 |
| 每相峰数 | 150 | 100（按强度取前 K） |
| 定位 | 检索 + 精修一体化 | GUI 交互检索 |

> 同一份 CIF 在两个库里算出的绝对强度体系不同。**混用会得到系统性偏差**，
> 不要拿 xinDB 的强度去校验 polyDB，反之亦然。

---

## 3. 从零构建

### 3.1 环境

| 依赖 | 要求 | 说明 |
|---|---|---|
| Python | 3.10+ | 实测 3.10.11 |
| numpy | 必需 | 反射计算与列式打包 |
| gemmi | **必需** | CIF 解析（`peak_compute` 硬依赖） |
| 磁盘 | ≈ 2× 库体积 | PCOD 库 1.7 GB，构建峰值约 4 GB |
| 内存 | 8 GB 级 | 分块写入，不吃满反射表 |

```bash
pip install numpy gemmi
```

### 3.2 CIF 归档与无机子集

代码不下载数据。链路是：**原始归档 → 无机子集 → `.xinDB`**。

**源数据不在本仓库**，2026-10-23 起统一放在 `D:/Project/XRD/sync-COD`（与
[PolyXRD-DataBase](https://github.com/PolyXRD/PolyXRD-DataBase) 共用，
数据侧工具与布局说明见 [PolyXRD/sync-COD](https://github.com/PolyXRD/sync-COD)）：

```
D:/Project/XRD/sync-COD/       ← 数据根，可用环境变量 XRD_DATA_ROOT 覆盖
├── COD/  PCOD/  TCOD/          # 原始归档（层级 3 / 2 / 3，脚本按此解析 id）
├── inorg_cif/COD|PCOD|TCOD/    # ← 第 ① 步产物：判为无机的 CIF 副本
└── inorg_manifest.sqlite       # ← 第 ① 步产物：每文件一行清单
```

第 ① 步（无机判定 = "不含 C，或含 C 但不含 H" + 复制 + 清单）：

```bash
python tools/scan_inorg.py --workers 6          # 缺省就地读写上面的数据根
XRD_DATA_ROOT=/other/root python tools/scan_inorg.py --workers 6   # 换盘时覆盖
```

- 判据、id 偏移（PCOD +2e7 / TCOD +1e8）与本仓库已发布的三库**逐条一致**（脚本原样沿用）。
- 已发布库用的 `inorg_cif/` 现状就在 `D:\Project\XRD\sync-COD\inorg_cif\`（1,188,103 个 CIF，
  与 polyDB 共用）。本仓库不带这批数据，也不必带 —— 只读引用即可。
- 建库时 `--cif-root` 缺省指向 `xinmatch/core/db/cif_scan.py` 的 `DEFAULT_INORG_*_ROOT`，
  后者由 `XRD_DATA_ROOT` 推导（缺省即 `D:\Project\XRD\sync-COD\inorg_cif\*`）。
  换机器**用环境变量或 `--cif-root` 指过去**，不必改代码。

### 3.3 建库

```bash
python -m tools.build_cod_inorg_xindb --source COD  --data-date 20261020 --workers 14 -v
python -m tools.build_cod_inorg_xindb --source PCOD --data-date 20261020 --workers 14 -v
python -m tools.build_cod_inorg_xindb --source TCOD --data-date 20261020 --workers 14 -v
```

- 输出默认 `$XINMATCH_XINDB_DIR/<kind>-<data_date>.xinDB`（未设环境变量则落 `<repo>/data/`），
  可用 `--out` 指定。
- `--limit N` 先试产（实测 300 条约 11 s）；`--resume` 断点续跑；`--start` 从指定号段起。
- **先跑切片再跑全量**：`--limit 500` 出来的库必须 `validate_xindb` 全 PASS 才值得跑全量。
- 源码默认 `max_peaks=150`、`dmin=0.7703`（= λ/2），要改峰口径请同时改
  `peak_compute.PEAK_CONVENTION`，否则 `meta` 与实际不符、校验器会拦。

### 3.4 实测规模（本机 4 核 / 机械盘）

| 库 | 相数 | 反射数 | 平均峰/相 | 文件 |
|---|---|---|---|---|
| COD-Inorg | 78,036 | 9,394,774 | 120.4（min 1 / max 150） | 113.2 MiB |
| PCOD-Inorg | 1,061,953 | 159,274,265 | 150.0（全为上限截断） | 1694.1 MiB |
| TCOD-Inorg | 46,846 | 5,874,790 | 125.4（min 3 / max 150） | 68.0 MiB |

本仓库代码的构建速率实测约 **1,590 相/分钟（单进程，COD 切片）**，`--workers` 线性摊薄。

---

## 4. 校验

```bash
python -m tools.validate_xindb COD-Inorg-20261020.xinDB            # 25 项
python -m tools.validate_xindb PCOD-Inorg-20261020.xinDB --json out.json
```

只读、只用公开 schema，六组判据：`S` 表/索引齐全，`M` meta 完整且与文件名一致，
`C` 块记账（`Σpeak_n` = `Σn_peaks` = `entry_count`、块数 = `ceil(N/4096)`、`peak_off` 连续），
`E` 晶胞物理量 / 元素非空 / `n_peaks ≤ 150`，`P` 抽样块解码 + norm-1000 + 零填充 + `top⊆all`，
`X` 元素倒排值域/升序/自洽。

本次三个发布库的实测结果：

| 库 | 检查项 | FAIL | WARN |
|---|---|---|---|
| COD-Inorg | 26 | 0 | **1**（`E1.lattice`：1 条退化晶胞，如 `96-210-4629`） |
| PCOD-Inorg | 25 | 0 | 0 |
| TCOD-Inorg | 25 | 0 | 0 |

另做 SQLite `PRAGMA integrity_check` = ok、`entries` 行数与 `meta.entry_count` 逐库相符、
`element_inv` 覆盖元素数与 `entries.elements` 实际元素集合完全一致（缺失 0）。

发布前另可跑两份专项复核：`_compare_xindb_peaks.py`（逐 `peak_chunks` 行比对 BLOB 摘要，
证明"只改元数据没动物理"）与 `_verify_resume_equivalence.py`（证明 `--resume` 等价于一次跑完）。

---

## 5. 目录结构

```
.
├── _atomic_scattering_params.json Cromer-Mann 系数（fastxrd 按 tools/../ 找它）
├── data/                          **本地产物目录**（不入仓）：三个 .xinDB + _backup/
├── xinmatch/                      **构库子集**（完整程序见 PolyXRD/XinMatch）
│   ├── __init__.py                空壳，避免把 scipy/GUI 拖进构建链路
│   └── core/db/
│       ├── cif_meta.py            CIF 元数据提取（化学式/空间群/晶胞/密度/元素）
│       ├── cif_scan.py            CIF 归档扫描与 id 分配（含 3 库偏移）
│       ├── peak_compute.py        反射计算 + 峰口径常量（含 Debye-Waller）
│       ├── schema_v3.py           v3 DDL / 必需表 / 必需 meta 键
│       ├── writer_v3.py           分块列式写入器 + element_inv 重建
│       └── data/atomic_scattering_params.json   Cromer-Mann 系数（99 元素）
├── tools/                         四步链路，索引见 tools/README.md
│   ├── scan_inorg.py              ① 无机分类 + 副本 + 清单（XRD_DATA_ROOT 决定根目录）
│   ├── cif_common.py              ① 的判据实现（parse_elements / is_inorganic / formula_mass）
│   ├── fastcif.py  fastxrd.py     ① 依赖的 CIF 解析 / 反射计算
│   ├── build_cod_inorg_xindb.py   ② 三库构建入口（**主入口**）
│   ├── build_cod_xindb.py         ② 底层构建器（分块写入 / 断点续跑 / 记账）
│   ├── validate_xindb.py          ③ 只读验收器（S/M/C/E/P/X 六组）
│   ├── _compare_xindb_peaks.py    ④ 峰载荷逐块比对（证明只改元数据没动物理）
│   └── _verify_resume_equivalence.py ④ 证明 --resume 等价于一次跑完
└── docs/
    └── xinDB格式规范-v3.md        字段级规范（容器/寻址/各段/meta/校验）
```

`xinmatch/` 是**裁剪子集**：只保留构库用到的 5 个模块，`__init__.py` 换成空壳
（完整版会 `import scipy` 等，构库不需要）。其余检索/精修/GUI 模块见完整版仓库。

不入仓（`.gitignore`）：`data/ *.xinDB *.db *.sqlite _release/ __pycache__/`。
本机 `data/` 里放着三个已发布库（118.7 MB / 1.78 GB / 71.4 MB）与 `_backup/`
（COD 库改 meta 之前的两个历史副本，各 118 MB，确认无用可直接删）。

---

## 6. 已知缺口

| 项 | 现状 | 影响 |
|---|---|---|
| 退化晶胞 | COD 库 1 条（`96-210-4629`：a=b=c=1、体积=1） | 该相检索会给出无意义 d；占 0.0013%，源 CIF 本身的占位值 |
| 无空间群 | COD 库 127 条（`E3` 仅提示） | 不影响峰位；按空间群过滤时这类相会被排除 |
| 仅无机子集 | 只发布 Inorg；纯有机相不在三个库里 | 有机样品检索需另建 Full 库 |
| PCOD 元素面窄 | 只有 19 种元素（O/Si/P/Al/V/S/B/F…） | GRINSP 预测集的客观构成，不是建库丢数据 |
| `top_i` 为 float32 | 强度用 `<f4` 而非量化整数 | 体积换精度；norm-1000 下有效位足够 |

---

## 7. 数据来源与许可

- **代码**：沿 XinMatch 的双许可（学术非商业 / 商业另议），见 `LICENSE`、`LICENSE-CN`。
- **数据**：CIF 内容版权归 COD 及各原始文献作者。使用本仓库产物请遵守
  [COD 的使用条款](http://www.crystallography.net/cod/) 并引用 COD 与相应原始文献；
  PCOD（GRINSP 预测）、TCOD（理论）条目另有其原始出处。本仓库只分发**派生索引**
  （反射位置/强度/元数据），不改写也不替代原始 CIF。
- 强度为**计算值**，含 Debye-Waller 但无仪器/样品展宽，与实验卡片（ICDD PDF）仍有偏差；
  用于物相识别足够，用于定量须自行标定。

## 8. 版本

| 日期 | 事件 |
|---|---|
| 2026-10-21 / 22 | `COD-Inorg`、`TCOD-Inorg`、`PCOD-Inorg` 三库建成（xinDB v3） |
| 2026-10-23 | 本仓库建立并发布 `20261020` release（代码 + 三个 zip）；COD 库 `journal_mode` 由 WAL 归一为 DELETE，`integrity_check` 复核通过 |
| 2026-10-23 | 三个库从 `PolyXRD-DataBase/xindb/` 迁到本仓库 `data/`；补入无机分类（①）与两份专项复核（④）工具，链路自足 |
| 2026-10-23 | **源数据迁出** → `D:/Project/XRD/sync-COD`（两库共用），`XRD_DATA_ROOT` 缺省改为该路径；`build_dbs.py` 剥成纯函数模块 `cif_common.py`（旧 3COD sqlite 打包器删除） |
| 2026-10-23 | 数据侧工具独立成公开仓库 [**PolyXRD/sync-COD**](https://github.com/PolyXRD/sync-COD)（数据仍不入库） |
