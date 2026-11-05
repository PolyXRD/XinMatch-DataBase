# `.xinDB` v3 格式规范（XinMatch 参考库）

> 状态：**现行（schema_version = 3）**　实现：`xinmatch/core/db/schema_v3.py`、`writer_v3.py`
> 设计目标：把"每相 × 每条反射一行"的旧 §2 布局，改成**列式分块**存储 —— 反射数据只存一份，
> 预筛矩阵与全反射表分开放。实测原因见 `schema_v3.py` 顶部注释（旧库 95% 体积是峰表 + 3 个索引，
> 同一份数据被存了三遍）。

---

## 1. 容器

**单文件 SQLite**（扩展名 `.xinDB`），**无侧车、无外部依赖**：

```
COD-Inorg-20261020.xinDB
├── meta          自描述 key/value（建库参数、峰口径、追溯信息）
├── entries       每相一行（标量元数据 + 指向反射区的 CSR 偏移）
├── peak_chunks   分块列式反射数据（每 4096 相一个块）
└── element_inv   元素倒排索引（element → 0 起相序号数组）
```

- 索引 4 个（`space_group` / `crystal_system` / `formula` / `tier`），均为 B-tree 普通索引。
- 日志模式：**`journal_mode = DELETE`**（分发版）。不要用 WAL 分发 —— 只读介质
  （NAS 只读共享、光盘镜像）上 WAL 库打不开，因为读也需要 `-shm` 写权限。
- 全部数值 little-endian（`<f4` / `<u4`）。

---

## 2. 峰口径（`meta.peak_convention`，冻结字符串）

```
CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000
```

| 项 | 取值 | 说明 |
|---|---|---|
| 波长 λ | 1.5406 Å | Cu Kα1 |
| d 范围 | **0.7703 – 30.0 Å** | 下限 = λ/2（2θ = 180°），上限 30 Å |
| 等价反射 | 按 2θ 合并，计入多重性 | — |
| 温度因子 | **有**，`T = exp(-B s²)`，`B = 8π²U_iso` | CIF 给了 `U_iso` 才生效，逐位点不同的 B |
| 散射因子 | Cromer-Mann `f(s) = Z - 41.78214 s² Σ aᵢ exp(-bᵢ s²)` | 系数在 `xinmatch/core/db/data/atomic_scattering_params.json` |
| 强度归一 | I/I₀ × 1000 | 每相最强峰 = 1000 |
| 每相最多反射 | `max_peaks = 150` | 超出按强度截断 |
| 预筛矩阵宽 | `top_k = 32` | 每相最强 32 条 |

> ⚠️ **与 polyDB 口径不同，不可混用**：polyDB 是 `d ≥ 1.0 Å`、2θ 5–100.76°、**无**温度因子、
> 按强度取前 100 条（见 PolyXRD-DataBase 仓库）。同一份 CIF 在两个库里算出的强度体系不同，
> 库里各自的 `meta` 会记下口径，reader 发现未知口径必须拒绝加载而不是混算。

### 2.1 实测峰表（`meta.peak_origin = observed`）

上面那条冻结串描述的是**由 CIF 算峰**的库。ICDD PDF-2 这类**实测卡片库**没有算峰这一步，
峰表是几十年文献里测出来的相对强度，套用计算口径就是让文件说谎。规范为此留了一个显式的
开关：`meta.peak_origin`。

| 项 | `computed`（缺省） | `observed` |
|---|---|---|
| 谁用 | COD / PCOD / TCOD / MP | ICDD PDF-2 |
| `peak_convention` | 必须是 §2 的冻结串 | 必须以 `observed-` 开头，且含 `norm1000` |
| 强度尺度 | I/I₀ × 1000 | **同左**（不共用一把尺子，两类库的打分权重就不可比） |
| `entries` 里 NULL 晶胞 | 构建缺陷 | 数据事实（PDF-2 有 18.25% 的老卡片本就无晶胞） |
| 峰数 > `max_peaks=150` | 引擎安全阀被绕过 | 数据事实（PDF-2 有 14.9% 超 150，最多 199） |
| d 范围 | 由 §2 决定 | 源数据决定（PDF-2 实测 0.36–99.2 Å） |

`peak_origin` 缺省按 `computed` 解释，所以既有库的判据不受影响。读取侧**不需要**关心这个
字段（`reader_v3` 只校验 `schema_version`），只有验收器按它切换判据。

---

## 3. `entries`（每相一行）

```sql
CREATE TABLE entries (
    entry_id       TEXT PRIMARY KEY,   -- "96-XXX-XXXX…" 见 §6
    formula        TEXT,
    space_group    TEXT,
    crystal_system TEXT,
    mineral_name   TEXT,
    elements       TEXT,               -- 逗号分隔符号串，如 "Al,Ca,O,Si"
    elem_mask      INTEGER,            -- IUPAC 元素位图，仅低 63 位可存（见 §7）
    a REAL, b REAL, c REAL,
    alpha REAL, beta REAL, gamma REAL,
    volume  REAL,
    density REAL,
    year    INTEGER, author TEXT, journal TEXT, quality TEXT,
    n_peaks   INTEGER,                 -- 该相反射数（≤ 150）
    peak_off  INTEGER,                 -- 首条反射在 all_d/all_i 中的全局下标
    peak_n    INTEGER,                 -- = n_peaks（冗余，读取时校验用）
    tier      INTEGER DEFAULT 1        -- 1=实验 2=理论/预测
);
```

**寻址铁律**：`rowid` **显式从 0 开始写入**，即 `rowid` 就是 0 起相序号。
`element_inv.eids` 存的就是这个 `rowid`，`peak_chunks` 的块号 = `rowid // chunk_size`。
**不能用 `entry_id` 反推序号**（COD 编号从 1,000,000 起，PCOD/TCOD 还带偏移）。

---

## 4. `peak_chunks`（列式反射块）

```sql
CREATE TABLE peak_chunks (
    chunk_id INTEGER PRIMARY KEY,
    n        INTEGER,   -- 本块实际相数（最后一块可能 < chunk_size）
    top_d BLOB, top_i BLOB,   -- 预筛矩阵 (n, 32) <f4
    all_d BLOB, all_i BLOB    -- 全反射 CSR (Σpeak_n,) <f4
);
```

| 段 | 形状/长度 | 类型 | 用途 |
|---|---|---|---|
| `top_d` / `top_i` | `n × 32`，行内按强度降序，空槽填 0 | `<f4` | **预筛**：一次载入，与测量峰窗口做 NumPy 交叠，生成候选，不碰全反射数据 |
| `all_d` / `all_i` | `Σ peak_n`（CSR） | `<f4` | **打分**：候选相的全部反射，按 `peak_off/peak_n` 切片取出 |

- `chunk_size = 4096`；块数 = `ceil(entry_count / 4096)`。
- 第 `i` 相（`rowid = i`）的第 `k` 条反射 = `all_d[peak_off + k]`。
- `top_*` **只用于剪枝，不参与排序**：最终名次一律用 `all_*` 全反射重算，
  避免"预筛截断改变排序"。
- float32 的 d 精度在相关量级约 1e-7 Å，比 0.01–0.02 Å 的匹配窗低 3 个数量级，
  不会改变匹配判定（有 parity 测试背书）。

---

## 5. `element_inv`（元素倒排）

```sql
CREATE TABLE element_inv (element TEXT PRIMARY KEY, eids BLOB);
```

- `eids` = **升序、去重、little-endian `uint32`** 的相序号数组（`<u4`），值域 `0 ≤ eid < entry_count`。
- 覆盖**全部**元素（含 `elem_mask` 装不下的重元素），用于"允许元素池 ⊆ 相元素集"的精确过滤。

---

## 6. `entry_id` 编号规则

```
96-<移位后数字的前 3 位>-<其余位>
```

| 库 | 原始号 | 加偏移 | entry_id 示例 |
|---|---|---|---|
| COD | 1,000,005 | +0 | `96-100-0005` |
| PCOD | 1,000,000 | +20,000,000 | `96-210-00000` |
| TCOD | 10,000,001 | +100,000,000 | `96-110-000001` |
| MP（外部库） | 862,690 | +0 | `mp-862690` |

偏移不是装饰：PCOD 与 COD 共用 7 位数字空间，不偏移会让 `96-100-0000` 在两库中互相遮蔽。

**外部库的 id 约定**：COD/PCOD/TCOD 之外的库可以保留自己的稳定标识，不必套 `96-` 前缀。
Materials Project 库使用其原生 `mp-<number>`（`meta.library_kind = "MP"`、`data_date` 与文件名一致），
号段 `[0, 20,000,000)`；`accept_xindb` 对 MP 要求 **id 集合与归档 `# mp-id:` 注释集合完全相等**
（构建/验收见 [`MP-库构建与验收.md`](MP-库构建与验收.md)，代码 `mp_scan.py` / `build_mp_xindb.py`）。

---

## 7. `meta`（自描述）

必需键（`schema_v3.META_KEYS`，缺一个 reader 就应拒绝）：

| 键 | 示例 | 说明 |
|---|---|---|
| `schema_version` | `3` | 变化即不兼容 |
| `peak_convention` | `CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000` | §2 冻结串 |
| `top_k` | `32` | 预筛矩阵宽 |
| `chunk_size` | `4096` | 每块相数 |
| `dtype` | `f32` | 反射数组类型 |
| `entry_count` | `78036` | 相数 |
| `build_time` | `2026-10-21T17:25:59+08:00` | 计算时刻 |
| `cif_source` | `…\inorg_cif\COD` | CIF 归档根 |
| `build_tool` | `tools.build_cod_inorg_xindb` | 构建器模块名 |
| `library_kind` | `COD-Inorg` | 与文件名前缀一致 |

附加键（本仓库产物均带）：

| 键 | 示例 | 说明 |
|---|---|---|
| `data_date` | `20261020` | **源数据快照日**，与 `build_time` 分开记 |
| `tier_default` | `1` | 建库时写入的 tier 默认值 |
| `element_bits` | `H:0;He:1;…` | 元素符号 → 位号（`iupac128-v1` 方案） |
| `element_mask_bits` | `63` | `entries.elem_mask` 实际可存位宽（SQLite INTEGER 为有符号 64 位） |
| `element_mask_full_bits` | `128` | 方案总位宽 |
| `element_mask_scheme` | `iupac128-v1` | 位图方案名 |
| `element_mask_unmapped` | `Ac:88;Au:78;…` | 位号 ≥ 63、装不进 `elem_mask` 的元素（仍完整存在于 `elements` 与 `element_inv`） |

实测库（`peak_origin = observed`）另带：

| 键 | 示例 | 说明 |
|---|---|---|
| `peak_origin` | `observed` | 峰是实测值而非计算值，见 §2.1 |
| `peak_max` | `199` | 实际最大峰数；验收器用它自查，替代引擎的 `max_peaks=150` 上限 |
| `peak_d_min` / `peak_d_max` | `0.36` / `99.2` | 实际 d 范围（源数据决定，不是 §2 的 0.7703–30） |
| `dataset` | `PDF2-2004` | 数据集标识 |
| `raw_file` / `raw_bytes` / `raw_sha256` | `C:\pdf2\pdf2.dat` / `587999120` / `e25b521c…` | 原始二进制溯源 |
| `name_policy` | `mineral>compound` | `mineral_name` 列的填充策略（PDF-2 只有 10.8% 有矿物名） |
| `reference_policy` | `coden-year-author (卷/页无列位，未落库)` | 文献字段的拆分与丢失声明 |

---

## 8. 校验（`tools/validate_xindb.py` / `tools/accept_xindb.py`）

只读、只用公开 schema，可对"别处建好的库"复核。`validate_xindb` 分组：

| 组 | 判据 |
|---|---|
| `S` schema | 3 张必需表 + 4 个索引存在 |
| `M` meta | 必需键齐全；`peak_convention` 与 `peak_origin` 自洽（`computed` 必须等于代码里的常量，`observed` 必须自报 `observed-*`+`norm1000`）；`library_kind`/`data_date` 与文件名一致 |
| `C` chunk 记账 | `Σpeak_chunks.n` = `Σentries.n_peaks` = `entry_count`；块数 = `ceil(N/4096)`；`peak_off` 阶梯连续 |
| `E` entry 合理性 | 晶胞/体积物理量合法；每相至少一个元素；`peak_n == n_peaks`；`n_peaks ≤ 150` |
| `E5` 单峰坏相 | 穷举每个 chunk 的 `top_d[:,0]`（最强反射）：**最强反射 d 落在 λ/2 边界窗口、且其余峰强度和 < 5/1000（整谱已塌缩）→ FAIL**；`n_peaks==1` 且 d≈1.0 Å（gemmi 占位晶胞）→ FAIL；其它单峰正常 d → WARN |
| `E6` 空相 | `n_peaks <= 0` 或 NULL → FAIL |
| `P` 反射载荷 | 抽样块解码形状正确；反射有限且为正；强度在 norm-1000 尺度；`n_peaks` 之后是零填充；`top` 块是全反射的子集；`--all-chunks` 可全块解码 |
| `X` 元素倒排 | `eids` 值域合法、严格升序去重；抽样相能在自己 `elements` 的每个元素下找到自己 |

```bash
python -m tools.validate_xindb COD-Inorg-20261020.xinDB --json out.json
python -m tools.validate_xindb COD-Inorg-20261020.xinDB --all-chunks   # 发布前全块
```

`tools/accept_xindb.py` 是面向发布产物的**通用验收**（6 库一条命令）：全量 id 号段、
来源归属（源目录 / `inorg_manifest.sqlite` 抽查）、`element_inv` 全量一致性、
峰载荷（默认全量，可抽样）、坏相穷举；`FAIL` 退出码 1。构建/验收全流程见
[`构建与验收手册.md`](构建与验收手册.md)。
