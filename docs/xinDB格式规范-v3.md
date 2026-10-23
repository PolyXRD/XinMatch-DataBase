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

偏移不是装饰：PCOD 与 COD 共用 7 位数字空间，不偏移会让 `96-100-0000` 在两库中互相遮蔽。

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

---

## 8. 校验（`tools/validate_xindb.py`）

只读、只用公开 schema，可对"别处建好的库"复核。分组：

| 组 | 判据 |
|---|---|
| `S` schema | 3 张必需表 + 4 个索引存在 |
| `M` meta | 必需键齐全；`peak_convention` 等于代码里的常量（不是副本）；`library_kind`/`data_date` 与文件名一致 |
| `C` chunk 记账 | `Σpeak_chunks.n` = `Σentries.n_peaks` = `entry_count`；块数 = `ceil(N/4096)`；`peak_off` 阶梯连续 |
| `E` entry 合理性 | 晶胞/体积物理量合法；每相至少一个元素；`peak_n == n_peaks`；`n_peaks ≤ 150` |
| `P` 反射载荷 | 抽样块解码形状正确；反射有限且为正；强度在 norm-1000 尺度；`n_peaks` 之后是零填充；`top` 块是全反射的子集 |
| `X` 元素倒排 | `eids` 值域合法、严格升序去重；抽样相能在自己 `elements` 的每个元素下找到自己 |

```bash
python -m tools.validate_xindb COD-Inorg-20261020.xinDB --json out.json
```
