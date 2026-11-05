# PDF2-20040317.xinDB 构建与验收（ICDD PDF-2 2004 外部库）

> 源：`C:\pdf2\pdf2.dat`（587,999,120 B，80 字符定宽流，7,349,989 行；mtime 2004-03-17）
> 产物：`D:\Project\XRD\XinMatch-DataBase\data\PDF2-20040317.xinDB`，按 `.xinDB` v3 外部库格式构建；
> `meta.library_kind = PDF2`、`meta.data_date = 20040317`，文件名与 meta 一致（`M3.identity` 通过）。
>
> **这份产物不入 git、不进 Release**：ICDD PDF-2 是商品数据。它落在 `XinMatch-DataBase/data/`
> （该目录已被仓库 `.gitignore` 的 `/data/` 屏蔽），PolyXRD 仓库的 `.gitignore` 也加了
> `**/PDF2*.xinDB` 兜底，`tools/package_release.py` 在默认 glob 上也会主动跳过 `PDF2*`。

## 1. 源数据

| 项 | 值 |
|---|---|
| 文件 | `C:\pdf2\pdf2.dat`，587,999,120 B（560.8 MiB），SHA-256 `e25b521c…930a4dbd` |
| 形态 | **无换行的 80 字符定宽流**；子记录标记固定在第 71 列，`<数据集字母><6 位卡片号><序列字母><子记录类型>` |
| 解析记录数 | **163,835** |
| 入库相数 | **163,834**（唯一未入库：卡片 `531651`，无峰表） |
| 反射总数 | **10,445,852**（每相 1~199 条，均值 63.8） |
| 有晶胞 | 133,935（81.75%，类型 E 的**常规胞**；类型 D 的约化胞不入库） |
| 有空间群号 | 119,197（72.76%） |
| 有矿物名 | 17,744（10.83%） |
| d 范围 | **0.36 ~ 99.2 Å** |
| 辐射 | CuKa / CuKa1 / MoKa / MoKa1 / FeKa1 / DD …（峰表是 **d 值**，与辐射无关，仅作溯源） |

**entry_id 约定**：`PDF2-<6 位卡片号>`（如 `PDF2-010837` = ICDD `01-0837`），沿用 PolyXRD
旧索引的 `ref_id` 形态，同一张卡片在两个体系里是同一标识，迁移不用建映射表。

## 2. 与「CIF → 算峰」库的根本差别

这是本次最需要说清的一点：**PDF-2 不是计算库，是实测卡片库**。

| 维度 | COD / PCOD / TCOD / MP | PDF-2 |
|---|---|---|
| 峰从哪来 | 由 CIF 的原子坐标算出 | 几十年来文献里实测的相对强度 |
| `meta.peak_origin` | `computed`（缺省） | **`observed`** |
| `meta.peak_convention` | `CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000` | **`observed-ICDD-PDF2-2004-norm1000`** |
| `entries.tier` | 1（COD）/ 2（MP，DFT） | **0（实验测定）** |
| 晶胞缺失 | 构建缺陷 | 数据事实（18.25%，老 Hanawalt 卡片本就无晶胞） |
| 峰数 > 150 | 算峰引擎安全阀触发，属异常 | 数据事实（14.9% 超 150，最多 199） |
| 建库耗时 | 小时级（要算结构因子） | **41 秒**（只解析 + 落盘） |

强度仍在 `norm1000` 尺度：ICDD 原表按 I/I0=100 归一，`XinDBWriter` 写入时按每相自身
最大值乘到 1000。**两个体系共用同一把强度尺子**，否则同一套打分权重在两类库上不可比。

## 3. 构建

```powershell
cd D:\Project\XRD\XinMatch-DataBase
python -m tools.build_pdf2_xindb -v            # 全量，默认输出到 XinMatch-DataBase\data
python -m tools.build_pdf2_xindb --limit 3000 --out tmp/pdf2/PDF2-20040317.xinDB   # 切片试产
```

| 指标 | 值 |
|---|---|
| 相数 / 反射数 | **163,834 / 10,445,852**（40 chunk = `ceil(163834/4096)`） |
| 耗时 | **41.1 s**（解析 13 s + 落库 28 s，单进程） |
| 产物 | **172,302,336 B（164.3 MiB）**，`journal_mode=DELETE` |
| SHA-256 | `f7bd3df1e4d09f40c94f335e1322a0800c1d36c25e102a7110bbe577be279938` |
| 体积对比 | 旧 `PDF2_2004.sqlite` 214,454,272 B（行式 + 3 索引）→ **小 19.7%** |

解析规则不是这次猜出来的：`xinmatch/core/db/pdf2_scan.py` 移植自 PolyXRD
`services/pdf2_database.py`，那套列位经 `tests/test_pdf2_space_group.py` 的 85,013 条
逐条核对、并由 `PDF2-2004重建验收报告-20261028.md` 的 51 项验收背书。移植后用 3,000 张
卡片与 PolyXRD 解析器**逐字段比对，零差异**（卡片号 / 化学式 / 主名称 / 矿物名 /
空间群 / 空间群号 / 晶胞 / 峰表 / 报告体积）。

### 3.1 写入策略

- **晶胞缺失写 NULL**，不填哨兵值 —— 填 `(1,1,1,90,90,90)` 会造出 d=1.0 Å 的假峰，
  这正是 2026-10-23 在 COD 上踩过的坑。
- **密度不写**：PDF-2 的 Z 与 Dx 在同一张卡片内可以互相矛盾（刚玉报 Z=2 但 Dx 按 Z=6
  算），算出来的密度只会误导筛选。
- **峰不截断**：`peak_compute.DEFAULT_MAX_PEAKS = 150` 是算峰引擎的安全阀，对实测峰表
  没有意义，199 条原样保留。
- **名称位**：`entries.mineral_name` 优先放 `X6` 矿物名（10.83%），缺失时退回主名称
  （化合物名，约 100%）。该列是 GUI 的显示名，只放矿物名会让 89% 的卡片只能显示化学式；
  `meta.name_policy = "mineral>compound"` 记录了这个取舍。
- **文献**：`X9` 拆成 CODEN → `journal`、年份 → `year`、年份之后 → `author`；卷/页在 v3
  的 `entries` 里没有列位，**未落库**（`meta.reference_policy` 已声明）。`journal` 与
  `year` 分别覆盖 100% / 68.3%。

## 4. 验收

### 4.1 结构校验（`tools.validate_xindb`）

```
python -m tools.validate_xindb "D:\Project\XRD\XinMatch-DataBase\data\PDF2-20040317.xinDB" ^
       --sample-chunks 20 --sample-entries 2000
```

**30 项，0 FAIL / 2 WARN**（详见 §5）。关键项：| 检查 | 结果 |
|---|---|
| `M2.convention` | `observed-ICDD-PDF2-2004-norm1000 (observed)` |
| `M3.identity` | `kind=PDF2 date=20040317`（文件名与 meta 一致） |
| `C1`~`C6` | 163,834 相 / 10,445,852 反射 / 40 chunk / `peak_off` 阶梯连续 |
| `E1.lattice` | 29,899 条无晶胞（18.25%）按数据事实保留 NULL；**有值的全部物理** |
| `E4.max_peaks` | `observed library, max 199`（不套算峰引擎上限） |
| `P1.payload` | 20 chunk / 81,920 相解码干净，norm-1000 成立 |
| `X1`~`X4` | 96 个元素倒排项；抽样 2,000 相的每条元素都能找到自己 |

> 验收器为此新增两条口径：`meta.peak_origin`（`computed` / `observed`，缺省 `computed`，
> 老库判据一字未改）与「实测库的 NULL 晶胞不计缺陷」。三个既存库（COD-Inorg / TCOD-Inorg /
> MP）复跑结论不变：**0 fail / 1 warn**，检查项数一致。

全量解码版（`--all-chunks`，36.6 s，40 个 chunk 逐块解码，共 **31 项**）同样是
**0 FAIL / 2 WARN**：`P1.all_chunks decoded all 40 chunks` + `P1.payload 163,834 entries
decode cleanly, norm-1000 held`。

### 4.2 库内自召回（`tools/_probe_selfrecall.py`）

把相自己的峰表反算成 2θ 再检索回本库 —— d 或 I 在任一侧错位都排不到第一。

| 库 | 抽样 | top-1 | top-5 | 未命中 |
|---|---|---|---|---|
| **PDF2-20040317** | 120 | **119** | 1 | **0** |
| COD-Inorg-20261020（对照） | 120 | 118 | 2 | 0 |

实测库与计算库的自召回表现同一水平，说明「观测峰表 + norm1000」这条路走通了。

### 4.3 ICDD 已知卡片点检

与 PolyXRD 验收报告 §6.1 逐字一致：

| 卡片 | 空间群 | 晶胞 | 最强线 |
|---|---|---|---|
| `010837` 方解石 | `R-3c` / 167 | 4.983 / 4.983 / 17.019，90/90/120 | d=3.040 @ I=100 |
| `461045` 石英 | `P3221` / 154 | 4.913 / 4.913 / 5.405，90/90/120 | d=3.3435 @ I=100（58 峰） |
| `271402` Si | `Fd-3m` / 227 | 5.431³ | d=3.135 @ I=100 |
| `461212` 刚玉 | `R-3c` / 167 | 4.759 / 4.759 / 12.993 | d=2.551 @ I=100 |
| `010646` 石墨 | `P63/mmc` / 194 | 2.47 / 2.47 / 6.8，90/90/120 | d=3.370 @ I=100 |

### 4.4 应用侧挂载

```python
from xinmatch.core.db import multi_db
multi_db.add_external_library(r"D:\Project\XRD\XinMatch-DataBase\data\PDF2-20040317.xinDB",
                              label="ICDD PDF-2 2004 (本地授权)")
```

`probe()` → `ok=True`；`list_libraries()` 显示 `available=True entries=163834`。
**读取侧零改动**：`reader_v3` 只校验 `schema_version == 3`，PDF-2 与 COD 家族走同一条
检索链路，差异全部写在 `meta` 里。

## 5. 已知 WARN（不阻塞）

| 检查 | 数量 | 性质 |
|---|---|---|
| `E2.elements` | 93（0.06%） | 化学式为空或残缺（如 `)x ]n`），元素列表无从生成；这些相在元素过滤检索里不可见 |
| `E5.boundary_peak` | 2 | `PDF2-290085`(d=4.04) / `PDF2-371460`(d=3.95) 只有一条线，是 ICDD 的真实单线卡，不是引擎假峰 |

另有一项**不影响正确性但要知道**的现象：**43,782 相含有 bit ≥ 63 的元素**（Au / Pb / Bi /
Pt 等重元素），`entries.elem_mask` 是单个 SQLite INTEGER，只能装 63 位，所以快路径位掩码
表达不了它们。这些元素的精确信息保留在 `entries.elements` 与 `element_inv` 里
（`meta.element_mask_unmapped` 已列出全部未映射元素），MP 库有同样现象。元素过滤按
「必需元素」语义仍然正确；按「元素池子集」语义需走精确倒排。

## 6. 复现与产物

```powershell
# 构建（41 s）
cd D:\Project\XRD\XinMatch-DataBase
python -m tools.build_pdf2_xindb -v

# 验收（只读）
python -m tools.validate_xindb "D:\Project\XRD\XinMatch-DataBase\data\PDF2-20040317.xinDB"
cd D:\Project\XRD\XinMatch
python -m tools._probe_selfrecall --db "D:/Project/XRD/XinMatch-DataBase/data/PDF2-20040317.xinDB" --sample 120

# 回归测试（19 项，含真实文件的 ICDD 已知卡片）
cd D:\Project\XRD\XinMatch-DataBase
python -m pytest tests/test_pdf2_scan.py -q
```

| 项 | 值 |
|---|---|
| 产物 | `D:\Project\XRD\XinMatch-DataBase\data\PDF2-20040317.xinDB`（164.3 MiB） |
| SHA-256 | `f7bd3df1e4d09f40c94f335e1322a0800c1d36c25e102a7110bbe577be279938` |
| 建库工具 | `tools/build_pdf2_xindb.py` + `xinmatch/core/db/pdf2_scan.py` |
| 旧产物 | `PDF2_2004.sqlite`（205 MiB）已不再被任何代码读取，可随时删除 |
