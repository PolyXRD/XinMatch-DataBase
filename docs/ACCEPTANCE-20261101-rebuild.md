# 三库重建验收（2026-11-01，峰引擎修复后）

> 背景：2026-10-24 的评估在三个已发布 Inorg 库中确认了两类「假峰相」——
> ① 2θ=180° 边界（`d=λ/2`）反射被浮点舍入放进峰表、LP 发散、归一化后整谱塌缩；
> ② gemmi 读不到晶胞时用默认 `(1,1,1,90,90,90)` 算出 `d=1.0` 假峰。
> 两类根因已在 `peak_compute.py` 修复（见 [`CODE-REVIEW-v3.0.0-20261024.md`](../../XinMatch/docs/CODE-REVIEW-v3.0.0-20261024.md)
> 与 `peak_compute.py` 注释）；本次用修复后的引擎**逐库重建**三个 Inorg 库并验收。
>
> 约定：数据库在本仓库构建，产物写 `data/`；CIF 源为上级 `sync-COD` 的
> `inorg_cif/<LIB>`（`XRD_DATA_ROOT` 可覆盖）；源快照仍为 **`20261020`**（`meta.data_date` 不变）。
> 旧库在重建前已移到 `data/_prev/`，用于逐条目对比；Release 附件（zip）本轮未重打。

**构建环境**：Python 3.10.11 / numpy 2.2.6 / gemmi 0.7.5；构建器 commit `80a7c25`（六库统一入口 + 单次解析 + 自适应结构因子分块）。
**构建命令**（每库单独执行，`XINMATCH_XINDB_DIR` 已清空以确保写本仓库 `data/`）：

```powershell
cd D:\Project\XRD\XinMatch-DataBase
$PY = D:\Project\XRD\python310\python.exe
& $PY -m tools.build_cod_xindb --source COD  --set inorg --data-date 20261020 --workers 8  -v --out data\COD-Inorg-20261020.xinDB
& $PY -m tools.build_cod_xindb --source TCOD --set inorg --data-date 20261020 --workers 14 -v --out data\TCOD-Inorg-20261020.xinDB
& $PY -m tools.build_cod_xindb --source PCOD --set inorg --data-date 20261020 --workers 14 -v --out data\PCOD-Inorg-20261020.xinDB
```

---

## 1. 结果总表

| 库 | 相数（旧→新） | 反射数（旧→新） | 体积（旧→新，字节） | 新 SHA-256 | 构建耗时 | 验收 |
|---|---|---|---|---|---|---|
| **COD-Inorg** | 78,036 → **78,035** | 9,394,774 → **9,395,681** | 118,710,272 → **118,640,640** | `7C08B34D…D7C40E3` | 3,693 s（61.6 min） | **PASS**（validate 0F/1W，accept 0F/1W） |
| **TCOD-Inorg** | 46,846 → **46,846** | 5,874,790 → **5,874,790** | 71,352,320 → **71,352,320** | `8F2A25D5…921C4DC` | 914 s（15.2 min） | **PASS**（validate 0F/1W，accept 0F/1W） |
| **PCOD-Inorg** | 1,061,953 → **1,061,953** | 159,274,265 → **159,275,318** | 1,776,353,280 → **1,777,357,376** | `ABB7EC8C…A1B823` | 47,781 s（13.3 h） | **PASS**（validate 0F/1W，accept 0F/1W） |
| **TCOD-Full**（新库） | — → **361,701** | — → **53,099,238** | 0 → **623,173,632**（594.3 MiB） | `36838FBE…C91238` | 43,469 s（12.1 h） | **PASS**（validate 0F/2W，accept 0F/1W） |
| **COD-Full**（新库） | <!-- CODFULL-ENTRIES --> | <!-- CODFULL-PEAKS --> | <!-- CODFULL-BYTES --> | <!-- CODFULL-SHA --> | <!-- CODFULL-TIME --> | <!-- CODFULL-VERDICT --> |
| **PCOD-Full** | 不构建（与 PCOD-Inorg 等价，见 §5.5） | — | — | — | — | — |

- COD 相数 −1：唯一一条无可用晶胞的源 CIF `96-210-4629` 现在被正确判为「不可计算」并在建库时跳过
  （`skipped_no_peaks` 1,162 → 1,163）。
- COD 反射数 +907：5 条边界坏相恢复完整峰表（122/150/150/150/79 条）+ 2 条占位晶胞坏相恢复（150/114 条）；
  同时不再包含 180° 假反射。
- TCOD 与旧库**逐条目标识**（46,846/46,846 条 d、i 字节完全相同）；SHA 变化仅来自 `meta.build_time` /
  `cif_source`（迁移前路径 → `sync-COD`）/ `build_tool` 等 provenance 字段。

---

## 2. COD-Inorg 验收明细

### 2.1 构建

| 指标 | 值 |
|---|---|
| 源 | `sync-COD/inorg_cif/COD`（79,198 个 CIF） |
| 入库 | **78,035 相 / 9,395,681 反射 / 20 chunk**；`skipped_no_peaks=1,163`，`errors=0` |
| 耗时 / 吞吐 | 3,693 s（61.6 min）/ **1,267.8 相/分钟**（`--workers 8`） |
| 产物 | 118,640,640 B（113.14 MiB），`journal_mode=delete`，`PRAGMA integrity_check=ok` |
| meta | `library_kind=COD-Inorg`、`data_date=20261020`、`build_tool=tools.build_cod_xindb`、`cif_source=…\sync-COD\inorg_cif\COD` |

### 2.2 坏相修复（旧库登记的 8 条）

| entry_id | 旧库峰数 | 新库峰数 | 说明 |
|---|---:|---:|---|
| `96-100-8057` | 1（d=0.7703 假峰） | **122** | 180° 边界塌缩 → 完整粉图恢复 |
| `96-151-1248` | 1 | **150** | 同上（存储上限 150） |
| `96-153-3251` | 1 | **150** | 同上 |
| `96-154-5882` | 1 | **150** | 同上 |
| `96-210-1175` | 1 | **79** | 同上 |
| `96-156-3735` | 1（d=1.0 假峰） | **150** | gemmi 占位晶胞 → 从 CIF 文本恢复晶胞后正常计算 |
| `96-156-3736` | 1 | **114** | 同上 |
| `96-210-4629` | 1 | **不在库中** | 源 CIF 无任何晶胞标签 → 判为不可计算并跳过（预期行为） |

- 全库扫描：`n_peaks <= 1` 的相 **0** 条（旧库 8 条）；「最强反射落在 λ/2 边界窗口且整谱塌缩」的相 **0** 条。

### 2.3 结构/记账验收

```powershell
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --all-chunks --json tmp/rebuild-20261101/COD-validate.json
python -m tools.accept_xindb data/COD-Inorg-20261020.xinDB --manifest D:\Project\XRD\sync-COD\inorg_manifest.sqlite --json tmp/rebuild-20261101/COD-accept.json
```

| 工具 | 结果 | 说明 |
|---|---|---|
| `validate_xindb --all-chunks` | **0 FAIL / 1 WARN / 30 项**（8.7 s） | 全 20 chunk 解码；`E5` 单峰坏相 0；WARN 为 2,251 条「半塌缩」口径边缘相（见 §4，不阻塞） |
| `accept_xindb`（全量载荷 + manifest） | **0 FAIL / 1 WARN / 19 项**（9.3 s） | 全量 id 号段、来源归属 40/40、`element_inv` 全量一致、坏相 0；WARN 同为 2,251 条半塌缩相 |

**前后逐条目对比**（`tmp/rebuild-20261101/compare_entries.py COD`）：

| 项 | 结果 |
|---|---|
| 两库共有条目 | 78,035 |
| `(d, i)` 字节完全相同 | **78,028 / 78,035** |
| 发生变化的条目 | **恰好 7 条**：`96-100-8057`、`96-151-1248`、`96-153-3251`、`96-154-5882`、`96-210-1175`、`96-156-3735`、`96-156-3736` |
| 仅存在于旧库 | 1 条：`96-210-4629`（应剔除） |
| 仅存在于新库 | 0 条 |

> 即：**修复只影响坏相本身，其余 78,028 相的峰数据逐字节未变**。

### 2.4 端到端检索（重建后库）

`tmp/rebuild-20261101/self_recall.py COD 20`（从库内取 20 个随机相反推实验谱，走 `search_match_v3`）：

| Top-1 | Top-5 | 检索中位耗时 |
|---|---|---|
| **20/20** | **20/20** | 0.98 s（测量时 PCOD 正在构建，机器满负载；空闲时应更低） |

---

## 3. TCOD-Inorg 验收明细

### 3.1 构建

| 指标 | 值 |
|---|---|
| 源 / 入库 | `sync-COD/inorg_cif/TCOD`（46,943 CIF）/ **46,846 相 / 5,874,790 反射 / 12 chunk**；`skipped_no_peaks=97`，`errors=0` |
| 耗时 / 吞吐 | 914 s（15.2 min）/ **3,075.7 相/分钟**（`--workers 14`） |
| 产物 | 71,352,320 B（68.05 MiB），`journal_mode=delete`，`integrity_check=ok` |

### 3.2 验收

| 工具 | 结果 |
|---|---|
| `validate_xindb --all-chunks` | **0 FAIL / 1 WARN / 29 项**（全 12 chunk，WARN 为半塌缩口径边缘相） |
| `accept_xindb --manifest` | **0 FAIL / 1 WARN / 19 项** |
| 旧/新逐条目对比 | **46,846 / 46,846 条 `(d, i)` 字节完全相同，0 条变化，0 条增删** —— TCOD 无坏相，重建产物与旧库物理完全一致 |
| 自重召回 20 相 | **Top-1 18/20、Top-5 20/20**，中位 0.282 s（Top-1 未中的 2 例为库内同结构重复相，Top-5 均命中） |

> TCOD 文件 SHA 与旧库不同（`5CD6516E…` → `8F2A25D5…`），但**峰载荷逐字节相同**：
> 差异只来自 `meta.build_time`、`meta.cif_source`（旧路径 → `sync-COD`）和 `meta.build_tool`。

---

## 4. 关于 2,251 / 5,598 / 1,670 条「半塌缩」WARN

`E5` 的 WARN 类指：某相**最强反射**恰落在 `d_min ≈ λ/2`（2θ→180°）窗口、但其余峰仍然存在
（强度和 ≥5/1000）。重建后 COD 库仍有 **2,251** 条（TCOD-Inorg 1,670；PCOD-Inorg 5,598；TCOD-Full 3,008）。

- 这些**不是假峰**：它们是真实反射，只是 2θ 极接近 180°，Lorentz-polarisation 因子在该处发散
  （口径把 `d_min` 定在精确的 λ/2）。引擎修复（剔除 180° 反射、给 LP 三角因子取 floor）**不能也不应**
  删除这些真实反射，所以重建后 WARN 数量与旧库基本一致。
- 若产品希望消除这类极端谱，需要改**峰口径**（例如 `d_min` 抬到 0.78–0.79 Å，或剔除 2θ>179°），
  这会改变 `PEAK_CONVENTION` 并需要重建全部库 —— 属 v4 决策，本次不做；当前仅记录 WARN、不阻塞发布。

---

## 5. PCOD-Inorg 验收明细

### 5.1 构建

| 指标 | 值 |
|---|---|
| 源 / 入库 | `sync-COD/inorg_cif/PCOD`（1,061,962 CIF）/ **1,061,953 相 / 159,275,318 反射 / 260 chunk**；`skipped_no_peaks=9`，`errors=0` |
| 耗时 / 吞吐 | 47,781 s（**13.3 h**）/ **1,333.5 相/分钟**（`--workers 14`；与历史 12.7 h 同量级，期间本机另有外部任务抢占 CPU） |
| 产物 | 1,777,357,376 B（1694.1 MiB），`journal_mode=delete`，`integrity_check=ok` |
| 反射数变化 | +1,053（7 条坏相恢复真峰；另有 10 条相的 180° 边界反射被剔除后由后续峰补位，n_peaks 仍为 150） |

### 5.2 坏相修复（旧库登记的 7 条）

| entry_id | 旧库峰数 | 新库峰数 |
|---|---:|---:|
| `96-231-03437` | 1（d=0.7703 假峰） | **150** |
| `96-232-31453` | 1 | **150** |
| `96-235-07656` | 1 | **150** |
| `96-254-30089` | 2（180° 假峰 + 真峰 I=0.013，塌缩双峰） | **150** |
| `96-254-30092` | 1 | **150** |
| `96-256-30724` | 1 | **150** |
| `96-264-40783` | 1 | **150** |

全库扫描：`n_peaks <= 1` 的相 **0** 条（旧库 6 条）；「最强反射在 λ/2 窗口且整谱塌缩」的相 **0** 条。

### 5.3 结构/记账验收

| 工具 | 结果 |
|---|---|
| `validate_xindb --all-chunks`（260 块全解码，114.7 s） | **0 FAIL / 1 WARN / 29 项**；`C7.peaks_per_phase` min=20 max=150 mean=150.0（旧库 min=1）；WARN 为 5,598 条半塌缩口径边缘相 |
| `accept_xindb --manifest`（全量 260 块 + 全量 id/元素，121 s） | **0 FAIL / 1 WARN / 19 项**；`bad_phase_ids` 空；id 范围 `[1000000, 9585139]`、元素 19 种 |

**前后逐条目对比**（`compare_entries.py PCOD`）：

| 项 | 结果 |
|---|---|
| 共有条目 | 1,061,953 |
| `(d,i)` 字节完全相同 | **1,061,936 / 1,061,953（99.9984%）** |
| 发生变化 | **17 条**（7 条已登记坏相 + `96-244-04910` 等 10 条边界反射被剔除/补位的相） |
| 增删条目 | **0 条**（PCOD 无不可计算相变化） |

### 5.4 端到端检索（重建后库）

`self_recall.py PCOD 10`：**Top-1 10/10、Top-5 10/10**，检索中位 **1.90 s**（106 万相库，加载预筛矩阵后）。

### 5.5 PCOD-Full：不构建（与 PCOD-Inorg 等价）

- `inorg_manifest.sqlite`：PCOD 共 **1,061,962** 个 CIF，**`is_inorg=1` 的比例为 100%**（0 条有机）；
- PCOD 原始归档 **1,061,964** 个 CIF，仅 2 个未进入 ① 清单/无机子集（不足 0.0002%）；
- 因此 `PCOD-Full` 与 `PCOD-Inorg` 在相集合/峰口径上实质等价（差 ≤2 相），按用户决定**不单独构建**。

---

## 6. TCOD-Full 验收明细（新库）

### 6.1 构建

| 指标 | 值 |
|---|---|
| 源 / 入库 | `sync-COD/TCOD`（362,293 个清单 CIF / 362,294 个原始 CIF）/ **361,701 相 / 53,099,238 反射 / 89 chunk**；`skipped_no_peaks=592`，`errors=0` |
| 耗时 / 吞吐 | 43,469 s（**12.1 h**）/ **499.3 相/分钟**（`--workers 14`） |
| 产物 | **623,173,632 B（594.3 MiB）**，`journal_mode=delete`，`integrity_check=ok` |
| SHA-256 | `36838FBE1F18EDBF568DEDC86484FCD58138549258AE2978F658EC3589C91238` |

### 6.2 验收

| 工具 | 结果 |
|---|---|
| `validate_xindb --all-chunks`（89 块全解码，127 s） | **0 FAIL / 2 WARN / 29 项**：`C7` min=3 max=150 mean=146.8；WARN①=3,008 条半塌缩口径边缘相；WARN②=`E1` 1 条源晶胞角 `alpha=1°`（源数据异常，非建库错误） |
| `accept_xindb`（全量载荷，147 s） | **0 FAIL / 1 WARN / 19 项**：`bad_phase_ids=0`、`boundary_phases=0`、`placeholder_phases=0`、`single_peak_phases=0`；id 范围 `[10000001, 40354755]`、元素 **91** 种；WARN 为 3,008 条半塌缩相 |
| 与 `TCOD-Inorg` 公共相逐条目对比（`compare_full_inorg.py TCOD`） | **46,846 / 46,846 条 `(d,i)` 逐字节相同，0 条差异、0 条缺失** —— Full 只"加相"，没有改动 Inorg 已覆盖相的任何物理量 |
| 自重召回 10 相 | **Top-1 10/10、Top-5 10/10**，中位 3.46 s（测量时 COD-Full 正在构建，机器满负载；空闲时更快） |

### 6.3 已知提示（不阻塞）

- 构建日志提示 **48,732 条相的某个元素在 `entries.elem_mask`（63 位）里没有位可放**：
  TCOD-Full 含 91 种元素，超出 63 位掩码容量。这些相的元素列表仍完整保存在 `entries.elements`
  与 `element_inv` 中，按元素检索不受影响；只有「用 63 位掩码做快速预剪枝」的场景会退化。
  若要让掩码覆盖全部元素，需要扩宽 `elem_mask`（schema v4 决策）。
- `E1` 的 `alpha=1°` 是一条源 CIF 的异常晶胞值（原样入库），不影响峰位/检索。

---

## 7. COD-Full

<!-- CODFULL-SECTION -->

---

## 8. 复现与产物

| 文件 | 说明 |
|---|---|
| `data/{COD,TCOD,PCOD}-Inorg-20261020.xinDB` | 重建产物 |
| `data/_prev/*.xinDB` | 重建前旧库（用于逐条目对比；确认后可删除，PCOD 旧库 1.7 GB） |
| `data/_build-20261101-{COD,TCOD,PCOD}.log` | 构建日志（进度行、吞吐、构建 summary） |
| `tmp/rebuild-20261101/inspect_rebuild.py` | 单库身份 + 坏相 + 全库边界扫描 |
| `tmp/rebuild-20261101/compare_entries.py` | 旧/新逐条目 `(d,i)` 字节对比 |
| `tmp/rebuild-20261101/self_recall.py` | 自重召回（真实 `search_match_v3`） |
| `tmp/rebuild-20261101/*-{inspect,validate,accept,selfrecall}.json` | 原始结果 |

```powershell
# 完整复跑（只读验收）
python tmp/rebuild-20261101/inspect_rebuild.py COD
python tmp/rebuild-20261101/compare_entries.py COD
python tmp/rebuild-20261101/self_recall.py COD 20
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --all-chunks
python -m tools.accept_xindb  data/COD-Inorg-20261020.xinDB --manifest D:/Project/XRD/sync-COD/inorg_manifest.sqlite
```
