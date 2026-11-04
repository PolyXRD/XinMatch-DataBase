# MP-20261031.xinDB 构建与验收（Materials Project 外部库）

> 源：`D:\Project\MP-summary\cif`（Materials Project summary，CIF 下载日 **2026-10-31**）
> 产物：`data/MP-20261031.xinDB`，按 `.xinDB` v3 外部库格式构建；`meta.library_kind = MP`、
> `meta.data_date = 20261031`，文件名与 meta 一致（`validate_xindb` 的 `M3.identity` 通过）。

## 1. 源数据

| 项 | 值 |
|---|---|
| CIF 数 | 155,204（平铺目录，文件名 = 化学式） |
| 带 `# mp-id:` 头 | **155,201**（唯一，范围 `mp-2` … `mp-3347425`） |
| 无 mp-id（跳过） | 3：`Al(CuO2)3.cif`、`Ca2YAgO5.cif`、`ZnCoF5.cif` |
| 真实空间群 | 来自 CIF 注释 `# mp-id: mp-... \| spacegroup: ...`（pymatgen 把结构写成 P1，`_symmetry_space_group_name_H-M` 恒为 `P 1`，不能直接用）；共 **225** 个不同符号 |
| 能量 | 注释里的 `energy_above_hull` 仅用于核对；schema 无对应列，未入库 |

**entry_id 约定**：`mp-<number>`（MP 原生 id，如 `mp-862690`），不套用 COD 的 `96-XXX-XXXX`。
号段 `[0, 20,000,000)`（当前最大 3,347,425）；`accept_xindb` 的 `BANDS["MP"]` 记录了该约定。

## 2. 构建

```powershell
cd D:\Project\XRD\XinMatch-DataBase
python -m tools.build_mp_xindb --data-date 20261031 --workers 14 -v
```

| 指标 | 值 |
|---|---|
| 相数 / 反射数 | **155,201 / 20,900,193**（38 chunk） |
| 耗时 / 吞吐 | 537.7 s（**9.0 min**）/ **17,318 相/分钟**（14 workers） |
| 产物 | **251,633,664 B（240.0 MiB）**，`journal_mode=DELETE`，`integrity_check=ok` |
| SHA-256 | `26D7A9BED83A57C0AC63AF422330210A075C03660B96AD7D4358778F5289EE20` |
| meta | `tier_default=2`（DFT 计算结构）、`dataset=MP`、`mp_archive_cifs=155204`、`mp_with_id=155201`、`mp_without_id=3` |

构建管线复用六库的 `build()`：单 CIF 单次解析（`compute_phase`）+ 自适应结构因子分块；
`id_mode="mp"` 只改变 id 生成、任务遍历与 MP 元数据补全，物理计算与 COD 族完全一致。
每个 MP 相额外写入：

- `space_group` / `crystal_system`：取自 CIF 注释的真实空间群（225 个）及由其 H-M 符号推导的晶系（7 种，无空值）；
- `journal = "Materials Project"`、`author = "Materials Project"`、`quality = "DFT"`、`year = 2025`（summary 日期）。

## 3. 验收

### 3.1 结构校验（`validate_xindb --all-chunks`，19 s）

**0 FAIL / 1 WARN / 29 项**。全 38 chunk 解码；`C7.peaks_per_phase` min=8、max=150、mean=134.7；
WARN 为 5,660 条「最强反射落在 d≈0.7703 Å 窗口但谱仍完整」的半塌缩边缘相（真实反射，非假峰）。

### 3.2 产品验收（`accept_xindb --mp-root`，31.8 s）

**0 FAIL / 1 WARN / 19 项**：

| 检查 | 结果 |
|---|---|
| `R1.id_band` | 155,201/155,201 个 id 均为 `mp-<n>` 规范形式，数值范围 `[2, 3,347,425]`，落在 MP 号段 |
| `A1.ownership` | **id 集合与归档完全相等**：`155,201` vs `155,201`，`missing=0`、`extra=0`；3 个无 mp-id 文件按约定跳过 |
| `B1.single_peak` | `bad_phase_ids = []`、`n_peaks≤1` 相 **0** 条；`n_peaks=1` 也无 |
| 元素 | 91 种；`element_inv` 全量一致（含 `elem_mask` 63 位装不下的大元素，如 `Ac/U/W`） |
| 载荷 | 38 chunk 全量解码，norm-1000、零填充、`top⊆all` 通过 |
| WARN | 5,660 条半塌缩边界相（同上，属峰口径边缘现象） |

### 3.3 空间群/晶系一致性（全量复核）

重新扫描 155,201 个 CIF 注释并逐相统计：库内 **225** 个空间群与归档**集合完全一致、每符号计数也完全一致**；
7 个晶系、无 `NULL`。说明没有把 pymatgen 的 `P 1` 误存为真实空间群。

### 3.4 端到端检索（自召回）

`tmp/rebuild-20261101/self_recall.py data/MP-20261031.xinDB 20`：
**Top-1 20/20、Top-5 20/20**，中位 **0.52 s**（15.5 万相）。

### 3.5 应用侧挂载

```powershell
python -m xinmatch.cli.main libraries --add D:\Project\XRD\XinMatch-DataBase\data\MP-20261031.xinDB --label "Materials Project"
```

挂载后 `libraries --json`：`ext_567de5034377`，`available=true`，`entries=155201`。

## 4. 已知 WARN（不阻塞）

- **5,660 条半塌缩边界相**：最强反射恰在 `d≈λ/2`（2θ→180°）、其余峰仍存在；与 COD 族同类，属峰口径边缘现象，不是假峰，重建不会消除。
- **3 个 CIF 无 `# mp-id:` 头**：无法给出稳定 MP id，按约定跳过并在构建日志中列出。
- MP 是 DFT 计算结构（`tier=2`），强度不含实验仪器/样品展宽；用于物相识别可以，定量需自行标定。

## 5. 复现与产物

| 文件 | 说明 |
|---|---|
| `data/MP-20261031.xinDB` | 产物（251,633,664 B，SHA-256 见 §2） |
| `data/_build-20261104-MP.log` | 构建日志 |
| `tmp/release_20261104/MP-{validate,accept}.json` | 校验/验收原始 JSON |
| `xinmatch/core/db/mp_scan.py` | MP 归档扫描（mp-id / 空间群注释 / 晶系） |
| `tools/build_mp_xindb.py` | MP 构建入口（`--data-date 20261031`） |
| `tools/accept_xindb.py` | 新增 `MP` 号段、`mp-` id 解析与集合级归属校验（`--mp-root`） |
| `tests/test_mp_scan.py` | MP 扫描 / id / 晶系 / 构建+验收回归 |
| Release | [`20261031-mp`](https://github.com/PolyXRD/XinMatch-DataBase/releases/tag/20261031-mp)：`MP-20261031.xinDB.zip`（197,925,329 B，SHA-256 `d20c4a61…`）+ `SHA256SUMS.txt` |

```powershell
# 复跑验收（只读）
python -m tools.validate_xindb data/MP-20261031.xinDB --all-chunks
python -m tools.accept_xindb  data/MP-20261031.xinDB --mp-root D:\Project\MP-summary\cif
```
