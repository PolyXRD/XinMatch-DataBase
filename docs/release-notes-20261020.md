XinMatch 参考物相库（`.xinDB` v3）。三个库**互不包含，按需装**（XinMatch 支持多库挂载）：COD 实验库是日常交互检索的主力，PCOD（GRINSP 预测）/ TCOD（理论）用于实验库查不到时兜底。

| 库 | entry_id 段 | 相数 | 反射数 | 块数 | 元素 | 解压后 | zip |
|---|---|---|---|---|---|---|---|
| **COD-Inorg** | `96-100-*` | 78,036 | 9,394,774 | 20 | 96 种 | 118,710,272 B（113.2 MiB） | 87.2 MiB |
| **PCOD-Inorg** | `96-210-*` | 1,061,953 | 159,274,265 | 260 | 19 种 | 1,776,353,280 B（1694.1 MiB） | 1.35 GiB |
| **TCOD-Inorg** | `96-110-*` | 46,846 | 5,874,790 | 12 | 91 种 | 71,352,320 B（68.0 MiB） | 54.1 MiB |
| **合计** | — | **1,186,835** | **174,543,829** | 292 | — | 1.83 GiB | ≈1.52 GiB |

- 命名里的 `20261020` 是**源数据快照日**（`meta.data_date`），不是构建日；构建时刻见 `meta.build_time`（COD 2026-10-21 17:25、TCOD 2026-10-22 23:29、PCOD 2026-10-22 23:46）。
- 三个库全部 `journal_mode=DELETE`，**无 `-wal`/`-shm` 侧车**，可直接放只读介质（NAS 只读共享、光盘镜像）。

## 用法

```bash
unzip COD-Inorg-20261020.xinDB.zip                          # → COD-Inorg-20261020.xinDB（单文件）
python -m tools.validate_xindb COD-Inorg-20261020.xinDB     # 可选：自检 25 项
```

> 文件名必须是 `<library_kind>-<data_date>.xinDB`，否则验收器的 `M3.identity` 会判 FAIL（它拿文件名与 `meta` 对账）。

## 校验

```bash
sha256sum -c SHA256SUMS.txt
```

开源验收器 `tools/validate_xindb.py` 只读、只用公开 schema，六组判据：表/索引齐全（S）、meta 完整且与文件名一致（M）、块记账 `Σpeak_n = Σn_peaks = entry_count` 且 `peak_off` 连续（C）、晶胞物理量合法且 `n_peaks ≤ 150`（E）、抽样块解码 + norm-1000 + 零填充 + `top ⊆ all`（P）、元素倒排值域/升序/自洽（X）。

本次实测：

| 库 | 检查项 | FAIL | WARN |
|---|---|---|---|
| COD-Inorg | 26 | 0 | 1（`E1.lattice`：1 条退化晶胞，如 `96-210-4629`） |
| PCOD-Inorg | 25 | 0 | 0 |
| TCOD-Inorg | 25 | 0 | 0 |

另：三库 `PRAGMA integrity_check` 均为 `ok`，`entries` 行数与 `meta.entry_count` 逐库相符，`element_inv` 覆盖的元素集合与 `entries.elements` 实际元素集合完全一致（缺失 0）。

## 峰口径

```
CuKa1-1.5406-d0.7703-30.0-merge-mult-DW-norm1000
```

λ=1.5406 Å（Cu Kα1）｜d 0.7703–30.0 Å｜等价反射按 2θ 合并（含多重性）｜**含 Debye-Waller**（`B = 8π²U_iso`，逐位点）｜强度归一 1000｜每相最多 150 条反射｜预筛矩阵宽 32（只剪枝、不参与排序）

> ⚠️ 与姊妹仓库 [PolyXRD-DataBase](https://github.com/PolyXRD/PolyXRD-DataBase) 的 polyDB **口径不同，不可互换**：polyDB 是 d ≥ 1.0 Å、**无**温度因子、每相 100 峰。同一份 CIF 在两库里算出的绝对强度体系不同。

## 格式

单文件 SQLite，四张表：`meta`（自描述）、`entries`（每相一行 + CSR 偏移）、`peak_chunks`（每 4096 相一块，`top_d/top_i` 预筛矩阵 + `all_d/all_i` 全反射）、`element_inv`（element → 0 起相序号 `<u4` 数组）。寻址铁律：`entries.rowid` 显式从 0 写入 = 相序号，**不能**用 `entry_id` 反推。字段级规范见 [`docs/xinDB格式规范-v3.md`](https://github.com/PolyXRD/XinMatch-DataBase/blob/main/docs/xinDB格式规范-v3.md)。

## 构建信息

| 库 | 构建时刻 | 相数 | 平均峰/相 |
|---|---|---|---|
| COD-Inorg | 2026-10-21 17:25 | 78,036 | 120.4（min 1 / max 150） |
| TCOD-Inorg | 2026-10-22 23:29 | 46,846 | 125.4（min 3 / max 150） |
| PCOD-Inorg | 2026-10-22 23:46 | 1,061,953 | 150.0（全为上限截断） |

构建器：`python -m tools.build_cod_inorg_xindb --source {COD|PCOD|TCOD} --data-date 20261020 --workers 14`，依赖 numpy + gemmi。

## 已知缺口

- COD 库 1 条退化晶胞（`96-210-4629`，a=b=c=1、体积=1，源 CIF 占位值），127 条无空间群。
- 只发布**无机子集**（不含 C，或含 C 但不含 H）；纯有机样品需另建 Full 库。
- PCOD 只有 19 种元素（O/Si/P/Al/V/S/B/F…），是 GRINSP 预测集的客观构成，不是建库丢数据。

---

数据源为 [Crystallography Open Database](http://www.crystallography.net/cod/)，CIF 内容版权归 COD 及各原始文献作者，使用请遵守 COD 条款并引用 COD 与相应文献。本仓库只分发派生索引（反射位置/强度/元数据），不改写也不替代原始 CIF。构建代码见仓库 [README](https://github.com/PolyXRD/XinMatch-DataBase#readme)。
