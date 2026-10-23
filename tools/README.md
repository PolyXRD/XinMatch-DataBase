# tools/ 索引

构库分四步，脚本也按这四组排：

```
原始 CIF 归档 ──①无机分类/复制──> inorg_cif/{COD,PCOD,TCOD}
                                  └──②算反射+打包──> *.xinDB ──③验收──> PASS
                                                                 └──④专项复核（可选）
```

## ① 无机子集提取（`scan_inorg.py` 组）

| 脚本 | 作用 |
|---|---|
| `scan_inorg.py` | 扫 COD/PCOD/TCOD，按"**不含 C，或含 C 但不含 H**"判无机；把无机条目复制到 `<root>/inorg_cif/{LIB}/`；生成清单 `<root>/inorg_manifest.sqlite`（每文件一行：统一 id、化学式、元素、是否无机） |
| `build_dbs.py` | 提供 `parse_elements` / `is_inorganic` / `cell_volume`（被 ① 复用）。**文件里还有旧的 3COD `.sqlite` 打包器 `main()`，本仓库不用，忽略即可** |
| `fastcif.py` | 快速 CIF 解析 + 对称展开（① 判元素时要读化学式） |
| `fastxrd.py` | 反射计算（polyDB 侧的峰引擎；① 只在 import 链上被带入） |

数据根目录由 **`XRD_DATA_ROOT`** 决定（缺省 = 本仓库根）。本仓库不带 200 万条原始 CIF，实际用法：

```bash
XRD_DATA_ROOT=D:/Project/XRD/PolyXRD-DataBase python tools/scan_inorg.py --workers 6
# 读 <root>/COD|PCOD|TCOD，写 <root>/inorg_cif/ 与 <root>/inorg_manifest.sqlite
```

`_atomic_scattering_params.json` 放在**仓库根**（`fastxrd.py` 按 `tools/../` 找它）。

> 这一步沿用 PolyXRD-DataBase 的原实现**原样复制**，没有重写 —— 判据、id 偏移（PCOD +2e7 / TCOD +1e8）
> 与目录层级必须与已发布的三库逐条一致，重写会引入不可比对的偏差。

## ② 构建 `.xinDB`

| 脚本 | 作用 |
|---|---|
| `build_cod_inorg_xindb.py` | **主入口**。`--source {COD,PCOD,TCOD}` 固定源归档与 `library_kind`，套用 id 偏移 |
| `build_cod_xindb.py` | 底层构建器：分块列式写入（chunk=4096）、`--resume` 断点续跑、峰值记账 |

```bash
python -m tools.build_cod_inorg_xindb --source COD --data-date 20261020 --workers 14 -v
python -m tools.build_cod_inorg_xindb --source COD --limit 300 --out /tmp/COD-Inorg-20261020.xinDB  # 先试产
```

输出默认 `$XINMATCH_XINDB_DIR/<kind>-<data_date>.xinDB`，否则 `<repo>/data/`。
**文件名必须**是 `<library_kind>-<data_date>.xinDB`，否则 ③ 的 `M3.identity` 判 FAIL。

## ③ 验收

| 脚本 | 作用 |
|---|---|
| `validate_xindb.py` | 只读验收器，六组 25–26 项：`S` 表/索引、`M` meta 与文件名一致、`C` 块记账、`E` entry 合理性、`P` 反射载荷、`X` 元素倒排 |

```bash
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --json out.json
```

## ④ 专项复核（可选，但发布前建议跑）

| 脚本 | 证明什么 |
|---|---|
| `_compare_xindb_peaks.py` | **只改元数据没动物理**：逐 `peak_chunks` 行比较 `top_d/top_i/all_d/all_i` 的 BLOB 摘要，按 `entry_id` 关联（容许行序不同）。例：`… .xinDB.prev` vs `… .xinDB` |
| `_verify_resume_equivalence.py` | **`--resume` 等价于一次跑完**：先建参考库，再逐个构造"半途被杀"的残局续跑，逐表比对。反射偏移错位不会报错、只会静默串位，所以必须实证 |
