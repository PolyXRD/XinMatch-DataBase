# tools/ 索引

构库/验收分四步，脚本按四组排（**数据库在本仓库构建，CIF 源在上级目录 `sync-COD`**，
完整约定与发布检查清单见 [`../docs/构建与验收手册.md`](../docs/构建与验收手册.md)）：

```
sync-COD（CIF 源） ──①无机分类/复制──> inorg_cif/{COD,PCOD,TCOD}（Inorg 库才需要）
                                      └──②算反射+打包──> data/*.xinDB ──③结构验收──> PASS
                                                                     └──④来源/载荷/坏相验收 + 专项复核
```

支持 **6 个库**：`COD` / `PCOD` / `TCOD` × `Inorg`（无机子集）/ `Full`（全库）。

## ① 无机子集提取（`scan_inorg.py` 组）

| 脚本 | 作用 |
|---|---|
| `scan_inorg.py` | 扫 `sync-COD/{COD,PCOD,TCOD}`，按"**不含 C，或含 C 但不含 H**"判无机；复制到 `<root>/inorg_cif/{LIB}/`；生成清单 `<root>/inorg_manifest.sqlite`（每文件一行：统一 id、化学式、元素、is_inorg） |
| `cif_common.py` | `parse_elements` / `is_inorganic` / `cell_volume` / `formula_mass`（被 ① 复用；前身 `build_dbs.py` 的纯函数部分） |
| `fastcif.py` | 快速 CIF 解析 + 对称展开（① 判元素时读化学式） |
| `fastxrd.py` | 反射计算（polyDB 侧峰引擎；① 只在 import 链上被带入） |

```bash
python tools/scan_inorg.py --workers 6          # 缺省读写 D:/Project/XRD/sync-COD
XRD_DATA_ROOT=/other/root python tools/scan_inorg.py --workers 6   # 换盘时覆盖
```

`_atomic_scattering_params.json` 放在**仓库根**（`fastxrd.py` 按 `tools/../` 找它）。
① 的判据/id 偏移/目录层级必须与已发布三库逐条一致，**不要重写**。

## ② 构建 `.xinDB`（6 库 + 2 个外部库）

| 脚本 | 作用 |
|---|---|
| `build_cod_xindb.py` | **唯一构建入口**：`--source {COD,PCOD,TCOD}` × `--set {inorg,full}`，自动推导 CIF 根 / `library_kind` / id 偏移；分块列式写入（chunk=4096）、`--resume` 断点续跑、峰值记账 |
| `build_mp_xindb.py` | **Materials Project 外部库**：平铺 CIF，`entry_id=mp-<n>`（取文件头 `# mp-id:`），真实空间群/晶系取自注释，`tier=2`；`--data-date` 缺省 `20261031` |
| `build_pdf2_xindb.py` | **ICDD PDF-2 外部库**：读 `C:\pdf2\pdf2.dat` 的 80 字符定宽流（解析器在 `../xinmatch/core/db/pdf2_scan.py`），**不算峰**——PDF-2 本身就是实测峰表；`entry_id=PDF2-<6 位卡片号>`、`tier=0`、`meta.peak_origin=observed`。约 41 s 建 163,834 相；默认输出到 `PolyXRD\cod_data\`（ICDD 数据不得进发布目录） |
| `build_cod_inorg_xindb.py` | 薄兼容包装（默认 `--set inorg`，保留旧 `--source` 用法与 `meta.build_tool`） |
| `../xinmatch/core/db/build_pipeline.py` | 构建热路径：`compute_phase()` 每个 CIF **只解析一次**，同时产出 `PhaseMeta` 与 `(d, i)`（与旧两步逐位等价） |

```bash
# 六个库（默认输出 <repo>/data/<kind>-<date>.xinDB；--out 可指定）
python -m tools.build_cod_xindb --source COD  --set inorg --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source COD  --set full  --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source PCOD --set inorg --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source PCOD --set full  --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source TCOD --set inorg --data-date 20261020 --workers 14 -v
python -m tools.build_cod_xindb --source TCOD --set full  --data-date 20261020 --workers 14 -v
# Materials Project 外部库（155k CIF，约 9 min）
python -m tools.build_mp_xindb --data-date 20261031 --workers 14 -v
# ICDD PDF-2 外部库（163,834 相实测峰表，约 41 s；产物不进 git/Release）
python -m tools.build_pdf2_xindb -v

# 先试产切片（永远先跑这个）
python -m tools.build_cod_xindb --source COD --set inorg --limit 1000 --workers 4 --out tmp/slice.xinDB
# 中断续跑：原命令加 --resume，不要删库重来
python -m tools.build_cod_xindb --source PCOD --set full --resume -v --workers 14
```

输出优先级：`--out` > `XINMATCH_XINDB_DIR` > `<repo>/data`。
**文件名必须**是 `<library_kind>-<data_date>.xinDB`，否则 ③ 的 `M3.identity` 判 FAIL。
`--kind full|inorg` 仅为 COD 的向后兼容别名（与 `--set` 同给会报错）。

## ③ 结构验收

| 脚本 | 作用 |
|---|---|
| `validate_xindb.py` | 只读结构验收器：`S` 表/索引、`M` meta 与文件名一致、`C` 块记账、`E` entry 合理性（**E5 单峰坏相**、**E6 空相**）、`P` 反射载荷、`X` 元素倒排；`--all-chunks` 可全块解码 |

判据按 `meta.peak_origin` 分两套：缺省 `computed`（峰由 CIF 算出，`peak_convention` 必须
等于引擎常量、NULL 晶胞算缺陷、峰数受 `DEFAULT_MAX_PEAKS` 约束）；`observed`（PDF-2 这类
实测卡片库，`peak_convention` 必须自报 `observed-*`+`norm1000`、NULL 晶胞是数据事实、峰数
上限交给 `meta.peak_max` 自查）。缺省值让老库的判据一字不变。

```bash
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --json tmp/validate.json
python -m tools.validate_xindb data/COD-Inorg-20261020.xinDB --all-chunks   # 发布前全量峰载荷
```

## ④ 通用验收 + 专项复核

| 脚本 | 作用 |
|---|---|
| `accept_xindb.py` | **6 库一条命令**：全量 id 号段、来源归属（源目录 / `inorg_manifest.sqlite` 抽查）、`element_inv` 全量一致、峰载荷（默认全量，可抽样）、坏相穷举；`FAIL` 退出码 1 |
| `_compare_xindb_peaks.py` | **只改元数据没动物理**：逐 `peak_chunks` 行比 BLOB 摘要（按 `entry_id` 关联，容许行序不同） |
| `_verify_resume_equivalence.py` | **`--resume` 等价于一次跑完**：构造残库续跑，逐表比对 |
| `_bench_build.py` | 构建吞吐基准（旧两步 vs 单次解析、workers/chunksize），产物写 `tmp/bench/` |

```bash
python -m tools.accept_xindb data/COD-Inorg-20261020.xinDB \
    --manifest D:/Project/XRD/sync-COD/inorg_manifest.sqlite --json tmp/accept.json
python -m tools._verify_resume_equivalence          # 默认合成数据，秒级
```

## ⑤ 发布打包

| 脚本 | 作用 |
|---|---|
| `package_release.py` | 把 `data/*.xinDB`（或指定文件名）打成 zip（Deflate level 1，~1.2×）并生成 `SHA256SUMS.txt`；配 `gh release create` 使用，步骤见构建手册 §5 |

```bash
python -m tools.package_release --out _release/20261104
```

## 测试

```bash
cd D:/Project/XRD/XinMatch-DataBase
python -m pytest tests -q          # 不依赖真实 1.8 GiB 库；用 tests/conftest.py 的合成 CIF/归档
```
