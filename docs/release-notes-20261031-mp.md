# XinMatch xinDB MP-20261031（Materials Project 外部库）

> 数据源：**Materials Project summary**（`D:\Project\MP-summary`，CIF 下载日 **2026-10-31**）
> 格式：`.xinDB` v3；`meta.library_kind = MP`、`meta.data_date = 20261031`。
> 这是 XinMatch 的**外挂（外部）数据库**：单文件，按需挂载，不随应用打包。

## 库概览

| 项 | 值 |
|---|---|
| CIF 源 | 155,204 个（平铺目录，文件名为化学式） |
| 入库相数 | **155,201**（3 个 CIF 无 `# mp-id:` 头，按约定跳过） |
| 反射数 / 块数 | **20,900,193** / 38 chunk |
| 元素 | 91 种 |
| 空间群 / 晶系 | **225** 个真实空间群（取自 CIF 头注释）/ 7 个晶系，无空值 |
| 解压后 | **251,633,664 B（240.0 MiB）**，`journal_mode=DELETE`，`integrity_check=ok` |
| zip | **197,925,329 B** |
| SHA-256（zip） | `d20c4a61112bcc2efe094cf11929f1020f58642c2e108328b9edb042580cc6bc` |
| SHA-256（.xinDB） | `26D7A9BED83A57C0AC63AF422330210A075C03660B96AD7D4358778F5289EE20` |

## 与 COD/PCOD/TCOD 库的差异

- **entry_id 用 MP 原生 id**：`mp-<number>`（如 `mp-862690`），不套 `96-XXX-XXXX`；
  号段 `[0, 20,000,000)`（当前最大 3,347,425）。`accept_xindb` 的 `BANDS["MP"]` 记录该约定。
- **真实空间群来自 CIF 注释**：pymatgen 导出的结构写成 `P 1`（`_symmetry_space_group_name_H-M 'P 1'`），
  真实群只在首行注释 `# mp-id: ... | spacegroup: P6_3/mmc | ...` 里；构建器用该注释覆盖
  `entries.space_group`，并由 H-M 符号推导 `crystal_system`（225 个符号 / 7 个晶系，全量核对与归档一致）。
- **tier = 2**（DFT 计算结构）；`journal/author = "Materials Project"`、`quality = "DFT"`、`year = 2025`。
- 强度为**计算值**（无实验展宽），用于物相识别；定量需自行标定。

## 构建与验收

```powershell
cd D:\Project\XRD\XinMatch-DataBase
python -m tools.build_mp_xindb --data-date 20261031 --workers 14 -v
```

- 构建：**537.7 s（9.0 min）**，**17,318 相/分钟**（14 workers），`errors=0`。
- `validate_xindb --all-chunks`：**0 FAIL / 1 WARN / 29 项**（38 块全解码；
  `C7` min=8 max=150 mean=134.7；WARN 为 5,660 条半塌缩边界相）。
- `accept_xindb --mp-root`：**0 FAIL / 1 WARN / 19 项**：
  - `R1.id_band`：155,201 个 `mp-<n>` id 全部规范，数值范围 `[2, 3,347,425]`；
  - `A1.ownership`：**id 集合与归档 155,201 个 `mp-id` 完全相等**（`missing=0`、`extra=0`）；
  - `B1.single_peak`：`bad_phase_ids=[]`，`n_peaks≤1` 相 **0** 条；
  - 91 种元素、`element_inv` 全量一致、载荷 norm-1000/零填充/`top⊆all` 通过。
- 自重召回 20 相：**Top-1 20/20、Top-5 20/20**，中位 **0.52 s**。

## 安装 / 挂载

```powershell
unzip MP-20261031.xinDB.zip          # 得到 MP-20261031.xinDB
# 方式一：CLI 挂载外部库
python -m xinmatch.cli.main libraries --add D:\path\to\MP-20261031.xinDB --label "Materials Project"
# 方式二：GUI「设置 → 参考数据库 → 添加外部库」选择该文件
```

- 文件名保持 `MP-20261031.xinDB`（`validate_xindb` 的 `M3.identity` 会核对）；
- 需要 XinMatch v3.0.0+（`.xinDB` v3）；
- 可与 COD/PCOD/TCOD（Inorg/Full）多库同时挂载。

## 校验

```bash
sha256sum -c SHA256SUMS.txt
certutil -hashfile MP-20261031.xinDB.zip SHA256   # Windows
```

## 已知 WARN（不阻塞）

- **5,660 条半塌缩边界相**：最强反射恰在 `d≈λ/2`（2θ→180°）但其余峰仍存在；与 COD 族同类，
  属峰口径边缘现象（真实反射，非假峰），重建不会消除。
- 3 个 CIF（`Al(CuO2)3.cif`、`Ca2YAgO5.cif`、`ZnCoF5.cif`）没有 `# mp-id:` 头，无法分配稳定 id，按约定跳过。

## 代码与验收记录

- 新增 `xinmatch/core/db/mp_scan.py`（MP 归档扫描 / id / 空间群注释 / 晶系）、
  `tools/build_mp_xindb.py`（MP 构建入口）、`tools/accept_xindb.py` 的 MP 号段与全集归属校验、
  `tests/test_mp_scan.py`（回归）。
- 完整验收见仓库 `docs/MP-库构建与验收.md`；发布说明与代码提交在同一仓库历史中。
