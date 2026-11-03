# XinMatch 参考物相库 20261020 重建版（2026-11-01/04 构建）

> 数据快照仍是 **`20261020`**（`meta.data_date` 不变）；本 Release 是**修复峰引擎后的重建版**，
> 取代 `20261020` Release 中的三个 Inorg 库。旧的三个 Inorg zip 仍保留，但含已知坏相，**不建议再使用**。

## 本次修复

1. **2θ=180° 边界假峰**：`d=λ/2=0.7703 Å` 的反射因浮点舍入进入峰表、LP 因子发散，归一化后整谱塌缩成 1 条假峰。现已剔除该边界反射并给 LP 三角因子取 floor。
2. **gemmi 占位晶胞假峰**：读不到晶胞时不再用默认 `(1,1,1,90,90,90)` 算 `d=1.0`，改为从 CIF 文本恢复晶胞；确实无晶胞的相在建库时跳过。
3. 构建引擎优化（单 CIF 单次解析、结构因子分块按原子数自适应，输出逐位不变），全库吞吐约 +10–15%。

## 五个库

| 库 | entry_id 段 | 相数 | 反射数 | 元素 | 解压后 | zip | SHA-256 |
|---|---|---:|---:|---:|---:|---:|---|
| **COD-Inorg** | `96-100-*` … `96-901-*` | 78,035 | 9,395,681 | 96 | 118,640,640 B（113.1 MiB） | 92,590,746 B | `3da3bcbd88ff322126ce011bb40f30b3196f4bd817126de8dcc8f6ffb4232a56` |
| **TCOD-Inorg** | `96-110-*` … | 46,846 | 5,874,790 | 91 | 71,352,320 B（68.0 MiB） | 57,326,743 B | `3cf0111500cb814332502aadc1751ccf742920a5f3d08c7c891525c92bf42e53` |
| **PCOD-Inorg** | `96-210-*` … `96-295-*` | 1,061,953 | 159,275,318 | 19 | 1,777,357,376 B（1694.1 MiB） | 1,458,071,868 B | `024e46a2c4298e878661fdd69d566e261540a5992868ef486cbaa606c7910193` |
| **TCOD-Full** | `96-110-*` … | 361,701 | 53,099,238 | 91 | 623,173,632 B（594.3 MiB） | 505,902,529 B | `a14d41bc946d94d20f0563ec3188d0b516f2ab111a8bf5206b9c4e650d9164c8` |
| **COD-Full** | `96-100-*` … `96-901-*` | 532,038 | 77,458,154 | 96 | 922,120,192 B（879.4 MiB） | 738,273,981 B | `b8d0a425ffa1c3ee5035be80ea11e80b13d9fb740e8364788373637f6c842768` |

- 库**互不包含，按需装**：日常交互检索用 COD-Inorg；查不到时用 PCOD/TCOD-Inorg 兜底；Full 库覆盖全部 COD/TCOD 归档（含有机物/理论结构）。
- **PCOD-Full 不提供**：PCOD 归档 1,061,962 条 CIF **全部为无机**（`is_inorg=1` 100%），与 PCOD-Inorg 实质等价（原始归档仅差 2 个未入库 CIF），单独构建无意义。
- 随附 `SHA256SUMS.txt`（5 个 zip 的 SHA-256）。

## 验收摘要（全部 0 FAIL）

| 库 | `validate_xindb --all-chunks` | `accept_xindb` | 与旧库 / Inorg 对比 | 自重召回 |
|---|---|---|---|---|
| COD-Inorg | 0 FAIL / 1 WARN（30 项） | 0 FAIL / 1 WARN（19 项），坏相 0 | 仅 7 条坏相变化 + `96-210-4629`（无晶胞）剔除；78,028/78,035 逐字节相同 | Top-1 20/20、Top-5 20/20 |
| TCOD-Inorg | 0 FAIL / 1 WARN（29 项） | 0 FAIL / 1 WARN（19 项），坏相 0 | 与旧库 **46,846/46,846 逐字节相同** | Top-1 20/20、Top-5 20/20 |
| PCOD-Inorg | 0 FAIL / 1 WARN（29 项） | 0 FAIL / 1 WARN（19 项），坏相 0 | 1,061,936/1,061,953 逐字节相同；17 条边界相变化（7 条登记坏相 + 10 条边界反射剔除/补位） | Top-1 10/10、Top-5 10/10 |
| TCOD-Full | 0 FAIL / 2 WARN（29 项） | 0 FAIL / 1 WARN（19 项），坏相 0 | 与 TCOD-Inorg 公共相 **46,846/46,846 逐字节相同** | Top-1 10/10、Top-5 10/10 |
| COD-Full | 0 FAIL / 2 WARN（30 项） | 0 FAIL / 2 WARN（20 项），坏相 0 | 与 COD-Inorg 公共相 **78,035/78,035 逐字节相同** | Top-1 10/10、Top-5 10/10 |

**坏相清零**：COD 旧库 8 条（7 条恢复真峰、`96-210-4629` 无晶胞剔除，相数 78,036→78,035）、PCOD 旧库 7 条（含 `96-254-30089` 双峰塌缩）全部恢复；三个 Inorg 库 `n_peaks≤1` 均为 **0**。

## 已知 WARN（不阻塞，重建不会消除）

- **半塌缩边界相**：最强反射恰落在 `d≈λ/2`（2θ→180°）但其余峰存在的真实谱 —— COD-Inorg 2,251 / PCOD-Inorg 5,598 / TCOD-Inorg 1,670 / TCOD-Full 3,008 / COD-Full 2,368 条。属峰口径边缘现象（LP 在 180° 发散），不是假峰。
- **COD-Full 2 条普通单峰相**（`96-200-1543` d=12.15 Å、`96-200-2110` d=11.23 Å）：唯一反射是真实低角度线，非边界假峰。
- **TCOD-Full 1 条源晶胞 `alpha=1°`**：源 CIF 异常值原样入库，不影响峰位。
- **元素掩码**：TCOD-Full 48,732 相、COD-Full 81,284 相的元素超出 63 位 `elem_mask`；元素列表完整保存在 `entries.elements` / `element_inv`，按元素检索不受影响（扩宽掩码属 schema v4 决策）。
- COD 126 条 / COD-Full 135 条无空间群（磁结构/超空间群源数据），仅影响按空间群过滤。

## 安装

```bash
unzip COD-Inorg-20261020.xinDB.zip        # 得到 COD-Inorg-20261020.xinDB
```

- 直接替换 `data/` 下同名旧库即可；**文件名必须保持 `<kind>-<date>.xinDB`**，否则 `validate_xindb` 的 `M3.identity` 会报 FAIL。
- 全部 `journal_mode=DELETE`，无 `-wal`/`-shm` 侧车，可直接放 NAS 只读共享/光盘镜像。
- 需要 XinMatch v3.0.0+（`.xinDB` v3 格式）；多库可同时挂载（Inorg + Full 并存）。
- 自检：`python -m tools.validate_xindb COD-Inorg-20261020.xinDB --all-chunks`。

## 校验

```bash
sha256sum -c SHA256SUMS.txt        # Linux/macOS/Git Bash
certutil -hashfile COD-Full-20261020.xinDB.zip SHA256   # Windows
```

## 构建与验收记录

- 构建/验收代码：`XinMatch-DataBase` 仓库 `80a7c25`（六库统一入口 + 单次解析 + 自适应分块 + `E5/E6` 验收器）。
- 验收报告：`XinMatch-DataBase/docs/ACCEPTANCE-20261101-rebuild.md`（含每库命令、坏相清单、逐条目对比与 WARN 说明）。
- 文档提交：`8679bf0`（COD/TCOD-Inorg）、`1bb9c63`（PCOD-Inorg）、`a4b36f0`（TCOD-Full）、`13e85b9`（COD-Full）。
- 数据快照：`sync-COD`（`20261020`），`meta.build_time` 分别记录各库构建时刻。
