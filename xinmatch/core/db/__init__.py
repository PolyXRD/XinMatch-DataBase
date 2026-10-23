"""XinMatch 参考库层 —— **仅构库所需子集**。

本仓库（XinMatch-DataBase）只带 ``.xinDB`` 构建链路需要的模块：

``cif_meta``       CIF 元数据提取（化学式/空间群/晶胞/密度/元素）
``cif_scan``       CIF 归档扫描与 id 分配（COD / PCOD / TCOD 无机子集）
``peak_compute``   反射计算与峰口径（λ=1.5406，含 Debye-Waller 温度因子）
``schema_v3``      ``.xinDB`` v3 表结构、必需键与冻结的峰口径常量
``writer_v3``      分块列式写入器（entries + peak_chunks + element_inv）

原始 XinMatch 还有 ``schema`` / ``build_cod`` / ``index_build`` / ``query_peaks`` /
``query_filter`` / ``manager`` / ``user_db`` / ``search_v3`` 等检索与 GUI 侧模块，
**不在本仓库**，见 https://github.com/PolyXRD/XinMatch 。
"""
