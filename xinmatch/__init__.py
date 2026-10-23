"""XinMatch — 粉末 XRD 物相分析平台（**构库子集**）。

本仓库只包含参考库 ``.xinDB`` 的构建与校验链路，不包含 GUI、寻峰、检索打分、
Le Bail / Rietveld 精修、定量等模块 —— 那些在完整版仓库里。

完整程序入口（``load_and_preprocess`` / ``identify_phase`` / ``full_analysis``）
由完整版 ``xinmatch/__init__.py`` 提供；这里的 ``__init__`` 故意保持空壳，
避免构库时把 scipy / GUI 等无关依赖一并拖进来。

完整版：https://github.com/PolyXRD/XinMatch
格式规范：``docs/xinDB格式规范-v3.md``
"""

__all__: list[str] = []
