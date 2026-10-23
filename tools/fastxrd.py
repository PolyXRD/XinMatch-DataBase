#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""矢量化 X 射线粉末衍射峰计算（与 pymatgen XRDCalculator 等价，但快 1~2 个量级）。

等价性依据：完全照搬 pymatgen 的 physics
  - 原子散射因子 Cromer-Mann 形式
        f(s) = Z - 41.78214 * s^2 * Σ_i a_i exp(-b_i s^2),  s = sinθ/λ = |g|/2
    系数取自 pymatgen 自带的 atomic_scattering_params.json（与官方库同源）
  - 结构因子 F(hkl) = Σ_j f_j · occ_j · exp(2πi g·r_j) · exp(-B_j s^2)
  - 强度 I = |F|^2 · Lorentz-Polarization，(1 + cos²2θ)/(sin²θ cosθ)
  - 相同 2θ 的衍射累加（多重性），容差 1e-5
  - 相对强度低于 0.001（归一化到 100 后）的弱线丢弃
  - 输出按 d 降序（等价于 2θ 升序），归一化到 max = 1000（PolyXRD 无机库约定）

与 pymatgen 的差异（仅实现层面，物理等价）:
  - hkl 用 numpy 网格枚举 + 球壳筛选，取代 get_points_in_sphere 的 Python 循环
  - 结构因子整块矩阵运算，取代逐 hkl 循环
  - 可选只算半个倒易球（Friedel 对 |F|^2 相等，归一化后强度一致，省一半算力）
"""
from __future__ import annotations

import json
import math
import os

import numpy as np

_TWO_THETA_TOL = 1e-5
_SCALED_INTENSITY_TOL = 0.001  # 相对 100 的百分比阈值

# 原子散射因子 Cromer-Mann 系数（与 pymatgen 同一份数据）
_PARAMS_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "_atomic_scattering_params.json",
)

_CACHE: dict = {}


def _load_params() -> dict:
    if "p" in _CACHE:
        return _CACHE["p"]
    src = None
    try:  # 优先用 pymatgen 自带的数据，保证与官方库完全一致
        from pymatgen.analysis.diffraction.xrd import ATOMIC_SCATTERING_PARAMS
        src = ATOMIC_SCATTERING_PARAMS
    except Exception:
        if os.path.exists(_PARAMS_PATH):
            with open(_PARAMS_PATH, "r", encoding="utf-8") as fh:
                src = json.load(fh)
    if src is None:
        raise RuntimeError("无法加载原子散射因子参数表")
    _CACHE["p"] = src
    return src


def reciprocal_matrix(a, b, c, al, be, ga):
    """由晶胞参数算倒易基矢矩阵（行向量 a*, b*, c*，crystallographic 约定不含 2π，
       故 |h a* + k b* + l c*| = 1/d）。"""
    al, be, ga = (math.radians(float(x)) for x in (al, be, ga))
    ca, cb, cg = math.cos(al), math.cos(be), math.cos(ga)
    sg = math.sin(ga)
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * cg, b * sg, 0.0])
    vcz = math.sqrt(max(1e-12, 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg))
    vc = np.array([c * cb, c * (ca - cb * cg) / sg, c * vcz / sg])
    vol = float(np.dot(va, np.cross(vb, vc)))
    if abs(vol) < 1e-12:
        return None
    return np.array([np.cross(vb, vc), np.cross(vc, va), np.cross(va, vb)]) / vol


def _gen_hkl(recip, max_r, min_r, half_sphere, max_hkl, abc):
    """枚举并筛选倒易球内的 hkl。

    逐 h 分块生成（而不是一次性 meshgrid 全空间）: 全空间网格在大胞上可达
    上亿点, 光生成就要好几 GB 内存和数秒时间。分块后峰值内存降到 (2K+1)(2L+1),
    并可按块提前剪枝。
    """
    abc_star = np.linalg.norm(recip, axis=1)
    lims = [int(math.ceil(max_r / max(1e-12, abc_star[i]) * 1.0001)) + 1
            for i in range(3)]
    lims = [max(lims[i], int(math.ceil(max_r * abc[i])) + 1) for i in range(3)]
    lims = [min(l, 400) for l in lims]
    if lims[0] * lims[1] * lims[2] > 4.0e7:      # 网格过大 → 放弃（超大胞）
        return None

    ks = np.arange(-lims[1], lims[1] + 1, dtype=np.float64)
    ls = np.arange(-lims[2], lims[2] + 1, dtype=np.float64)
    K, L = np.meshgrid(ks, ls, indexing="ij")
    kl = np.stack([K.ravel(), L.ravel()], axis=1)          # (M,2)

    # 每个 h 块内 k,l 平面的可行性剪枝阈值
    base = kl @ recip[1:]                                   # (M,3) k*b* + l*c*
    hs = np.arange(-lims[0], lims[0] + 1, dtype=np.float64)

    H_parts, G_parts = [], []
    for h in hs:
        g = h * recip[0] + base
        nrm = np.linalg.norm(g, axis=1)
        keep = (nrm <= max_r) & (nrm >= min_r) & (nrm > 1e-12)
        if not keep.any():
            continue
        H_parts.append(np.column_stack([np.full(keep.sum(), h), kl[keep]]))
        G_parts.append(g[keep])
    if not H_parts:
        return None
    H = np.concatenate(H_parts, axis=0)
    G = np.concatenate(G_parts, axis=0)
    gabs = np.linalg.norm(G, axis=1)

    if half_sphere:  # Friedel 对 |F|^2 相等，归一化后强度一致，省一半算力
        hh, kk, ll = H[:, 0], H[:, 1], H[:, 2]
        m = (hh > 0) | ((hh == 0) & (kk > 0)) | ((hh == 0) & (kk == 0) & (ll > 0))
        H, G, gabs = H[m], G[m], gabs[m]
    if H.shape[0] > max_hkl:        # 按 d 降序截断（保留低角度主区）
        order = np.argsort(gabs, kind="stable")[:max_hkl]
        H, G, gabs = H[order], G[order], gabs[order]
    return H, G, gabs


def compute_pattern_arrays(a, b, c, al, be, ga, frac, symbols, occus,
                           wavelength: float = 1.5406, tth_range=(5.0, 115.0),
                           max_hkl: int = 200_000, chunk: int = 8192,
                           half_sphere: bool = True, scale: float = 1000.0,
                           max_peaks: int = 500, return_hkl: bool = False):
    """纯数组版峰计算（不构造 pymatgen Structure，省掉每原子 PeriodicSite 的开销）。

    Args:
        a,b,c,al,be,ga: 晶胞参数
        frac: (N,3) 分数坐标（已对称展开、已去重）
        symbols: 长度 N 的元素符号列表
        occus: 长度 N 的占据数
        return_hkl: 为 True 时额外返回每条峰的**代表** Miller 指数 (N,3) 列表。
            等价反射按 2θ 合并后，取该组内 2θ 升序的第一个 (h,k,l) 作代表
            （与 pymatgen 对合并线给出的首个指数同策略）。XinMatch 的
            March-Dollase 择优取向拟合要用 peaks.h/k/l，规范 §9 要求回填。
    Returns:
        (d_list, i_list[, hkl_list]) —— d 降序，强度归一化到 scale
    """
    def _empty():
        return ([], [], []) if return_hkl else ([], [])

    recip = reciprocal_matrix(a, b, c, al, be, ga)
    if recip is None:
        return _empty()
    latt_abc = (a, b, c)
    tth_min, tth_max = tth_range
    max_r = 2.0 * math.sin(math.radians(tth_max / 2.0)) / wavelength
    min_r = 2.0 * math.sin(math.radians(tth_min / 2.0)) / wavelength
    if max_r <= 0 or frac is None or len(frac) == 0:
        return _empty()

    got = _gen_hkl(recip, max_r, min_r, half_sphere, max_hkl, latt_abc)
    if got is None:
        return _empty()
    H, G, gabs = got
    n_hkl = H.shape[0]
    if n_hkl == 0:
        return _empty()

    frac = np.asarray(frac, dtype=float)
    occus = np.asarray(occus, dtype=float)
    symbols = list(symbols)

    params = _load_params()
    uniq = sorted(set(symbols))
    try:
        coeff = np.array([params[s] for s in uniq], dtype=float)
    except KeyError as e:
        raise ValueError("缺少散射因子系数: %s" % e)
    zs = np.array([float(_Z[s]) for s in uniq], dtype=float)
    el_of_atom = np.array([uniq.index(s) for s in symbols], dtype=np.int64)

    s_arr = gabs / 2.0
    s2 = s_arr * s_arr
    theta = np.arcsin(np.clip(wavelength * gabs / 2.0, -1.0, 1.0))
    two_theta = np.degrees(2.0 * theta)
    lp = (1.0 + np.cos(2.0 * theta) ** 2) / (np.sin(theta) ** 2 * np.cos(theta))

    intens = np.empty(n_hkl, dtype=float)
    TWOPI = np.float32(2.0 * np.pi)
    # 一次性把常数数组降到 float32: 内层循环的每个临时数组都是
    # (n_hkl × n_atoms) 规模, 全程 f64 会多花 3~4 倍内存带宽。
    # 精度: cos/sin(f32) + f32 成对求和的 |F|² 相对误差 ~1e-5,
    # 比峰表落盘分辨率 (强度 1 位小数 = 1e-4 相对) 还低一个量级。
    H32 = H.astype(np.float32)
    frac32 = np.ascontiguousarray(frac, dtype=np.float32)
    occ32 = occus.astype(np.float32)
    # 散射因子: 全量一次性算 (exp 中间数组只有 n_hkl×n_uniq×4, n_uniq 通常 <10,
    # 很小)。实测 A/B (bench_v3_v4.py, 交替计时消负载漂移): 逐块算 fs 比全量
    # 一次算慢 14% (0.86x), 故保留全量版本。w_all 一次 gather 成连续 f32 大块,
    # 循环内 w_all[sl] 是零拷贝视图; re/im 以 f64 累加, 与 pymatgen 逐点比对
    # max|Δd| = 5.3e-15 Å, max|ΔI| = 0.0033 (1000 归一化)。
    s2c = s2[:, None]
    fs = zs[None, :] - 41.78214 * s2c * np.sum(
        coeff[None, :, :, 0] * np.exp(-coeff[None, :, :, 1] * s2c[:, :, None]),
        axis=2)
    fs32 = fs.astype(np.float32)
    # 每个原子位置的散射因子 × 占据数 —— 每个结构只算一次
    w_all = fs32[:, el_of_atom]
    w_all *= occ32[None, :]
    for start in range(0, n_hkl, chunk):
        sl = slice(start, min(start + chunk, n_hkl))
        ph = H32[sl] @ frac32.T
        ph *= TWOPI
        w = w_all[sl]
        re = (w * np.cos(ph)).sum(axis=1, dtype=np.float64)
        im = (w * np.sin(ph)).sum(axis=1, dtype=np.float64)
        intens[sl] = (re * re + im * im) * lp[sl]

    order = np.argsort(two_theta, kind="stable")
    tt_sorted = two_theta[order]
    it_sorted = intens[order]
    grp = np.empty(len(tt_sorted), dtype=np.int64)
    grp[0] = 0
    np.cumsum(np.diff(tt_sorted) > _TWO_THETA_TOL, out=grp[1:])
    new_grp = np.concatenate(([True], np.diff(tt_sorted) > _TWO_THETA_TOL))
    tt_u = tt_sorted[new_grp]
    it_u = np.zeros(tt_u.shape[0], dtype=float)
    np.add.at(it_u, grp, it_sorted)
    # 合并线的**代表**指数: 每组在 2θ 升序下的第一个 (h,k,l)。
    # new_grp 为 True 的位置即每组首元素(grp 非降, 与 np.unique(return_index) 等价)。
    hkl_u = None
    if return_hkl:
        hkl_u = np.asarray(H)[order][np.flatnonzero(new_grp)]

    imax = float(it_u.max()) if it_u.size else 0.0
    if imax <= 0:
        return _empty()
    rel = it_u / imax * 100.0
    keep = rel > _SCALED_INTENSITY_TOL
    tt_u, it_u = tt_u[keep], it_u[keep]
    if return_hkl:
        hkl_u = hkl_u[keep]
    if tt_u.size == 0:
        return _empty()
    d_u = wavelength / (2.0 * np.sin(np.radians(tt_u / 2.0)))
    scaled = it_u / imax * scale
    if max_peaks and len(d_u) > max_peaks:
        # 超限按**强度**保留，不是按角度截取前 N 条。
        # 按角度截断会把 2θ 高的强线整段丢掉（大胞条目尤其明显），
        # 而物相识别消费的正是 I 最高的那批线；PolyXRD 自带的
        # `regenerate_reference_db_peaks.py` 也是按强度截断后再按 2θ 排序。
        sel = np.argsort(-scaled, kind="stable")[:max_peaks]
        sel.sort()          # 回到 2θ 升序 (= d 降序)
        d_u = d_u[sel]
        scaled = scaled[sel]
        if return_hkl:
            hkl_u = hkl_u[sel]
    if return_hkl:
        return (d_u.tolist(), scaled.tolist(),
                [tuple(int(v) for v in row) for row in hkl_u])
    return d_u.tolist(), scaled.tolist()


def top_peaks_of(d_vals, i_vals, top_n: int = 64):
    """取 I 最高的前 ``top_n`` 个峰，按 I 降序返回 (d, i)。

    与 PolyXRD ``cif_database._top_peaks_of`` 逐元素等价:
    ``order = argsort(-i, kind="stable")[:top_n]``。建库期算好存成 float64
    BLOB，检索时零拷贝直读，省掉每相全峰串解析 + argsort。
    """
    d = np.asarray(d_vals, dtype=np.float64)
    i = np.asarray(i_vals, dtype=np.float64)
    if d.size == 0 or i.size == 0:
        return np.empty(0, dtype=np.float64), np.empty(0, dtype=np.float64)
    order = np.argsort(-i, kind="stable")[:top_n]
    if i.size > d.size:
        order = order[order < d.size]
    return d[order], i[order]


def compute_pattern(structure, wavelength: float = 1.5406,
                    tth_range=(5.0, 115.0), max_hkl: int = 200_000,
                    chunk: int = 8192, half_sphere: bool = True,
                    scale: float = 1000.0, max_peaks: int = 500):
    """计算粉末衍射峰（pymatgen Structure 版，供验证脚本使用）。

    Returns:
        (d_values, intensities) —— d 降序，强度归一化到 scale；无峰返回 ([], [])
    """
    latt = structure.lattice
    p = latt.parameters
    # 注意：混占位点（如 {Al:0.512, Si:0.488}）会展开成多个物种项，
    # 坐标必须**按物种项**逐条复制，否则三个数组长度对不上。
    symbols, occus, rows = [], [], []
    for site in structure:
        for sp, occu in site.species.items():
            sym = getattr(getattr(sp, "element", sp), "symbol", None) or str(sp)
            symbols.append(str(sym))
            occus.append(float(occu))
            rows.append(site.frac_coords)
    frac = np.asarray(rows, dtype=float) if rows else np.zeros((0, 3))
    return compute_pattern_arrays(p[0], p[1], p[2], p[3], p[4], p[5],
                                  frac, symbols, occus,
                                  wavelength=wavelength, tth_range=tth_range,
                                  max_hkl=max_hkl, chunk=chunk,
                                  half_sphere=half_sphere, scale=scale,
                                  max_peaks=max_peaks)
# 元素 → 原子序数（只需覆盖周期表常见元素）
_Z = {}
for _i, _s in enumerate(
        "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
        "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg "
        "Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm".split(), start=1):
    _Z[_s] = _i
_Z["D"] = 1
_Z["T"] = 1
