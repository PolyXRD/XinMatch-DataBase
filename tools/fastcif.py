#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CIF 快速解析器 —— 面向 165 万条 COD/PCOD/TCOD CIF 的批量建库。

为什么不用 pymatgen CifParser: 实测 97.6 ms/条 (COD 大结构), 116 万条无机
条目光解析就要 31 核心小时。本模块只取建库必需的字段, 用行扫描 + 向量化
对称展开, 目标 < 5 ms/条。

对称展开规则与 pymatgen io.cif 保持一致:
  - 优先取 _symmetry_equiv_pos_as_xyz / _space_group_symop_operation_xyz
  - 缺失则回退 _symmetry_space_group_name_H-M 查 SpaceGroup 表
  - 展开后坐标取小数部分 (i - floor(i)), 按 PBC 去重 (atol = 1e-4,
    与 CifParser.site_tolerance 默认值一致)
"""
from __future__ import annotations

import math
import os
import re

import numpy as np

_SITE_TOL = 1e-4

_RE_NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

_CACHE: dict = {}


def str2float(s: str):
    """'5.430(2)' → 5.430 ;  '?' / '.' → None"""
    if s is None:
        return None
    s = s.strip().strip("'\"")
    s = re.sub(r"\([^)]*\)", "", s)          # 去掉不确定度括号
    if s in ("?", ".", ""):
        return None
    try:
        return float(s)
    except ValueError:
        m = _RE_NUM.search(s)
        return float(m.group()) if m else None


def _field(text: str, *tags) -> str | None:
    """取第一个命中的 _tag 的值（单行键值形式）。"""
    for tag in tags:
        m = re.search(r"(?m)^\s*" + re.escape(tag) + r"\s+(.+?)\s*$", text)
        if m:
            v = m.group(1).strip()
            if v.startswith(";") or v.startswith("'"):
                v = v.strip(";'\"").strip()
            return v
    return None


def _loops(text: str):
    """抽取所有 loop_ 块: 返回 [(tags, rows)]，rows 为按空白切分的列表。

    对含空格的对称操作字符串做保留（引号内或整行）。
    """
    lines = text.splitlines()
    out = []
    i = 0
    n = len(lines)
    while i < n:
        if lines[i].strip().startswith("loop_"):
            tags = []
            j = i + 1
            while j < n and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            if not tags:
                i += 1
                continue
            rows = []
            while j < n:
                s = lines[j].strip()
                if not s or s.startswith("_") or s.startswith("loop_") or s.startswith("#"):
                    break
                if s.startswith(";"):        # 分号多行文本块（罕见，跳过）
                    j += 1
                    while j < n and not lines[j].strip().startswith(";"):
                        j += 1
                    j += 1
                    continue
                if "'" in s or '"' in s:
                    m = re.match(r"['\"](.+?)['\"]", s)
                    rows.append([m.group(1)] if m else s.split())
                else:
                    rows.append(s.split())
                j += 1
            out.append((tags, rows))
            i = j
            continue
        i += 1
    return out


_SYM_TAGS = ("_symmetry_equiv_pos_as_xyz", "_symmetry_equiv_pos_as_xyz_",
             "_space_group_symop_operation_xyz", "_space_group_symop_operation_xyz_")


def _symop_strings(text: str) -> list[str]:
    """取对称操作字符串列表。

    必须单独扫描而不能复用 _loops() 的 split(): 对称操作含空格 (如 'x, y, z'),
    按空白切分会把它拆成 3 段而彻底失效。行可能形如:
        'x, y, z'          (带引号, 单列)
        1 'x, y, z'        (带 site_id + 引号)
        2 -x,y+1/2,-z+1/2  (带 site_id, 无空格无引号)
    """
    lines = text.splitlines()
    i, n = 0, len(lines)
    while i < n:
        if lines[i].strip().startswith("loop_"):
            tags = []
            j = i + 1
            while j < n and lines[j].strip().startswith("_"):
                tags.append(lines[j].strip())
                j += 1
            hit = [k for k, t in enumerate(tags) if t in _SYM_TAGS]
            if hit:
                vals = []
                while j < n:
                    s = lines[j].strip()
                    if not s or s.startswith(("_", "loop_", "#")):
                        break
                    m = re.search(r"['\"]([^'\"]+)['\"]", s)
                    if m:
                        vals.append(m.group(1))
                    else:
                        vals.append(re.sub(r"^\d+\s+", "", s))
                    j += 1
                return [v for v in vals if v]
            i = j
            continue
        i += 1
    return []


def _atom_site_rows(text: str):
    """返回 (col_index_map, rows)，只处理含 fract_x 的 atom_site loop。"""
    for tags, rows in _loops(text):
        if not any(t.startswith("_atom_site_") for t in tags):
            continue
        if not any(t == "_atom_site_fract_x" for t in tags):
            continue
        cmap = {t: k for k, t in enumerate(tags)}
        return cmap, rows
    return None, []


def _elem_of(s: str) -> str:
    """'Fe3+' / 'Fe' / 'O2-' → 'Fe' / 'O'"""
    s = (s or "").strip().strip("'\"")
    m = re.match(r"([A-Z][a-z]?)", s)
    return m.group(1) if m else s


class CifData:
    __slots__ = ("a", "b", "c", "alpha", "beta", "gamma", "volume", "z",
                 "sg_name", "sg_number", "formula", "name",
                 "year", "author", "journal", "crystal_system",
                 "symops", "sites", "ok", "err")

    def __init__(self):
        self.ok = False
        self.err = ""


def _unquote(s):
    """去掉 CIF 文本值最外层的成对引号（PCOD 的 'formula' 常写成 "Si O2"）。"""
    if not s:
        return s
    s = s.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "\"'":
        s = s[1:-1].strip()
    return s


def _first_of_loop(text: str, tag: str):
    """取 loop_ 里某个 tag 的第一个值（如 _publ_author_name 的第一作者）。"""
    m = re.search(r"(?m)^\s*" + re.escape(tag) + r"\s*$", text)
    if not m:
        return None
    tail = text[m.end():]
    for line in tail.splitlines():
        s = line.strip()
        if not s or s.startswith("_") or s.startswith("loop_") or s.startswith("#"):
            continue
        return _unquote(s)
    return None


_CRYSTAL_BY_SG = (
    (2, "triclinic"), (15, "monoclinic"), (74, "orthorhombic"),
    (142, "tetragonal"), (167, "trigonal"), (194, "hexagonal"), (230, "cubic"),
)


def crystal_system_of(sg_number: int, setting: str = "", cell=None) -> str | None:
    """晶系：① CIF 自带 _symmetry_cell_setting ② 空间群号 ③ 由晶胞度量推断。

    ③ 是给 PCOD 用的：PCOD 的 CIF 只有 ``_symmetry_space_group_name_H-M "P1"``，
    既无 IT 号也无 cell_setting，若不按晶胞推断，106 万条 PCOD 的 crystal_system
    会全为 NULL（上一版就是这样）。
    """
    s = (setting or "").strip().lower()
    if s in ("triclinic", "monoclinic", "orthorhombic", "tetragonal",
             "trigonal", "hexagonal", "cubic", "rhombohedral"):
        return "trigonal" if s == "rhombohedral" else s
    if sg_number and 1 <= sg_number <= 230:
        for hi, name in _CRYSTAL_BY_SG:
            if sg_number <= hi:
                return name
    if cell:
        return _cs_from_cell(*cell)
    return None


def _cs_from_cell(a, b, c, al, be, ga, tol=0.01):
    """按晶胞度量判断晶系（只用于缺空间群号的条目）。"""
    try:
        a, b, c = float(a), float(b), float(c)
        al, be, ga = float(al), float(be), float(ga)
    except (TypeError, ValueError):
        return None
    if min(a, b, c) <= 0:
        return None

    def eq(x, y, t=tol):
        return abs(x - y) <= t * max(abs(x), abs(y))

    def deg(x):
        return abs(x - 90.0) <= 0.5

    if deg(al) and deg(be) and deg(ga):
        if eq(a, b) and eq(b, c):
            return "cubic"
        if eq(a, b) or eq(b, c) or eq(a, c):
            return "tetragonal"
        return "orthorhombic"
    if deg(al) and deg(be) and abs(ga - 120.0) <= 0.5 and eq(a, b):
        return "hexagonal"
    if eq(a, b) and eq(b, c) and eq(al, be) and eq(be, ga) and not deg(al):
        return "trigonal"
    if (deg(al) + deg(be) + deg(ga)) == 2:
        return "monoclinic"
    return "triclinic"



def parse_cif(text: str) -> CifData:
    d = CifData()
    try:
        a = str2float(_field(text, "_cell_length_a"))
        b = str2float(_field(text, "_cell_length_b"))
        c = str2float(_field(text, "_cell_length_c"))
        al = str2float(_field(text, "_cell_angle_alpha")) or 90.0
        be = str2float(_field(text, "_cell_angle_beta")) or 90.0
        ga = str2float(_field(text, "_cell_angle_gamma")) or 90.0
        d.a, d.b, d.c = a or 0.0, b or 0.0, c or 0.0
        d.alpha, d.beta, d.gamma = al, be, ga
        d.volume = str2float(_field(text, "_cell_volume")) or 0.0
        d.z = int(str2float(_field(text, "_cell_formula_units_Z")) or 0)
        d.sg_name = _unquote(_field(text, "_symmetry_space_group_name_H-M",
                                    "_symmetry_space_group_name_H-M_alt",
                                    "_space_group_name_H-M_alt") or "")
        d.sg_number = int(str2float(_field(text, "_symmetry_Int_Tables_number",
                                           "_space_group_IT_number")) or 0)
        d.formula = _unquote(_field(text, "_chemical_formula_sum",
                                    "_chemical_formula_analytical",
                                    "_chemical_formula_structural") or "")
        d.name = _unquote(_field(text, "_chemical_name_mineral",
                                 "_chemical_name_common",
                                 "_chemical_name_systematic") or "")

        # ---- 文献/分类元数据（XinMatch 的 year/author/journal/crystal_system 筛选要用）----
        d.journal = _unquote(_field(text, "_journal_name_full",
                                    "_journal_name_short",
                                    "_journal_abbrev") or "") or None
        _y = str2float(_field(text, "_journal_year", "_publ_year"))
        d.year = int(_y) if _y and 1000 <= _y <= 2100 else None
        d.author = _first_of_loop(text, "_publ_author_name")
        if d.author and len(d.author) > 120:
            d.author = d.author[:120]
        d.crystal_system = crystal_system_of(
            d.sg_number,
            _unquote(_field(text, "_symmetry_cell_setting",
                            "_space_group_crystal_system") or ""),
            (d.a, d.b, d.c, d.alpha, d.beta, d.gamma))

        syms = _symop_strings(text)
        if not syms:
            # 回退：用空间群符号查表（与 pymatgen 同策略）
            syms = _symops_from_sg(d.sg_name, d.sg_number)
        d.symops = syms or ["x, y, z"]

        cmap, rows = _atom_site_rows(text)
        sites = []
        if cmap is not None:
            ix = cmap.get("_atom_site_fract_x")
            iy = cmap.get("_atom_site_fract_y")
            iz = cmap.get("_atom_site_fract_z")
            ie = cmap.get("_atom_site_type_symbol")
            il = cmap.get("_atom_site_label")
            io = cmap.get("_atom_site_occupancy")
            iu = cmap.get("_atom_site_U_iso_or_equiv")
            ib = cmap.get("_atom_site_B_iso_or_equiv")
            for r in rows:
                if ix is None or len(r) <= max(ix, iy, iz):
                    continue
                x = str2float(r[ix])
                y = str2float(r[iy])
                z = str2float(r[iz])
                if x is None or y is None or z is None:
                    continue
                if ie is not None and ie < len(r):
                    el = _elem_of(r[ie])
                elif il is not None and il < len(r):
                    el = _elem_of(r[il])
                else:
                    continue
                occ = 1.0
                if io is not None and io < len(r):
                    occ = str2float(r[io]) or 1.0
                u = None
                if iu is not None and iu < len(r):
                    u = str2float(r[iu])
                elif ib is not None and ib < len(r):
                    bb = str2float(r[ib])
                    u = bb / (8.0 * math.pi ** 2) if bb is not None else None
                label = r[il] if (il is not None and il < len(r)) else el
                sites.append((label, el, x, y, z, occ, u))
        d.sites = sites
        d.ok = bool(d.a and d.b and d.c and sites)
        if not d.ok:
            d.err = "缺少晶胞或原子位点"
    except Exception as e:                     # noqa: BLE001
        d.ok = False
        d.err = str(e)
    return d


def _symops_from_sg(sg_name: str, sg_number: int) -> list[str]:
    """CIF 未给 symop 时的回退：按空间群符号/编号查表（与 pymatgen 同策略）。

    注意两点:
      · H-M 名可能带 setting 后缀 ('P 4/n m m :2') 与空格, 两种写法都要试
      · SpaceGroup.symm_ops 返回 SymmOp 对象, 必须转成 xyz 字符串
    """
    # 数据源与 pymatgen.io.cif 完全一致: pymatgen/symmetry/symm_ops.json
    # （该表每项为 {"hermann_mauguin": "F m -3 m", "symops": ["x,y,z", ...]}）
    if "symm_ops" not in _CACHE:
        try:
            import json
            import pymatgen.symmetry.groups as _pg   # pymatgen 是 namespace
            # package, 顶层 __file__ 为 None, 必须用子模块定位
            p = os.path.join(os.path.dirname(_pg.__file__), "symm_ops.json")
            with open(p, "r", encoding="utf-8") as fh:
                _CACHE["symm_ops"] = json.load(fh)
        except Exception:                      # noqa: BLE001
            _CACHE["symm_ops"] = []

    def _find(name: str):
        key = re.sub(r"\s+", "", name)
        for _d in _CACHE["symm_ops"]:
            if key == re.sub(r"\s+", "", str(_d.get("hermann_mauguin", ""))):
                return list(_d["symops"])
        return None

    if sg_name:
        s = re.sub(r"\s*:\s*[A-Za-z0-9]{1,2}\s*$", "", sg_name).strip()
        got = _find(s) or _find(sg_name)
        if got:
            return got

    if sg_number:
        try:
            from pymatgen.symmetry.groups import SpaceGroup
            ops = SpaceGroup.from_int_number(int(sg_number)).symm_ops
            out = []
            for op in ops:
                if isinstance(op, str):
                    out.append(op)
                else:
                    fn = getattr(op, "as_xyz_string", None) or getattr(op, "as_xyz_str")
                    out.append(fn())
            return out
        except Exception:                      # noqa: BLE001
            pass
    return []


# ── 对称操作字符串 → 3x4 仿射矩阵（不依赖 pymatgen，避免逐条构造对象） ──

_RE_TERM = re.compile(r"([+-]?)\s*(?:([0-9]+\.?[0-9]*)/([0-9]+)|([0-9]*\.?[0-9]+))?\s*([xyz])?")
_RE_FRAC = re.compile(r"([+-]?)([0-9]+)/([0-9]+)")


def symop_to_matrix(s: str):
    """'x, y+1/2, -z' → (3x3 rot, 3 trans)"""
    rot = np.zeros((3, 3))
    tr = np.zeros(3)
    parts = [p.strip() for p in s.replace("'", "").replace('"', "").split(",")]
    if len(parts) != 3:
        return None
    axes = {"x": 0, "y": 1, "z": 2}
    for row, p in enumerate(parts):
        for m in re.finditer(r"([+-]?)(?:([\d.]+)/(\d+)|([\d.]+))?([xyz])?", p):
            sign, num, den, dec, var = m.groups()
            if not any((sign, num, den, dec, var)):
                continue
            if var:
                coeff = -1.0 if sign == "-" else 1.0
                if dec:
                    coeff *= float(dec)
                elif num and den:
                    coeff *= float(num) / float(den)
                rot[row, axes[var]] += coeff
            else:
                val = 0.0
                if dec:
                    val = float(dec)
                elif num and den:
                    val = float(num) / float(den)
                if sign == "-":
                    val = -val
                tr[row] += val
    return rot, tr


def group_asymmetric_unit(sites, atol: float = _SITE_TOL):
    """把非对称单元里同一位置的原子行合并成混占位点。

    典型场景 (COD 1530381): Al1 (0, .1849, .3696) occ=0.512 与
    Si1 (0, .1849, .3696) occ=0.488 —— 同一晶体学位置被两种元素部分占据,
    必须合成一个 {Al:0.512, Si:0.488} 的混占位点, 否则展开去重时后一种
    元素会被整体判为重复坐标而丢失 (实测会丢掉全部 Si)。
    """
    groups: list[dict] = []
    for _label, el, x, y, z, occ, _u in sites:
        c = np.array([x, y, z], dtype=float) % 1.0
        hit = None
        if groups:
            arr = np.asarray([g["coord"] for g in groups])
            diff = arr - c
            diff -= np.round(diff)
            m = np.all(np.abs(diff) < atol, axis=1)
            if m.any():
                hit = groups[int(np.argmax(m))]
        if hit is None:
            hit = {"coord": c, "species": {}}
            groups.append(hit)
        hit["species"][el] = hit["species"].get(el, 0.0) + float(occ)
    return groups


def expand_sites(sites, symops, atol: float = _SITE_TOL):
    """对称展开 + PBC 去重（返回每个位点的物种 dict，支持混占）。

    Args:
        sites: [(label, el, x, y, z, occ, u_iso), ...] 非对称单元
        symops: 对称操作字符串列表
    Returns:
        (coords(N,3), species_list, uiso(N,))
        species_list[i] = {元素符号: 占据数}
    """
    mats = []
    for s in symops:
        m = symop_to_matrix(s)
        if m is not None:
            mats.append(m)
    if not mats:
        mats.append((np.eye(3), np.zeros(3)))

    groups = group_asymmetric_unit(sites, atol)

    out_c: list[np.ndarray] = []
    out_sp: list[dict] = []
    out_u: list[float] = []

    for g in groups:
        base = g["coord"]
        cand = np.empty((len(mats), 3), dtype=float)
        for k, (rot, tr) in enumerate(mats):
            cand[k] = rot @ base + tr
        cand -= np.floor(cand)
        cand[cand >= 1.0] -= 1.0            # 消除 1.0000000001 之类
        if out_c:
            arr = np.asarray(out_c)
            diff = cand[:, None, :] - arr[None, :, :]
            diff -= np.round(diff)
            match = np.all(np.abs(diff) < atol, axis=2)
            dup_idx = np.argmax(match, axis=1) if match.any() else None
            dup_mask = match.any(axis=1)
        else:
            dup_mask = np.zeros(len(mats), dtype=bool)
            dup_idx = None
        seen: list[np.ndarray] = []
        for k in range(len(mats)):
            c = cand[k]
            if seen:
                dd = np.asarray(seen) - c
                dd -= np.round(dd)
                if np.all(np.abs(dd) < atol, axis=1).any():
                    continue        # 同一非对称原子在本组内自重复
            seen.append(c)
            if dup_mask[k] and dup_idx is not None:
                target = out_sp[int(dup_idx[k])]
                for el, oc in g["species"].items():
                    target[el] = target.get(el, 0.0) + oc
                continue
            out_c.append(c)
            out_sp.append(dict(g["species"]))
            out_u.append(0.0)
    if not out_c:
        return np.zeros((0, 3)), [], []
    return np.asarray(out_c), out_sp, np.asarray(out_u, dtype=float)


def build_structure(d: CifData):
    """构造 pymatgen Structure（含对称展开）；失败返回 None。"""
    try:
        from pymatgen.core import Lattice, Structure
    except Exception:                          # noqa: BLE001
        return None
    if not d.ok:
        return None
    lat = Lattice.from_parameters(d.a, d.b, d.c, d.alpha, d.beta, d.gamma)
    coords, sp_list, uiso = expand_sites(d.sites, d.symops)
    if not sp_list:
        return None
    # 必须带 occupancy: 部分占据位点（溶剂/无序/混占）直接影响结构因子，
    # 全按 1.0 会把 H、C 等轻原子的散射权重算错（实测导致峰强偏差 ~70%）
    species = []
    for sp in sp_list:
        clean = {}
        for el, o in sp.items():
            try:
                o = float(o)
            except (TypeError, ValueError):
                o = 1.0
            if not (o > 0.0):
                o = 1.0          # '?' / 0 / 负值 → 退化为满占据
            elif o > 1.0:
                o = 1.0          # 个别 CIF 的异常值, pymatgen 会直接抛错
            clean[el] = clean.get(el, 0.0) + o
        species.append(clean)
    return Structure(lat, species, coords, coords_are_cartesian=False)
