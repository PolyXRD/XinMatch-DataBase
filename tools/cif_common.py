#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""polyDB 构建链共用的 CIF 侧工具函数（纯计算，不碰数据库）。

历史说明
--------
本模块的前身是 ``build_dbs.py`` —— 一个 657 行的旧版 **3COD 单体 sqlite**
构建器（产出 ``COD-PCOD-TCOD-{full,inorg}.sqlite``，PolyXRD 的
``phases`` / ``cod_entries`` 两种旧槽位格式）。该格式已于 2026-10-23 彻底
废弃，构建机随之删除；只有下面这几个纯函数被新链路（``scan_inorg.py``
无机分类、``build_polydb.py`` 元数据填充）复用，故单独抽出来。

对外只保证这些符号：

  · ``parse_elements``  —— 分子式 → (元素 CSV, 约简式, 元素计数)
  · ``is_inorganic``    —— 无机判定口径（与历史清单完全一致，不可改）
  · ``cell_volume``     —— 由六个晶胞参数算体积
  · ``list_dirs``       —— 按 ``depth`` 层展开目录（COD 的 3 级 / PCOD 的 2 级布局）
  · ``formula_mass``    —— 元素计数 → 摩尔质量 g/mol（用于算密度）

口径警告
--------
``is_inorganic`` 是**判定口径**，直接决定 ``inorg_cif`` 副本与
``inorg_manifest.sqlite`` 的 1,188,103 条无机记录，任何改动都会让新清单
与已发布的两个库对不上。改之前先想清楚。
"""
from __future__ import annotations

import math
import os
import re

__all__ = [
    "ISOTOPE_TO_H",
    "ATOMIC_MASS",
    "parse_elements",
    "is_inorganic",
    "cell_volume",
    "list_dirs",
    "formula_mass",
]

#: 氘/氚按氢计（CIF 里偶尔出现 D2O 之类写法，不归一到 H 会被判成未知元素）。
ISOTOPE_TO_H = {"D": "H", "T": "H"}

_RE_ELEM = re.compile(r"([A-Z][a-z]?)\s*([0-9]*\.?[0-9]*)")

#: IUPAC 2021 标准原子量（常规值，g/mol）。只覆盖晶态无机物里会出现的元素；
#: 超铀元素与 Es..Og 不在 COD 中出现，缺失时 ``formula_mass`` 返回 None。
#: 精度取到小数点后 3 位 —— 密度只需要 3 位有效数字，多余位数没有意义。
ATOMIC_MASS: dict[str, float] = {
    "H": 1.008, "He": 4.0026, "Li": 6.94, "Be": 9.0122, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998, "Ne": 20.180,
    "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.922, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.468, "Sr": 87.62, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.91,
    "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "I": 126.90, "Xe": 131.29, "Cs": 132.91,
    "Ba": 137.33, "La": 138.91, "Ce": 140.12, "Pr": 140.91, "Nd": 144.24,
    "Pm": 145.0, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05,
    "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21,
    "Os": 190.23, "Ir": 192.22, "Pt": 195.08, "Au": 196.97, "Hg": 200.59,
    "Tl": 204.38, "Pb": 207.2, "Bi": 208.98, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0, "Th": 232.04,
    "Pa": 231.04, "U": 238.03, "Np": 237.0, "Pu": 244.0, "Am": 243.0,
    "Cm": 247.0, "Bk": 247.0, "Cf": 251.0, "Es": 252.0, "Fm": 257.0,
    "Md": 258.0, "No": 259.0, "Lr": 266.0,
}


def parse_elements(formula: str):
    """``'C5 H17 Al N2 O8 P2'`` → ``("Al,C,H,N,O,P", "C5H17AlN2O8P2", {...})``。

    返回 ``(elements_csv, formula_red, counts)``；元素按符号字典序排列，
    约简式按同一顺序拼接（整数计数去掉小数点，非整数保留 ``%g``）。
    空输入返回 ``("", "", {})``。
    """
    counts: dict[str, float] = {}
    for sym, cnt in _RE_ELEM.findall(formula or ""):
        if not sym:
            continue
        sym = ISOTOPE_TO_H.get(sym, sym)
        try:
            counts[sym] = counts.get(sym, 0.0) + (float(cnt) if cnt else 1.0)
        except ValueError:
            counts[sym] = counts.get(sym, 0.0) + 1.0
    if not counts:
        return "", "", {}
    keys = sorted(counts)
    red = "".join(
        k + (str(int(round(counts[k]))) if counts[k] != 1.0
             and abs(counts[k] - round(counts[k])) < 0.05
             else ("" if counts[k] == 1 else "%g" % counts[k]))
        for k in keys)
    return ",".join(keys), red, counts


def is_inorganic(counts: dict) -> bool:
    """无机判定：不含 C，或含 C 但不含 H（碳酸盐/碳化物/CO2 等归无机）。"""
    if not counts:
        return False
    if "C" not in counts:
        return True
    return "H" not in counts


def cell_volume(a, b, c, al, be, ga) -> float:
    """六参数 → 晶胞体积（Å³）。任一轴缺失返回 0.0。"""
    if not (a and b and c):
        return 0.0
    ca = math.cos(math.radians(al))
    cb = math.cos(math.radians(be))
    cg = math.cos(math.radians(ga))
    v = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg
    return a * b * c * math.sqrt(max(1e-12, v))


def list_dirs(root: str, depth: int):
    """按 ``depth`` 层展开子目录，返回 ``[(相对路径, 绝对路径), ...]``。

    COD 是 3 级（``1/00/00``）、PCOD 是 2 级、TCOD 是 3 级 —— 布局不同，
    所以层级必须由调用方按库给定。目录不可读时静默跳过该分支。
    """
    level = [("", root)]
    for _ in range(depth):
        nxt = []
        for rel, p in level:
            try:
                names = os.listdir(p)
            except OSError:
                continue
            for nm in names:
                q = os.path.join(p, nm)
                if os.path.isdir(q):
                    nxt.append((os.path.join(rel, nm) if rel else nm, q))
        level = nxt
        if not level:
            break
    return level


def formula_mass(counts: dict):
    """元素计数 → 摩尔质量（g/mol）。有任一元素缺标准原子量就返回 ``None``。

    返回 None 而不是"跳过该元素继续加"是刻意的：少算一个元素的密度是**错的
    数字**，而 None 在库里是明确缺失，调用方（``build_polydb``）会写成
    ``entries.density IS NULL``，用户一眼能看出没算，不会被假数据误导。
    """
    if not counts:
        return None
    total = 0.0
    for sym, n in counts.items():
        m = ATOMIC_MASS.get(sym)
        if m is None:
            return None
        total += m * n
    return total
