#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""把本地 COD / PCOD / TCOD 三个 CIF 库构建成 PolyXRD 可挂载的 SQLite 数据库。

产物
----
1) COD-PCOD-TCOD-full-20261010.sqlite
   PolyXRD 槽位 cod_index —— 表 cod_entries + cod_atomic_sites + meta
   （db_import.DB_KINDS 要求: cod_entries 必须含 cod_id/formula/elements/cif_gz,
    特征列 cif_gz, 至少 1000 行）

2) COD-PCOD-TCOD-Inorg-20261010.sqlite
   PolyXRD 槽位 cod_inorganics —— 表 phases + cod_atomic_sites + parse_errors + meta
   （要求: phases 必须含 cod_id/formula/n_peaks/peaks_d/peaks_i, 特征列 cell_a,
    且**不能**出现 pearson 列, 否则会被误判成 PDF2 库）

3) inorg_cif/{COD,PCOD,TCOD}/...  无机 CIF 的独立副本（原件保留不动）

ID 唯一化
---------
三个库各自编号会撞车（实测 PCOD 上界 ~1.006e7 与 TCOD 下界 1.0e7 重叠约 6 万条）。
方案: COD / PCOD 沿用原编号; TCOD 原编号 + 100,000,000 (1e8)。
原始编号保留在 file / ref_id / display_id 字段, 可回溯。

用法
----
  python build_dbs.py --plan            # 只打印分片计划
  python build_dbs.py --limit 50        # 小样本试跑（每片最多 50 个 CIF）
  python build_dbs.py                   # 全量构建
  python build_dbs.py --merge-only      # 只做分片合并（worker 已跑完时）
"""
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import shutil
import sqlite3
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np  # noqa: E402
from fastcif import parse_cif, expand_sites  # noqa: E402
from fastxrd import compute_pattern_arrays, top_peaks_of  # noqa: E402

# PolyXRD ``cif_database._PEAK_TOP_N`` —— 预截断强峰列固定 64 个。
# 检索快路径只读这 64 个 float64 (打分最多用到 top-40 / top-12), 有这两列
# 才能走零拷贝直读; 否则每相都要解析全部峰串再 argsort, 在百万级库上很慢。
TOP_N = 64

# 数据根目录：优先环境变量 XRD_DATA_ROOT，否则取本仓库根（tools/ 的上一级）。
# 早期版本硬编码 "E:/rsync-COD"，换机器即失效；build_polydb 只用到本模块的
# parse_elements / formula_mass，并不依赖这些路径，改成可移植默认值不影响建库。
BASE = os.environ.get("XRD_DATA_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(BASE, "_parts")
INORG_CIF_DIR = os.path.join(BASE, "inorg_cif")
FULL_DB = os.path.join(BASE, "COD-PCOD-TCOD-full-20261010.sqlite")
INORG_DB = os.path.join(BASE, "COD-PCOD-TCOD-Inorg-20261010.sqlite")

# 库 → (根目录, 目录层级深度, ID 偏移)
# 实测三库编号互相重叠（COD 与 PCOD 都有 3000034；PCOD 上界与 TCOD 下界也相交），
# 必须错开到互不重叠的号段。原始编号保留在 file/ref_id/display_id 中可回溯。
#   COD  : 原号        ~1.0e6 - 1.0e7
#   PCOD : +2e7        ~2.1e7 - 3.1e7
#   TCOD : +1e8        ~1.1e8 - 1.4e8
LIBS = [
    ("COD", os.path.join(BASE, "COD"), 3, 0),
    ("PCOD", os.path.join(BASE, "PCOD"), 2, 20_000_000),
    ("TCOD", os.path.join(BASE, "TCOD"), 3, 100_000_000),
]

# 超大结构保护：这些在 COD 里是有机大分子/蛋白质，不是无机库的目标，
# 但单条就能吃掉几十秒，会把整库构建拖垮。超限只跳过峰计算（仍进 full 库）。
#: 全量条目数（COD 535,181 + PCOD 1,062,044 + TCOD 57,478 = 1,654,703，
#: os.walk 实测 1,654,588。仅用于进度条估算剩余时间，不影响入库逻辑。
EST_TOTAL_FILES = 1_654_588

MAX_SITES_FOR_PEAKS = 1500      # 对称展开后的位点数
MAX_CELL_VOL_FOR_PEAKS = 30000  # Å³

# 单个结构最多参与结构因子计算的衍射矢量数（按 d 降序截断）。
# 最终只保留前 500 条峰（d 降序），高角度 hkl 必然被截断，算它们是纯浪费；
# 大胞（V 上万 Å³）在 2θ≤115° 内可有数十万 hkl，是耗时的主要来源。
MAX_HKL = 30_000

GZ_LEVEL = 6
WL = 1.5406
TTH = (5.0, 115.0)
MAX_PEAKS = 500
SCALE = 1000.0

# 非金属元素（其余视为金属）——用于无机判定
NONMETALS = {
    "H", "He", "B", "C", "N", "O", "F", "Ne", "Si", "P", "S", "Cl", "Ar",
    "Ge", "As", "Se", "Br", "Kr", "Sb", "Te", "I", "Xe", "At", "Rn",
}
ISOTOPE_TO_H = {"D": "H", "T": "H"}

_RE_ELEM = __import__("re").compile(r"([A-Z][a-z]?)\s*([0-9]*\.?[0-9]*)")

FULL_SCHEMA = """
PRAGMA journal_mode=OFF;
PRAGMA synchronous=OFF;
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS cod_entries (
    cod_id INTEGER PRIMARY KEY, file TEXT NOT NULL, mineral_name TEXT,
    formula TEXT, formula_red TEXT, elements TEXT, space_group TEXT,
    space_group_number INTEGER DEFAULT 0,
    a REAL DEFAULT 0, b REAL DEFAULT 0, c REAL DEFAULT 0,
    alpha REAL DEFAULT 90, beta REAL DEFAULT 90, gamma REAL DEFAULT 90,
    volume REAL DEFAULT 0, z_value INTEGER DEFAULT 0,
    last_updated INTEGER DEFAULT 0, parse_ok INTEGER DEFAULT 1, cif_gz BLOB);
CREATE TABLE IF NOT EXISTS cod_atomic_sites (
    cod_id INTEGER NOT NULL, site_idx INTEGER NOT NULL, label TEXT,
    element TEXT NOT NULL, x REAL NOT NULL, y REAL NOT NULL, z REAL NOT NULL,
    occupancy REAL DEFAULT 1.0, u_iso REAL,
    PRIMARY KEY (cod_id, site_idx));
"""

FULL_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_formula_red ON cod_entries(formula_red);
CREATE INDEX IF NOT EXISTS idx_elements    ON cod_entries(elements);
CREATE INDEX IF NOT EXISTS idx_sg_number   ON cod_entries(space_group_number);
CREATE INDEX IF NOT EXISTS idx_sg_name     ON cod_entries(space_group);
CREATE INDEX IF NOT EXISTS idx_volume      ON cod_entries(volume);
CREATE INDEX IF NOT EXISTS idx_parse_ok    ON cod_entries(parse_ok);
CREATE INDEX IF NOT EXISTS idx_name_fts    ON cod_entries(mineral_name);
CREATE INDEX IF NOT EXISTS idx_cif_file    ON cod_entries(file);
CREATE INDEX IF NOT EXISTS idx_sites_cod_id ON cod_atomic_sites(cod_id);
"""

INORG_SCHEMA = """
PRAGMA journal_mode=OFF;
PRAGMA synchronous=OFF;
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS phases (
    cod_id INTEGER PRIMARY KEY, formula TEXT, space_group TEXT,
    cell_a REAL, cell_b REAL, cell_c REAL,
    cell_alpha REAL, cell_beta REAL, cell_gamma REAL,
    n_peaks INTEGER, peaks_d TEXT, peaks_i TEXT,
    ref_id TEXT, display_id TEXT, cif_gz BLOB,
    peaks_top_d BLOB, peaks_top_i BLOB,
    source TEXT, parse_ok INTEGER DEFAULT 1);
CREATE TABLE IF NOT EXISTS cod_atomic_sites (
    cod_id INTEGER NOT NULL, site_idx INTEGER NOT NULL, label TEXT,
    element TEXT NOT NULL, x REAL NOT NULL, y REAL NOT NULL, z REAL NOT NULL,
    occupancy REAL DEFAULT 1.0, u_iso REAL,
    PRIMARY KEY (cod_id, site_idx));
CREATE TABLE IF NOT EXISTS parse_errors (cod_id INTEGER PRIMARY KEY, reason TEXT);
"""

INORG_INDEXES = """
CREATE INDEX IF NOT EXISTS idx_phases_formula     ON phases(formula);
CREATE INDEX IF NOT EXISTS idx_phases_sg          ON phases(space_group);
CREATE INDEX IF NOT EXISTS idx_phases_display_id  ON phases(display_id);
CREATE INDEX IF NOT EXISTS idx_sites_cod_id       ON cod_atomic_sites(cod_id);
"""


# ──────────────────────────────────────────────────────────────
# 工具函数
# ──────────────────────────────────────────────────────────────

def parse_elements(formula: str):
    """'C5 H17 Al N2 O8 P2' → (elements_csv, formula_red, 元素计数dict)"""
    import re
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
    if not (a and b and c):
        return 0.0
    ca = math.cos(math.radians(al))
    cb = math.cos(math.radians(be))
    cg = math.cos(math.radians(ga))
    v = 1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg
    return a * b * c * math.sqrt(max(1e-12, v))


def list_dirs(root: str, depth: int):
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


# ──────────────────────────────────────────────────────────────
# Worker
# ──────────────────────────────────────────────────────────────

def _shard_stats(path: str):
    """只读探测某分片是否已跑完；跑完则返回上次统计 dict，否则 None。

    任何 SQLite 异常（文件不存在/不是库/上次被强杀导致损坏）都视为未完成，
    交给调用方重写，不会因为一个坏分片把整轮跑挂掉。
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        c = sqlite3.connect("file:%s?mode=ro" % path.replace("\\", "/"),
                            uri=True, timeout=10)
        row = c.execute(
            "SELECT value FROM meta WHERE key='shard_done'").fetchone()
        out = None
        if row and row[0] == "1":
            s = c.execute(
                "SELECT value FROM meta WHERE key='shard_stats'").fetchone()
            if s:
                out = json.loads(s[0])
        c.close()
        return out
    except sqlite3.Error:
        return None


def worker(args):
    shard_id, dir_items, limit, gz_level, copy_inorg, fresh = args
    t0 = time.time()
    fdb = os.path.join(OUT_DIR, "p%03d_full.sqlite" % shard_id)
    idb = os.path.join(OUT_DIR, "p%03d_inorg.sqlite" % shard_id)

    # ── 断点续跑 ──
    # 全量要跑十几小时, 中途被打断（断电/会话结束/异常）不该从头再来:
    # 分片跑完就在 meta 里留 shard_done 标记 + 统计 JSON, 重跑时命中即跳过。
    if not fresh:
        prev = _shard_stats(fdb)
        if prev is not None:
            return prev

    # ── 临时文件 + 原子替换 ──
    # journal_mode=OFF/synchronous=OFF 换来的速度, 代价是进程被强杀时正在写的库会
    # 直接损坏（实测 "database disk image is malformed"）。所以分片一律先写
    # .tmp, 收尾再 os.replace 覆盖正式文件: 中断只留下作废的 tmp, 上一版
    # 正式分片（若已跑完）始终完好, 下次重跑即可安全续上。
    ftmp = fdb + ".tmp"
    itmp = idb + ".tmp"
    for p in (ftmp, itmp):
        if os.path.exists(p):
            try:
                os.remove(p)
            except OSError:
                pass
    fc = sqlite3.connect(ftmp)
    fc.executescript(FULL_SCHEMA)
    ic = sqlite3.connect(itmp)
    ic.executescript(INORG_SCHEMA)

    ent_rows, site_rows = [], []
    ph_rows, ph_sites, ph_err = [], [], []
    made_dirs: set[str] = set()

    n_files = n_inorg = n_peak_ok = n_peak_fail = n_parse_fail = 0
    n_copied = n_copy_fail = n_nodigit = 0

    def flush():
        nonlocal ent_rows, site_rows, ph_rows, ph_sites, ph_err
        if ent_rows:
            fc.executemany(
                "INSERT OR REPLACE INTO cod_entries (cod_id,file,mineral_name,"
                "formula,formula_red,elements,space_group,space_group_number,"
                "a,b,c,alpha,beta,gamma,volume,z_value,last_updated,parse_ok,cif_gz)"
                " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", ent_rows)
            ent_rows = []
        if site_rows:
            fc.executemany(
                "INSERT OR REPLACE INTO cod_atomic_sites "
                "(cod_id,site_idx,label,element,x,y,z,occupancy,u_iso)"
                " VALUES (?,?,?,?,?,?,?,?,?)", site_rows)
            site_rows = []
        if ph_rows:
            ic.executemany(
                "INSERT OR REPLACE INTO phases (cod_id,formula,space_group,"
                "cell_a,cell_b,cell_c,cell_alpha,cell_beta,cell_gamma,"
                "n_peaks,peaks_d,peaks_i,ref_id,display_id,cif_gz,"
                "peaks_top_d,peaks_top_i,source,parse_ok)"
                " VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", ph_rows)
            ph_rows = []
        if ph_sites:
            ic.executemany(
                "INSERT OR REPLACE INTO cod_atomic_sites "
                "(cod_id,site_idx,label,element,x,y,z,occupancy,u_iso)"
                " VALUES (?,?,?,?,?,?,?,?,?)", ph_sites)
            ph_sites = []
        if ph_err:
            ic.executemany("INSERT OR REPLACE INTO parse_errors (cod_id,reason)"
                           " VALUES (?,?)", ph_err)
            ph_err = []

    for lib, depth, offset, rel_dir, abs_dir in dir_items:
        try:
            names = sorted(os.listdir(abs_dir))
        except OSError:
            continue
        cifs = [n for n in names if n.lower().endswith(".cif")]
        if limit:
            cifs = cifs[:limit]
        for nm in cifs:
            path = os.path.join(abs_dir, nm)
            n_files += 1
            try:
                with open(path, "rb") as fh:
                    raw = fh.read()
                mtime = int(os.path.getmtime(path))
            except OSError:
                n_parse_fail += 1
                continue
            text = raw.decode("utf-8", errors="replace")
            stem = os.path.splitext(nm)[0]
            digits = "".join(ch for ch in stem if ch.isdigit())
            orig_id = int(digits) if digits else 0
            if orig_id:
                cod_id = orig_id + offset
            else:
                # 文件名里没有数字（个别库有非数字命名的条目）。若一律给 0，
                # 会在 INTEGER PRIMARY KEY 上互相覆盖 → 静默丢条目。
                # 按 (分片号, 片内序号) 编到 8 亿号段，跨分片天然不冲突。
                n_nodigit += 1
                cod_id = 800_000_000 + shard_id * 1_000_000 + n_nodigit

            d = parse_cif(text)
            elements, formula_red, counts = parse_elements(d.formula)
            parse_ok = int(bool(d.ok and d.a and d.b and d.c))
            if not parse_ok:
                n_parse_fail += 1
            vol = d.volume or cell_volume(d.a, d.b, d.c, d.alpha, d.beta, d.gamma)
            rel = "%s/%s/%s" % (lib, rel_dir, nm) if rel_dir else "%s/%s" % (lib, nm)

            gz = gzip.compress(raw, compresslevel=gz_level)
            ent_rows.append((
                cod_id, rel, d.name[:250] if d.name else "", d.formula,
                formula_red, elements, d.sg_name, d.sg_number,
                d.a, d.b, d.c, d.alpha, d.beta, d.gamma, vol, d.z, mtime,
                parse_ok, gz))

            # 原子位点（非对称单元，与官方库一致）
            for idx, (label, el, x, y, z, occ, u) in enumerate(d.sites):
                site_rows.append((cod_id, idx, label, el, x, y, z,
                                  occ if occ else 1.0, u))
                if len(site_rows) > 10000:
                    flush()

            inorganic = bool(counts) and is_inorganic(counts)
            if inorganic:
                n_inorg += 1
                peaks_d = peaks_i = ""
                n_peaks = 0
                pok = 1
                if parse_ok:
                    try:
                        coords, sp_list, _ = expand_sites(d.sites, d.symops)
                        if len(sp_list) > MAX_SITES_FOR_PEAKS or vol > MAX_CELL_VOL_FOR_PEAKS:
                            n_peak_fail += 1
                            ph_err.append((cod_id, "too_large: sites=%d vol=%.0f"
                                           % (len(sp_list), vol)))
                        elif len(sp_list):
                            # 直接喂数组，不构造 pymatgen Structure：
                            # 上千原子的结构逐个构造 PeriodicSite 要数百毫秒，
                            # 在 116 万条无机条目上是纯浪费。
                            rows, syms, occs = [], [], []
                            for coord, sp in zip(coords, sp_list):
                                for el, o in sp.items():
                                    o = float(o)
                                    rows.append(coord)
                                    syms.append(el)
                                    occs.append(min(o, 1.0) if o > 0 else 1.0)
                            dd, ii = compute_pattern_arrays(
                                d.a, d.b, d.c, d.alpha, d.beta, d.gamma,
                                np.asarray(rows, dtype=float), syms, occs,
                                wavelength=WL, tth_range=TTH,
                                max_peaks=MAX_PEAKS, scale=SCALE,
                                max_hkl=MAX_HKL)
                            if dd:
                                n_peaks = len(dd)
                                peaks_d = ",".join("%.4f" % v for v in dd)
                                peaks_i = ",".join("%.1f" % v for v in ii)
                                n_peak_ok += 1
                            else:
                                n_peak_fail += 1
                        else:
                            n_peak_fail += 1
                    except Exception as e:      # noqa: BLE001
                        n_peak_fail += 1
                        ph_err.append((cod_id, "peak: %s" % str(e)[:180]))
                else:
                    n_peak_fail += 1
                if n_peaks == 0 and pok:
                    pok = 0
                    if not any(r[0] == cod_id for r in ph_err):
                        ph_err.append((cod_id, "no_peaks"))
                # 预截断强峰列: 与 PolyXRD `_top_peaks_of` 同算法, float64 BLOB
                if n_peaks:
                    td, ti = top_peaks_of(dd, ii, TOP_N)
                    top_d, top_i = td.tobytes(), ti.tobytes()
                else:
                    top_d = top_i = None
                ph_rows.append((
                    cod_id, d.formula, d.sg_name,
                    d.a, d.b, d.c, d.alpha, d.beta, d.gamma,
                    n_peaks, peaks_d, peaks_i,
                    "%s-%d" % (lib, orig_id), "%s%d" % (lib, orig_id),
                    gz, top_d, top_i, lib, pok))
                for idx, (label, el, x, y, z, occ, u) in enumerate(d.sites):
                    ph_sites.append((cod_id, idx, label, el, x, y, z,
                                     occ if occ else 1.0, u))
                # 无机 CIF 单独复制一份（原件保留不动）。raw 已在内存里，
                # 直接写出即可，省掉第二遍读盘。
                if copy_inorg:
                    dst_dir = os.path.join(INORG_CIF_DIR, lib, rel_dir)
                    if dst_dir not in made_dirs:
                        os.makedirs(dst_dir, exist_ok=True)
                        made_dirs.add(dst_dir)
                    try:
                        with open(os.path.join(dst_dir, nm), "wb") as fh:
                            fh.write(raw)
                        n_copied += 1
                    except OSError as e:                  # noqa: BLE001
                        n_copy_fail += 1
                        if n_copy_fail <= 3:
                            print("  [复制失败] %s: %s" % (path, e), flush=True)

            # 阈值调低到 500: 每条 CIF 压缩后约 50 KB, 2000 条就是 100 MB,
            # 5 个 worker 合计 500 MB; 本机只有 7.9 GB 内存且已接近用满,
            # 攒批太大会把系统缓存挤掉, 反而拖慢磁盘。
            if len(ent_rows) >= 500:
                flush()

    flush()
    stats = {
        "shard": shard_id, "files": n_files, "inorg": n_inorg,
        "peak_ok": n_peak_ok, "peak_fail": n_peak_fail,
        "parse_fail": n_parse_fail, "copied": n_copied,
        "copy_fail": n_copy_fail, "nodigit": n_nodigit,
        "seconds": round(time.time() - t0, 1),
    }
    # 先落盘再打标记: 只有数据真的写完才允许下次续跑跳过
    fc.execute("INSERT OR REPLACE INTO meta(key,value) VALUES('shard_stats',?)",
               (json.dumps(stats),))
    fc.execute("INSERT OR REPLACE INTO meta(key,value) VALUES('shard_done','1')")
    fc.commit()
    ic.commit()
    fc.close()
    ic.close()
    # 原子替换: 旧文件存在则被覆盖（包括被强杀写坏的）
    os.replace(ftmp, fdb)
    os.replace(itmp, idb)
    return stats


# ──────────────────────────────────────────────────────────────
# 合并
# ──────────────────────────────────────────────────────────────

def merge_parts(parts, target, schema, indexes, tables, label, drop=False):
    if os.path.exists(target):
        os.remove(target)
    conn = sqlite3.connect(target)
    conn.executescript(schema)
    total = {t: 0 for t in tables}
    conflicts: list[tuple] = []
    for i, part in enumerate(parts, 1):
        if not os.path.exists(part):
            continue
        conn.execute("ATTACH DATABASE ? AS p", (part,))
        for t in tables:
            try:
                conn.execute("INSERT INTO main.%s SELECT * FROM p.%s" % (t, t))
                conn.commit()
            except sqlite3.Error as e:
                # 主键冲突（同一 cod_id 出现在两个分片）时必须回滚才能继续,
                # 但绝不能静默丢弃 —— 记下来，最后统计
                conn.rollback()
                dup = conn.execute(
                    "SELECT COUNT(*) FROM (SELECT cod_id FROM p.%s)" % t
                ).fetchone()[0]
                print("  [%s] 分片 %s 的表 %s 合并失败: %s (分片内 %d 行)"
                      % (label, os.path.basename(part), t, e, dup), flush=True)
                conflicts.append((os.path.basename(part), t, str(e)))
        conn.execute("DETACH DATABASE p")
        if drop:
            # 合并完即删: 全量库 ~40 GB, 分片与成品同时在盘会白占一倍空间。
            # 逐个删可以让"成品增长"与"分片收缩"相抵, 峰值占用约等于成品大小。
            try:
                os.remove(part)
            except OSError as e:                          # noqa: BLE001
                print("  [%s] 分片删除失败(忽略): %s: %s"
                      % (label, os.path.basename(part), e), flush=True)
        if i % 20 == 0:
            print("  [%s] 已合并 %d/%d 分片%s"
                  % (label, i, len(parts), " (已释放)" if drop else ""),
                  flush=True)
    for t in tables:
        total[t] = conn.execute("SELECT COUNT(*) FROM %s" % t).fetchone()[0]
    conn.executescript(indexes)
    conn.commit()
    conn.execute("ANALYZE")
    conn.commit()
    conn.close()
    if conflicts:
        print("  [%s] !! 有 %d 次合并冲突: %s" % (label, len(conflicts), conflicts[:5]))
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0, help="每目录最多取 N 个 CIF（试跑）")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--shards", type=int, default=240)
    ap.add_argument("--plan", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    ap.add_argument("--out", default=BASE)
    ap.add_argument("--no-copy-inorg", action="store_true",
                    help="不复制无机 CIF 副本（只建库）")
    ap.add_argument("--keep-parts", action="store_true",
                    help="合并后保留分片文件（默认合并完即删，控制峰值磁盘占用）")
    ap.add_argument("--fresh", action="store_true",
                    help="忽略已有的分片进度，强制重跑（默认断点续跑）")
    a = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)

    # ── 枚举目录单元并分片 ──
    units = []
    for lib, root, depth, offset in LIBS:
        if not os.path.isdir(root):
            print("跳过（目录不存在）: %s" % root)
            continue
        dirs = list_dirs(root, depth)
        if not dirs:
            dirs = list_dirs(root, max(1, depth - 1))
        for rel, _abs in dirs:
            units.append((lib, depth, offset, rel, os.path.join(root, rel)))
    print("目录单元: %d 个" % len(units))

    shards = [[] for _ in range(a.shards)]
    for i, u in enumerate(units):
        shards[i % a.shards].append(u)

    copy_inorg = not a.no_copy_inorg
    tasks = [(i, shards[i], a.limit, GZ_LEVEL, copy_inorg, a.fresh)
             for i in range(a.shards) if shards[i]]
    if a.plan:
        print("分片数: %d, 每片平均 %.0f 个目录" % (len(tasks), len(units) / max(1, len(tasks))))
        return 0

    if not a.merge_only:
        t0 = time.time()
        done = 0
        agg = {"files": 0, "inorg": 0, "peak_ok": 0, "peak_fail": 0,
               "parse_fail": 0, "copied": 0, "copy_fail": 0, "nodigit": 0}
        with Pool(processes=a.workers) as pool:
            for st in pool.imap_unordered(worker, tasks):
                done += 1
                for k in agg:
                    agg[k] += st[k]
                # 每片都报: 全量要跑 ~20 h, 10 片一报太粗, 不便观察速率变化
                if True:
                    el = time.time() - t0
                    rate = agg["files"] / max(1e-9, el)
                    left = (EST_TOTAL_FILES - agg["files"]) / max(1e-9, rate)
                    print("[%d/%d] 文件 %d, 无机 %d, 峰成功 %d, 峰失败 %d, "
                          "解析失败 %d, 已复制 %d, 无数字ID %d, "
                          "用时 %.0f s (%.1f 条目/s, 剩余约 %.1f h)"
                          % (done, len(tasks), agg["files"], agg["inorg"],
                             agg["peak_ok"], agg["peak_fail"],
                             agg["parse_fail"], agg["copied"], agg["nodigit"],
                             el, rate, left / 3600.0), flush=True)
        print("worker 阶段完成: %s" % agg, flush=True)

    full_parts = [os.path.join(OUT_DIR, "p%03d_full.sqlite" % i) for i in range(a.shards)]
    inorg_parts = [os.path.join(OUT_DIR, "p%03d_inorg.sqlite" % i) for i in range(a.shards)]
    full_parts = [p for p in full_parts if os.path.exists(p)]
    inorg_parts = [p for p in inorg_parts if os.path.exists(p)]

    print("合并 full 库 (%d 分片) ..." % len(full_parts), flush=True)
    tf = merge_parts(full_parts, FULL_DB, FULL_SCHEMA, FULL_INDEXES,
                     ["cod_entries", "cod_atomic_sites"], "full",
                     drop=not a.keep_parts)
    print("合并 Inorg 库 (%d 分片) ..." % len(inorg_parts), flush=True)
    ti = merge_parts(inorg_parts, INORG_DB, INORG_SCHEMA, INORG_INDEXES,
                     ["phases", "cod_atomic_sites", "parse_errors"], "inorg",
                     drop=not a.keep_parts)

    for db, meta_extra in ((FULL_DB, "full"), (INORG_DB, "inorg")):
        conn = sqlite3.connect(db)
        conn.executemany("INSERT OR REPLACE INTO meta(key,value) VALUES(?,?)", [
            ("schema_version", "1"),
            ("build_time", time.strftime("%Y-%m-%d %H:%M:%S")),
            ("source_root", BASE),
            ("wavelength", str(WL)),
            ("two_theta_range", "%g-%g" % TTH),
            ("max_peaks", str(MAX_PEAKS)),
            ("peak_trim", "超限按强度保留, 再按 2θ 升序"),
            ("peaks_top_n", str(TOP_N)),
            ("intensity_scale", str(SCALE)),
            ("id_scheme", "COD 原编号; PCOD = 原编号 + 20000000; "
                          "TCOD = 原编号 + 100000000"),
            ("peak_engine", "fastxrd (与 pymatgen XRDCalculator 逐位等价, "
                            "实测 maxΔd≈1e-15 Å, maxΔI=0.00)"),
            ("gzip_level", str(GZ_LEVEL)),
            ("variant", meta_extra),
        ])
        conn.commit()
        conn.close()

    print("完成:")
    print("  full : %s  cod_entries=%d, cod_atomic_sites=%d"
          % (FULL_DB, tf["cod_entries"], tf["cod_atomic_sites"]))
    print("  inorg: %s  phases=%d, cod_atomic_sites=%d, parse_errors=%d"
          % (INORG_DB, ti["phases"], ti["cod_atomic_sites"], ti["parse_errors"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
