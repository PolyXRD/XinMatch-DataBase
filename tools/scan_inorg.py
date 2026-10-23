#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""3COD（COD+PCOD+TCOD）无机分类复核 + 副本复制 + 清单生成。

用途（用户要求的第 2、3 步）:
  1. 重新扫描三库，按既有无机判据（不含 C；或含 C 但不含 H）分类；
  2. 无机条目复制到 inorg_cif/{LIB}/{rel_dir}/{name}（已存在且大小一致则跳过，
     故为**增量**复制，可为同步后新增条目补齐）；
  3. 生成 inorg_manifest.sqlite 清单（每文件一行：统一 cod_id、formula、elements、
     是否无机），供打包与"是否漏条目"的集合级核对。

统一编号（三库原始编号互相重叠，必须错位）:
    COD  : 原号            1_000_000 –  9_999_999
    PCOD : +20_000_000    21_000_001 – 29_000_002
    TCOD : +100_000_000  110_000_001 –140_xxxxxx
保留 COD 原号，使 PolyXRD 的 COD REST 在线回退对 COD 条目仍然正确。

并发设计（重要）:
  worker **不写 SQLite**，只往自己的 JSONL 追加（append-only，无锁）。
  此前两版用"每 worker 一个 sqlite 分片"，都在 `con.commit()` 上撞
  "database is locked"（外部只读连接/杀毒扫描持有共享锁，journal_mode=OFF
  的提交要独占锁），整片前功尽弃。改为文本追加后彻底没有锁竞争；
  收尾由**单进程**把 JSONL 灌成清单库。

用法:
    python scan_inorg.py [--workers 6] [--limit N] [--fresh] [--merge-only]
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sqlite3
import sys
import time
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("OMP_NUM_THREADS", "1")

import build_dbs as B  # noqa: E402

# 数据根目录：优先环境变量 XRD_DATA_ROOT，否则取本仓库根（tools/ 的上一级）。
# 本仓库（XinMatch-DataBase）不带 COD/PCOD/TCOD 原始归档，实际使用时一般指到
# 归档所在目录，例如：
#   XRD_DATA_ROOT=D:/Project/XRD/PolyXRD-DataBase python tools/scan_inorg.py --workers 6
# 它会读 <root>/COD|PCOD|TCOD，把无机副本写到 <root>/inorg_cif/{LIB}/，
# 清单写到 <root>/inorg_manifest.sqlite。
BASE = os.environ.get("XRD_DATA_ROOT") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
PARTS = os.path.join(BASE, "_scan_parts")
MANIFEST = os.path.join(BASE, "inorg_manifest.sqlite")

LIBS = [
    ("COD", os.path.join(BASE, "COD"), 3, 0),
    ("PCOD", os.path.join(BASE, "PCOD"), 2, 20_000_000),
    ("TCOD", os.path.join(BASE, "TCOD"), 3, 100_000_000),
]

_RE_SUM = re.compile(r"^\s*_chemical_formula_sum\s+(.*?)\s*$", re.M)
_RE_SEMI = re.compile(r"^\s*;\s*(.*?)\s*^\s*;", re.M | re.S)
_RE_DIGIT = re.compile(r"(\d+)")

MANIFEST_SCHEMA = """
CREATE TABLE IF NOT EXISTS files (
    cod_id      INTEGER PRIMARY KEY,
    lib         TEXT NOT NULL,
    rel         TEXT NOT NULL,
    formula     TEXT,
    formula_red TEXT,
    elements    TEXT,
    is_inorg    INTEGER NOT NULL,
    size        INTEGER,
    mtime       INTEGER
);
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
"""


def read_text(raw: bytes) -> str:
    try:
        return raw.decode("utf-8")
    except UnicodeDecodeError:
        return raw.decode("latin-1", "replace")


def extract_formula(text: str) -> str:
    """只取 _chemical_formula_sum，不做全 CIF 解析。"""
    m = _RE_SUM.search(text)
    if m:
        v = m.group(1).strip()
        if v[:1] in ("'", '"'):
            q = v[0]
            end = v.find(q, 1)
            return v[1:end].strip() if end > 0 else v.strip(q)
        return v
    m = _RE_SEMI.search(text)
    return m.group(1).strip() if m else ""


def unit_jobs():
    jobs = []
    for name, root, depth, off in LIBS:
        for rel, abs_dir in B.list_dirs(root, depth):
            jobs.append((name, rel, abs_dir, off))
    return jobs


def worker(args):
    wid, jobs, fresh = args
    os.makedirs(PARTS, exist_ok=True)
    frows = os.path.join(PARTS, "w%02d.jsonl" % wid)
    funits = os.path.join(PARTS, "w%02d.units" % wid)
    if fresh:
        for p in (frows, funits):
            if os.path.exists(p):
                os.remove(p)

    # 续跑：已完成的目录单元（units 文件每行一个）
    done = set()
    if os.path.exists(funits):
        with open(funits, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    done.add(line)

    # 追加模式打开；若上次被杀留下半行，截掉它再续写
    if os.path.exists(frows):
        with open(frows, "rb") as fh:
            fh.seek(max(0, os.path.getsize(frows) - 65536))
            tail = fh.read()
        if tail and not tail.endswith(b"\n"):
            cut = tail.rfind(b"\n")
            if cut >= 0:
                with open(frows, "r+b") as fh:
                    fh.truncate(os.path.getsize(frows) - (len(tail) - cut - 1))

    out = open(frows, "a", encoding="utf-8")
    uf = open(funits, "a", encoding="utf-8")

    n_file = n_inorg = n_copy = n_skip = n_nodigit = n_unit = 0
    buf = []
    t0 = time.time()
    for lib, rel, abs_dir, off in jobs:
        key = ("%s/%s" % (lib, rel)) if rel else lib
        if key in done:
            continue
        try:
            names = sorted(x for x in os.listdir(abs_dir) if x.lower().endswith(".cif"))
        except OSError:
            names = []
        for nm in names:
            src = os.path.join(abs_dir, nm)
            try:
                st = os.stat(src)
                if st.st_size == 0:
                    continue
                with open(src, "rb") as fh:
                    raw = fh.read()
            except OSError:
                continue
            n_file += 1
            formula = extract_formula(read_text(raw))
            elements, formula_red, counts = B.parse_elements(formula)
            inorg = 1 if B.is_inorganic(counts) else 0
            m = _RE_DIGIT.search(nm)
            if m:
                cod_id = int(m.group(1)) + off
            else:
                cod_id = 800_000_000 + wid * 1_000_000 + n_nodigit
                n_nodigit += 1
            if inorg:
                n_inorg += 1
                dst_dir = os.path.join(BASE, "inorg_cif", lib, rel) if rel else \
                    os.path.join(BASE, "inorg_cif", lib)
                dst = os.path.join(dst_dir, nm)
                try:
                    need = True
                    if os.path.exists(dst):
                        need = os.stat(dst).st_size != st.st_size
                    if need:
                        if not os.path.isdir(dst_dir):
                            os.makedirs(dst_dir, exist_ok=True)
                        with open(dst, "wb") as o:
                            o.write(raw)
                        n_copy += 1
                    else:
                        n_skip += 1
                except OSError:
                    pass
            buf.append(json.dumps([cod_id, lib,
                                   ("%s/%s" % (rel, nm)) if rel else nm,
                                   formula[:250], formula_red[:250], elements,
                                   inorg, st.st_size, int(st.st_mtime)],
                                  ensure_ascii=False))
            if len(buf) >= 2000:
                out.write("\n".join(buf) + "\n")
                out.flush()
                buf.clear()
        if buf:
            out.write("\n".join(buf) + "\n")
            out.flush()
            buf.clear()
        uf.write(key + "\n")
        uf.flush()
        n_unit += 1
        if n_unit % 200 == 0:
            el = time.time() - t0
            print("  w%02d %d/%d 单元 | 文件 %d | 无机 %d | 新复制 %d | %.1f 文件/s"
                  % (wid, n_unit, len(jobs), n_file, n_inorg, n_copy,
                     n_file / max(el, 1e-9)), flush=True)

    if buf:
        out.write("\n".join(buf) + "\n")
        out.flush()
    out.close()
    uf.close()
    return (wid, n_file, n_inorg, n_copy, n_skip, round(time.time() - t0, 1))


def merge(final=MANIFEST):
    """单进程把各 worker 的 JSONL 灌成清单库（无并发写，不会锁）。

    注意：JSONL 里存在"双跑"叠加的**重复行**（同一单元被两个进程池各扫一次）。
    因此先灌进无主键的临时表 raw 做行数/冲突统计，再按 cod_id 折叠进 files。
    """
    import glob
    if os.path.exists(final):
        os.remove(final)
    con = sqlite3.connect(final)
    con.executescript(MANIFEST_SCHEMA)
    con.execute("PRAGMA journal_mode=OFF")
    con.execute("PRAGMA synchronous=OFF")
    con.execute("CREATE TABLE raw (cod_id INTEGER, lib TEXT, rel TEXT, "
                "formula TEXT, formula_red TEXT, elements TEXT, "
                "is_inorg INTEGER, size INTEGER, mtime INTEGER)")
    ins_raw = "INSERT INTO raw VALUES (?,?,?,?,?,?,?,?,?)"
    tot = 0
    for f in sorted(glob.glob(os.path.join(PARTS, "w*.jsonl"))):
        buf = []
        with open(f, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                try:
                    buf.append(json.loads(line))
                except ValueError:      # 半行残留，跳过
                    continue
                if len(buf) >= 20000:
                    con.executemany(ins_raw, buf)
                    tot += len(buf)
                    buf.clear()
        if buf:
            con.executemany(ins_raw, buf)
            tot += len(buf)
        con.commit()

    # 真实 ID 冲突：同一 cod_id 来自不同的 (lib, rel) —— 必须 0
    clash = con.execute(
        "SELECT COUNT(*) FROM (SELECT cod_id FROM raw "
        "GROUP BY cod_id HAVING COUNT(DISTINCT lib || '/' || rel) > 1)"
    ).fetchone()[0]

    con.execute("INSERT OR REPLACE INTO files "
                "SELECT cod_id, lib, rel, formula, formula_red, elements, "
                "is_inorg, size, mtime FROM raw GROUP BY cod_id")
    con.execute("DROP TABLE raw")
    con.commit()

    seen_units = set()
    for f in sorted(glob.glob(os.path.join(PARTS, "w*.units"))):
        with open(f, "r", encoding="utf-8") as fh:
            for line in fh:
                line = line.strip()
                if line:
                    seen_units.add(line)
    con.execute("CREATE INDEX IF NOT EXISTS idx_manifest_lib ON files(lib)")
    con.execute("CREATE INDEX IF NOT EXISTS idx_manifest_inorg ON files(is_inorg)")
    for k, v in (("rows_read", str(tot)), ("units", str(len(seen_units))),
                 ("id_clash", str(clash)),
                 ("built_at", time.strftime("%Y-%m-%d %H:%M:%S"))):
        con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (k, v))
    con.commit()

    uniq = con.execute("SELECT COUNT(*) FROM files").fetchone()[0]
    ino = con.execute("SELECT COUNT(*) FROM files WHERE is_inorg=1").fetchone()[0]
    per = con.execute("SELECT lib, COUNT(*), SUM(is_inorg) FROM files "
                      "GROUP BY lib ORDER BY lib").fetchall()
    for lib, n, ni in per:
        con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)",
                    ("n_%s" % lib, str(n)))
        con.execute("INSERT OR REPLACE INTO meta VALUES (?,?)",
                    ("n_inorg_%s" % lib, str(ni or 0)))
    con.commit()
    con.close()
    print("  灌库: 读到 %d 行 -> 唯一条目 %d (折叠重复 %d)"
          % (tot, uniq, tot - uniq))
    print("  真实 ID 冲突(同一 cod_id 不同来源): %d  (必须为 0)" % clash)
    return uniq, ino, clash, len(seen_units)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--fresh", action="store_true")
    ap.add_argument("--merge-only", action="store_true")
    args = ap.parse_args()

    if not args.merge_only:
        jobs = unit_jobs()
        if args.limit:
            jobs = jobs[:args.limit]
        n = args.workers
        chunks = [jobs[i::n] for i in range(n)]
        print("目录单元 %d 个, %d 进程并行" % (len(jobs), n), flush=True)
        with Pool(n) as pool:
            for wid, nf, ni, nc, ns, sec in pool.imap_unordered(
                    worker, [(i, c, args.fresh) for i, c in enumerate(chunks)]):
                print("  w%02d 完成: 文件 %d, 无机 %d, 新复制 %d, 已存在 %d, %d s"
                      % (wid, nf, ni, nc, ns, sec), flush=True)

    tot, ino, dup, units = merge()
    print("=" * 62)
    print("清单条目总数 : %d" % tot)
    print("无机条目     : %d (%.1f%%)" % (ino, 100.0 * ino / max(tot, 1)))
    print("已处理单元   : %d" % units)
    print("cod_id 重复数: %d  (必须为 0)" % dup)
    print("清单: %s" % MANIFEST)


if __name__ == "__main__":
    main()
