"""原地修正 PDF-2 库「X7 多行导致化学式只剩尾巴」的元数据缺陷。

背景
----
ICDD ``X7`` 长化学式会跨多行，旧 ``parse_record`` 用覆盖赋值，最终只留下
``)2`` / ``)0.43 ]n`` 这类尾巴；``display_formula`` 又优先返回它，于是这些卡
在 ``entries.elements`` / ``element_inv`` 里完全不可见（``validate_xindb`` 的
E2 WARN）。``pdf2_scan.py`` 已修：只在旧写法**抽不出元素**时退回 ``X8``
规范经验式或第一行，其余 163 K 张卡逐位不变。

本工具只改元数据，**不碰峰**：``formula`` / ``elements`` / ``elem_mask`` /
``element_inv``。源数据本来就没有化学式的 36 张卡保持 NULL，不编造。

用法
----
    python -m tools._patch_pdf2_formula --scan --db <xinDB> --diff-out <json>
    python -m tools._patch_pdf2_formula --apply --db <xinDB> --diff-out <json>
"""

from __future__ import annotations

import argparse
import json
import shutil
import sqlite3
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

from xinmatch.core.db.pdf2_scan import elements_from_pdf2_formula, iter_cards
from xinmatch.core.db.writer_v3 import pack_uint32

DEFAULT_SRC = Path(r"C:\pdf2\pdf2.dat")
_MASK63 = (1 << 63) - 1


def _ro(db: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{db}?immutable=1", uri=True)


def _element_bits(db: Path) -> dict[str, int]:
    con = _ro(db)
    text = con.execute(
        "SELECT value FROM meta WHERE key='element_bits'"
    ).fetchone()[0]
    con.close()
    return {
        part.split(":", 1)[0]: int(part.split(":", 1)[1])
        for part in str(text).split(";") if ":" in part
    }


def scan(db: Path, src: Path, diff_out: Path) -> int:
    """全量重解析，列出「旧式无元素、新式可解析」的可安全修补行。"""
    con = _ro(db)
    stored = {
        r[0].split("-", 1)[1]: (r[1], r[2])
        for r in con.execute("SELECT entry_id, formula, elements FROM entries")
    }
    total = len(stored)
    con.close()
    print(f"库内 entries: {total:,}　源: {src}")

    rows: list[dict[str, Any]] = []
    unresolved: list[str] = []
    t0 = time.perf_counter()
    for card in iter_cards(src):
        cid = card.card_id
        if cid not in stored:
            continue
        old_f, old_e = stored[cid]
        new_f = card.display_formula or None
        elems = elements_from_pdf2_formula(card.display_formula)
        new_e = ",".join(elems) or None
        if not elems:
            unresolved.append(cid)
        if (old_f or None) == new_f and (old_e or None) == new_e:
            continue
        # 护栏：只认「原本抽不出元素」→「现在能抽出」的卡。
        if old_e or not elems or (old_f or None) == new_f:
            print(f"[WARN] 跳过非预期变化 {cid}: {old_f!r}/{old_e!r} "
                  f"-> {new_f!r}/{new_e!r}")
            continue
        rows.append({
            "entry_id": f"PDF2-{cid}",
            "old_formula": old_f or None,
            "formula": new_f,
            "elements": new_e,
        })
    print(f"可安全修补 {len(rows):,} 行；源数据本就无组成 {len(unresolved):,} 行；"
          f"耗时 {time.perf_counter() - t0:.0f}s")
    diff_out.parent.mkdir(parents=True, exist_ok=True)
    diff_out.write_text(json.dumps({
        "db": str(db), "src": str(src), "total": total,
        "rows": rows, "unresolved": unresolved,
    }, ensure_ascii=False, indent=1), encoding="utf-8")
    print(f"diff -> {diff_out}")
    return 0


def apply(db: Path, diff_out: Path, backup: Path | None) -> int:
    payload = json.loads(diff_out.read_text(encoding="utf-8"))
    if Path(payload["db"]).resolve() != db.resolve():
        print("[ERR] diff 对应的库与 --db 不一致，拒绝写入")
        return 2
    rows: list[dict[str, Any]] = payload["rows"]
    if not rows:
        print("无可修补行。")
        return 0
    bits = _element_bits(db)
    if backup is not None:
        print(f"备份 {db} -> {backup} …")
        shutil.copy2(db, backup)
        if backup.stat().st_size != db.stat().st_size:
            print("[ERR] 备份大小不一致，拒绝写入")
            return 2
        print("备份完成。")

    con = sqlite3.connect(db)
    con.execute("PRAGMA journal_mode=WAL")
    try:
        con.execute("BEGIN IMMEDIATE")
        for item in rows:
            elems = str(item["elements"]).split(",")
            mask = 0
            for sym in elems:
                bit = bits.get(sym.strip())
                if bit is not None and bit < 63:
                    mask |= 1 << bit
            mask &= _MASK63
            cur = con.execute(
                "UPDATE entries SET formula=?, elements=?, elem_mask=? "
                "WHERE entry_id=? AND (elements IS NULL OR trim(elements)='') "
                "AND formula IS ?",
                (item["formula"], item["elements"], mask,
                 item["entry_id"], item["old_formula"]),
            )
            if cur.rowcount != 1:
                raise RuntimeError(
                    f"护栏触发：{item['entry_id']} rowcount={cur.rowcount}")
        # 从 entries 全表重建倒排（与 writer_v3._write_element_inv 同口径）。
        buckets: dict[str, list[int]] = defaultdict(list)
        for rowid, text in con.execute(
                "SELECT rowid, elements FROM entries ORDER BY rowid"):
            if not text:
                continue
            for sym in str(text).split(","):
                sym = sym.strip()
                if sym:
                    buckets[sym].append(int(rowid))
        con.execute("DELETE FROM element_inv")
        con.executemany(
            "INSERT INTO element_inv (element, eids) VALUES (?,?)",
            [(sym, pack_uint32(sorted(ids)))
             for sym, ids in sorted(buckets.items())],
        )
        count = con.execute("SELECT COUNT(*) FROM entries").fetchone()[0]
        if count != payload["total"]:
            raise RuntimeError(f"entries 行数变化 {payload['total']} -> {count}")
        con.commit()
        print(f"写入 {len(rows):,} 行 formula/elements/elem_mask；"
              f"element_inv 重建为 {len(buckets)} 个元素。")
    except Exception:
        con.rollback()
        raise
    finally:
        con.close()
    con = _ro(db)
    print("剩余无元素行:",
          con.execute("SELECT COUNT(*) FROM entries WHERE elements IS NULL "
                      "OR trim(elements)=''").fetchone()[0])
    print("剩余空式行:",
          con.execute("SELECT COUNT(*) FROM entries WHERE formula IS NULL "
                      "OR trim(formula)=''").fetchone()[0])
    con.close()
    return 0


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--src", type=Path, default=DEFAULT_SRC)
    ap.add_argument("--diff-out", type=Path,
                    default=Path(".scratch/pdf2_formula_diff.json"))
    ap.add_argument("--scan", action="store_true")
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--backup", type=Path, default=None)
    args = ap.parse_args(argv)
    if args.scan == args.apply:
        print("[ERR] 必须且只能选 --scan / --apply 之一")
        return 2
    if args.scan:
        return scan(args.db, args.src, args.diff_out)
    backup = args.backup or args.db.with_name(
        args.db.name + ".bak-" + time.strftime("%Y%m%d-%H%M%S"))
    return apply(args.db, args.diff_out, backup)


if __name__ == "__main__":
    raise SystemExit(main())
