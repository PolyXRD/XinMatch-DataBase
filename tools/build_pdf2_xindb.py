"""Build ``PDF2-<data-date>.xinDB`` — the ICDD PDF-2 2004 external library.

与 COD / PCOD / TCOD / MP 那几个库的根本差别
--------------------------------------------
那四个来源都是「CIF → 算峰」：库里的 d 与 I 是计算值，``meta.peak_convention``
记录的是算它们用的口径（λ、d 截止、是否 merge、强度归一化）。

PDF-2 不是。它本身就是一份**实测峰表**：没有 CIF、没有原子坐标，d 和 I 是几十
年文献里测出来的相对强度（I/I0 上限 100）。所以本工具做的是「读定宽流 → 按
``.xinDB`` v3 列存布局落盘」，中间**不做任何峰计算**，并且诚实地把这件事写进
``meta``：

===========================  ==================================================
``meta.peak_origin``         ``observed``（COD 家族是 ``computed``）
``meta.peak_convention``     ``observed-ICDD-PDF2-2004-norm1000``
``meta.dataset``             ``PDF2-2004``
``meta.name_policy``         ``mineral>compound``：``entries.mineral_name`` 优先放
                             ``X6`` 的矿物名（10.8%），缺失时退回主名称（化合物
                             名，约 100%）。列名是「名称位」，空着会让 GUI 只能
                             显示化学式。
===========================  ==================================================

强度口径仍是 ``norm1000``：ICDD 原表按 I/I0=100 归一，写入时
:class:`~xinmatch.core.db.writer_v3.XinDBWriter` 会按每相自身最大值乘到 1000，
与计算库共用同一把尺子，否则同一套打分权重在两个库上不可比。

不改写、不修补的数据
--------------------
* **晶胞 81.8%**：老 Hanawalt 卡片（1930~50 年代）本就没有晶胞，缺失是数据事实，
  不是解析缺陷。库内 ``a``~``gamma`` 与 ``volume`` 一起写 NULL，不填哨兵值。
* **空间群 72.8%**：同上，``space_group_num`` 为 0 的卡片只留 Pearson 符号。
* **密度**：不写。PDF-2 的 Z 与 Dx 在同一张卡片内可以互相矛盾（刚玉报 Z=2 但
  Dx 是按 Z=6 算的），算出来的密度会误导筛选，宁可 NULL。
* **峰数上限 199**：不用 ``peak_compute.DEFAULT_MAX_PEAKS``（150）截断 —— 那个
  上限是算峰引擎的安全阀，对实测峰表没有意义。14.9% 的卡片超过 150 条。

用法::

    python -m tools.build_pdf2_xindb                      # 全量，默认输出目录
    python -m tools.build_pdf2_xindb --limit 2000 --out tmp/PDF2-slice.xinDB
    python -m tools.build_pdf2_xindb -v                   # 打进度

产物默认落在 ``D:\\Project\\XRD\\PolyXRD\\cod_data\\``（与旧的
``PDF2_2004.sqlite`` 同处，该目录被 .gitignore 显式屏蔽）。
**ICDD 数据不得进入任何可发布目录**：``XinMatch-DataBase/data/`` 会被
``tools/package_release.py`` 默认 glob 进发布 zip，所以默认输出刻意不放在那里。
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from tools.build_cod_xindb import default_out_dir  # noqa: E402
from xinmatch.core.db.cif_meta import PhaseMeta, _crystal_system  # noqa: E402
from xinmatch.core.db.pdf2_scan import (  # noqa: E402
    Pdf2Card,
    elements_from_pdf2_formula,
    iter_cards,
)
from xinmatch.core.db.schema_v3 import XINDB_SUFFIX  # noqa: E402
from xinmatch.core.db.writer_v3 import XinDBWriter  # noqa: E402

logger = logging.getLogger("tools.build_pdf2_xindb")

#: 原始文件默认位置（ICDD PDF-2 2004 定宽流）。
DEFAULT_RAW = Path(r"C:\pdf2\pdf2.dat")

#: 默认输出目录：xinDB 统一落库点（与 COD/TCOD/PCOD/MP 同处），且被 /data/ 屏蔽。
DEFAULT_OUT_DIR = Path(r"D:\Project\XRD\XinMatch-DataBase\data")

#: 实测峰表的峰口径标签。刻意与 ``peak_compute.PEAK_CONVENTION`` 不同：
#: 那个常量描述「怎么从 CIF 算峰」，对实测峰表无意义；照抄就是让文件说谎。
OBSERVED_CONVENTION = "observed-ICDD-PDF2-2004-norm1000"

#: ``entries.tier``：0 = 实验测定，1 = 未标注，2 = 计算/理论（MP 用 2）。
TIER_EXPERIMENTAL = 0

#: Pearson 符号首字母 → 晶系。仅在没有空间群号时兜底使用。
_PEARSON_SYSTEM = {
    "a": "triclinic",
    "m": "monoclinic",
    "o": "orthorhombic",
    "t": "tetragonal",
    "h": "hexagonal",
    "c": "cubic",
}

#: 文献引用里的年份 token（``ANCHAM  10  457 1938 Hanawalt et al.``）。
_REF_YEAR_LOW = 1800
_REF_YEAR_HIGH = 2100


def _looks_like_coden(tok: str) -> bool:
    """期刊缩写（CODEN）形态：4~8 位大写字母/数字，且至少含一个字母。

    ``DANASG`` / ``ZEKRDZ`` / ``NBSMA6`` 是 CODEN；``990182`` 这种纯数字是 ICDD
    的原始文献编号，不是期刊。
    """
    if not 4 <= len(tok) <= 8:
        return False
    if not all(ch.isalnum() and (ch.isdigit() or ch.isupper()) for ch in tok):
        return False
    return any(ch.isalpha() for ch in tok)


def _split_reference(ref: str) -> tuple[str | None, int | None, str | None]:
    """把 ``X9`` 的文献串拆成 ``(journal, year, author)``.

    PDF-2 的 ``X9`` 是 CODEN 体例，实测形态只有三种::

        ZEKRDZ 198  177 1992                # CODEN 卷 页 年
        NBSMA6  13   35 1976                # CODEN 卷 页 年（无作者）
        DANASG                              # 只有 CODEN
        990182                              # 只有 ICDD 原始文献编号

    ``.xinDB`` 的 ``entries`` 只有 ``journal`` / ``year`` / ``author`` 三列，没有
    卷页字段，所以**卷与页不落库**（丢的是定位信息，不是内容判据）。拆不出年份
    时整串塞进 ``journal``，宁可字段语义粗糙，也不让溯源信息消失。

    Args:
        ref: ``X9`` 原始字符串（已去尾空格）。

    Returns:
        ``(journal, year, author)``；任一项拆不出时为 ``None``。
    """
    text = (ref or "").strip()
    if not text:
        return None, None, None
    tokens = text.split()
    year_idx = -1
    for i, tok in enumerate(tokens):
        digits = tok.strip(".,;")
        if digits.isdigit() and len(digits) == 4:
            value = int(digits)
            if _REF_YEAR_LOW <= value <= _REF_YEAR_HIGH:
                year_idx = i
                break

    if year_idx < 0:
        # 没有年份：整串（或那个孤零零的 CODEN）当来源文字保留。
        return text[:120], None, None

    head = tokens[:year_idx]
    journal = head[0] if head and _looks_like_coden(head[0]) else (" ".join(head) or None)
    author = " ".join(tokens[year_idx + 1 :]).strip() or None
    return journal, int(tokens[year_idx].strip(".,;")), author


def _crystal_system_of(card: Pdf2Card) -> str | None:
    """晶系：优先由空间群号推，其次用 Pearson 符号首字母。

    ``cif_meta._crystal_system`` 是七个晶系与 1~230 号段映射的**唯一**实现，
    这里直接复用而不是再抄一张表 —— 两张表迟早会分叉。
    """
    system = _crystal_system(card.space_group_num or None)
    if system is not None:
        return system
    pearson = (card.pearson or "").strip().lower()
    if pearson:
        return _PEARSON_SYSTEM.get(pearson[0])
    return None


def card_to_meta(card: Pdf2Card) -> PhaseMeta:
    """把一张卡片映射成 :class:`PhaseMeta`。

    ``entry_id`` 沿用 PolyXRD 的 ``PDF2-<6 位卡片号>`` ref_id 形态
    （``PDF2-010837``），这样同一张卡片在 PolyXRD 的旧索引与这里的新库里
    是同一个标识，迁移时不用另建映射表。
    """
    journal, year, author = _split_reference(card.reference)
    return PhaseMeta(
        entry_id=f"PDF2-{card.card_id}",
        formula=card.display_formula or None,
        space_group=card.space_group or None,
        crystal_system=_crystal_system_of(card),
        # 名称位：矿物名优先，缺失时退回化合物主名 —— 该列是 GUI 的显示名，
        # 只放 10.8% 有矿物名的卡片会让其余 89% 只能显示化学式。
        mineral_name=(card.mineral or card.name or None),
        elements=elements_from_pdf2_formula(card.display_formula),
        cell=card.cell,
        volume=card.cell_volume if card.cell is not None else None,
        density=None,  # Z 与 Dx 可互相矛盾，不臆造
        year=year,
        author=author,
        journal=journal,
        quality=(card.quality.strip() or None),
    )


def _sha256(path: Path) -> str:
    """算文件的 SHA-256（原始文件 563 MiB，约 1 秒）。"""
    digest = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def _default_data_date(raw: Path) -> str:
    """原始文件的日期（``YYYYMMDD``），作为 ``meta.data_date`` 默认值。

    ``pdf2.dat`` 的 mtime 就是 ICDD 2004 版的发布日期（2004-03-17），比「建库
    当天」更能回答「这份数据是哪一版」。两者语义不同，``build_time`` 另记。
    """
    return datetime.fromtimestamp(raw.stat().st_mtime).strftime("%Y%m%d")


@dataclass
class BuildStats:
    """一次建库的计数与哨兵。

    用 dataclass 而不是 ``dict[str, object]``：这些数字既要打印、又要写进
    ``meta``，用字典时每个取值点都得 ``int(...)`` 反推类型，mypy 也只能看到
    ``object``。
    """

    cards_read: int = 0
    cards_written: int = 0
    skipped_no_peaks: int = 0
    skipped_duplicate: int = 0
    no_elements: int = 0
    with_cell: int = 0
    with_space_group: int = 0
    vol_mismatch: int = 0
    max_peaks: int = 0
    total_peaks: int = 0
    #: 全部峰里的最小 / 最大 d（Å）；``d_min`` 初值取 ``inf`` 便于取最小值。
    d_min: float = 0.0
    d_max: float = 0.0
    elapsed_s: float = 0.0
    out: str = ""

    @property
    def reported_d_min(self) -> float:
        """写进 ``meta.peak_d_min`` 的值（一张卡都没有时给 0）。"""
        return 0.0 if self.d_min == float("inf") else round(self.d_min, 4)

    @property
    def reported_d_max(self) -> float:
        """写进 ``meta.peak_d_max`` 的值。"""
        return round(self.d_max, 4)

    def as_dict(self) -> dict[str, int | float | str]:
        """转成可直接打印 / 序列化的字典（数字保留 4 位的小数已取整）。"""
        return {
            "cards_read": self.cards_read,
            "cards_written": self.cards_written,
            "skipped_no_peaks": self.skipped_no_peaks,
            "skipped_duplicate": self.skipped_duplicate,
            "with_cell": self.with_cell,
            "with_space_group": self.with_space_group,
            "no_elements": self.no_elements,
            "vol_mismatch": self.vol_mismatch,
            "max_peaks": self.max_peaks,
            "total_peaks": self.total_peaks,
            "d_min": self.reported_d_min,
            "d_max": self.reported_d_max,
            "elapsed_s": self.elapsed_s,
            "out": self.out,
        }


def build(
    raw: Path,
    out: Path,
    *,
    data_date: str = "",
    limit: int | None = None,
    start: int = 0,
    log_every: int = 20000,
) -> dict[str, int | float | str]:
    """把 ``pdf2.dat`` 全量转成 ``.xinDB``，返回统计字典。

    Parameters
    ----------
    raw:
        原始定宽文件。
    out:
        目标 ``.xinDB`` 路径。
    data_date:
        快照日期 ``YYYYMMDD``；空串时取原始文件 mtime。
    limit, start:
        切片区间（验收用小库），``limit=None`` 表示全量。
    log_every:
        每处理多少张卡片打一行进度（``-v`` 时才可见）。

    Returns
    -------
    dict
        解析/落库统计，键名与 ``meta`` 里写回的哨兵一致。
    """
    raw = Path(raw)
    out = Path(out)
    if not raw.is_file():
        raise FileNotFoundError(f"PDF-2 原始文件不存在：{raw}")
    data_date = data_date or _default_data_date(raw)

    stats = BuildStats(d_min=float("inf"))
    t0 = time.time()

    writer = XinDBWriter(
        out,
        tier=TIER_EXPERIMENTAL,
        cif_source=str(raw),
        library_kind="PDF2",
        build_tool="tools.build_pdf2_xindb",
        data_date=data_date,
        extra_meta={
            # 这几条是本库与计算库的区别所在，必须写进文件本身。
            "peak_origin": "observed",
            "peak_convention": OBSERVED_CONVENTION,
            "dataset": "PDF2-2004",
            "pdf2_edition": "2004",
            "name_policy": "mineral>compound",
            # 卷/页在 v3 的 entries 里没有列位，明确写出这一层丢失，免得以后
            # 有人以为 journal 里那串数字是完整的文献。
            "reference_policy": "coden-year-author (卷/页无列位，未落库)",
            "raw_file": str(raw),
            "raw_bytes": raw.stat().st_size,
            "raw_sha256": _sha256(raw),
        },
    )

    try:
        for card in iter_cards(raw, limit=limit, start=start):
            stats.cards_read += 1
            if not card.peaks_d:
                # 无峰卡片对检索毫无意义，且 writer 不接受空 d 数组之外的
                # 语义 —— 直接跳过并计数（全库实测 1 张：531651）。
                stats.skipped_no_peaks += 1
                continue

            meta = card_to_meta(card)
            if not meta.elements:
                stats.no_elements += 1
            if card.cell is not None:
                stats.with_cell += 1
                reported = card.cell_vol_reported
                computed = card.cell_volume
                # 哨兵：报告体积与六参数自算体积差 10% 以上。PDF-2 里这类是
                # ICDD 的口径差异（菱面体原胞 / 超胞），列在这里是为了「数量
                # 对得上历史记录」，不是为了修数据。
                if reported and computed and abs(reported - computed) / reported > 0.10:
                    stats.vol_mismatch += 1
            if card.space_group_num:
                stats.with_space_group += 1

            writer.add(meta, np.asarray(card.peaks_d), np.asarray(card.peaks_i))

            n = card.n_peaks
            stats.total_peaks += n
            stats.max_peaks = max(stats.max_peaks, n)
            stats.d_min = min(stats.d_min, min(card.peaks_d))
            stats.d_max = max(stats.d_max, max(card.peaks_d))

            if log_every and stats.cards_read % log_every == 0:
                logger.info(
                    "解析 %s 张 / 已提交 %s 相 / %.0f s",
                    f"{stats.cards_read:,}",
                    f"{stats.cards_read - stats.skipped_no_peaks:,}",
                    time.time() - t0,
                )
    finally:
        # 这三条只有跑完才知道，但它们属于**文件自述**的一部分：验收器要靠
        # ``peak_max`` 判断「最大峰数是不是数据事实」（实测库不套算峰引擎的
        # DEFAULT_MAX_PEAKS），要 d 区间判断有没有越界。``extra_meta`` 在
        # close() 里落盘，所以在关闭前补进去。
        writer.extra_meta.update(
            {
                "peak_max": stats.max_peaks,
                "peak_d_min": stats.reported_d_min,
                "peak_d_max": stats.reported_d_max,
            }
        )
        closed = writer.close()

    stats.cards_written = int(closed.get("entries", 0))
    stats.skipped_duplicate = int(closed.get("skipped", 0))
    stats.elapsed_s = round(time.time() - t0, 1)
    stats.out = str(out)
    return stats.as_dict()


def build_arg_parser() -> argparse.ArgumentParser:
    """CLI 参数。"""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--raw", type=Path, default=DEFAULT_RAW, help=f"原始 pdf2.dat（默认 {DEFAULT_RAW}）")
    parser.add_argument("--out", type=Path, default=None, help="输出 .xinDB 路径")
    parser.add_argument("--data-date", default="", help="快照日期 YYYYMMDD（默认取原始文件 mtime）")
    parser.add_argument("--limit", type=int, default=None, help="只处理前 N 张卡片（切片验收）")
    parser.add_argument("--start", type=int, default=0, help="跳过前 N 张卡片")
    parser.add_argument("--log-every", type=int, default=20000, help="每 N 张打一行进度")
    parser.add_argument("-v", "--verbose", action="store_true", help="打开进度日志")
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI 入口。"""
    args = build_arg_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    raw = args.raw
    if not raw.is_file():
        print(f"[ERR] 找不到 PDF-2 原始文件：{raw}")
        return 2

    data_date = args.data_date or _default_data_date(raw)
    out = args.out
    if out is None:
        # 刻意不用 ``default_out_dir()``：那个函数会跟随 XINMATCH_XINDB_DIR
        # 落到发布用的 data/ 目录，而 ICDD 数据不能进发布包。
        out = DEFAULT_OUT_DIR / f"PDF2-{data_date}{XINDB_SUFFIX}"

    publishable = default_out_dir().resolve()
    if publishable in out.resolve().parents:
        print(
            f"[WARN] 输出目录 {out.parent} 是发布打包目录（tools/package_release.py "
            "默认 glob 全部 *.xinDB）。ICDD PDF-2 数据不得进入 git/Release，"
            "请改用本地目录。"
        )

    print(f"源文件 : {raw}  ({raw.stat().st_size / 1048576:.1f} MiB)")
    print(f"输出   : {out}")
    print(f"日期   : {data_date}")

    try:
        stats = build(
            raw,
            out,
            data_date=data_date,
            limit=args.limit,
            start=args.start,
            log_every=args.log_every,
        )
    except Exception as exc:  # pragma: no cover - CLI 兜底
        print(f"[ERR] 建库失败：{exc}")
        raise

    print("\n—— 建库统计 ——")
    for key in (
        "cards_read",
        "cards_written",
        "skipped_no_peaks",
        "skipped_duplicate",
        "with_cell",
        "with_space_group",
        "no_elements",
        "vol_mismatch",
        "max_peaks",
        "total_peaks",
        "d_min",
        "d_max",
        "elapsed_s",
    ):
        print(f"  {key:18} = {stats[key]}")
    size_mb = out.stat().st_size / 1048576 if out.is_file() else 0.0
    print(f"  {'file_size_mb':18} = {size_mb:.1f}")
    print("\n下一步：python -m tools.validate_xindb " + str(out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
