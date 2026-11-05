"""ICDD PDF-2 2004 定宽流解析（``C:\\pdf2\\pdf2.dat`` → :class:`Pdf2Card`）。

为什么需要这个模块
------------------
``.xinDB`` 的其它来源（COD / PCOD / TCOD / MP）都是「CIF → 算峰」，PDF-2 不是：
它本身就是一份**实测峰表**，没有 CIF、没有原子坐标，只有卡片号、化学式、晶胞
（81.8% 的卡片有）、空间群（72.8%）和 d-I 峰对。所以这里的活儿不是「算」，
而是「把定宽流读对」—— 一个字段读偏，整库的 d 值或强度就会系统性错位。

剖析口径的来源
--------------
下列列位不是猜出来的，是 PolyXRD 侧 ``services/pdf2_database.py`` 逐条核对
过的结论（配 ``tests/test_pdf2_space_group.py`` 的 85,013 条记录实测、以及
``docs/PDF2-2004重建验收报告-20261028.md`` 的 51 项验收）。本模块是那套规则
的移植，只把结果对象换成 dataclass、把「写 SQLite」的部分去掉 —— 建库目标从
``PDF2_2004.sqlite`` 换成了 ``.xinDB``，解析这一层不该有两份实现。

文件形态
--------
``pdf2.dat`` 是 **无换行的 80 字符定宽流**（587,999,120 B = 7,349,989 行）。
记录以 ``X1`` 子行开头，以 ``XK`` 结束；子记录标记固定落在第 71 列：

    P dddddd X <code> [nn]

第 1 位是数据集标识（P=在库 / D=已删除…），第 3 位随卡片递增（石英 470715
全是 X、方解石 010837 全是 R）—— **两个字母都会变**，硬编码 ``P...X...``
只会匹配到约 35% 的卡片，还会跨卡片串数据。

子记录类型（数据区见 :func:`_data_area`，默认截到第 67 列）：

===== ================================================================
code  内容
===== ================================================================
``1`` 主记录起始（卡片号 / 发布标记）
``3`` 空间群符号 + 空间群号 + Z/密度/**常规晶胞体积**（唯一截到 69 列的类型）
``4`` 分子量
``5`` 元素符号 / CAS 号
``6`` 名称，第 67~68 列是类别：``P``=主名称 ``M``=矿物名 ``C``=俗名 ``Z``=沸石名
``7`` 化学式（原始写法）  ``8`` 化学式（规范化）
``9`` 文献引用（CODEN 卷 页 年 作者）
``A`` Pearson 符号（``cF`` / ``hR`` / ``aP``…）
``D`` 约化晶胞 + 体积   ``E`` **常规晶胞**（对外暴露的是这个）
``F`` 辐射类型         ``G`` 质量标记
``I`` d-I 峰数据（每行 3 对，d=7 字符 + I=3 字符，步长 23）
``K`` 维护信息（记录结束）
===== ================================================================

用法::

    from xinmatch.core.db.pdf2_scan import iter_cards

    for card in iter_cards(r"C:\\pdf2\\pdf2.dat"):
        ...
"""

from __future__ import annotations

import math
import re
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path

from .cif_meta import IUPAC_ELEMENTS

#: 每行定宽（含行号与记录标记区）。
RECORD_LINE_SIZE = 80

#: ``XI`` 行里一对 d-I 的字段宽度与一行内的三对起始列。
_D_FIELD_SIZE = 7
_I_FIELD_SIZE = 3
_PAIR_OFFSETS = (0, 23, 46)

#: d 值的物理上限，超过它只能是读错列。
_D_MAX_SANE = 100.0

#: 原始文件里子记录标记的列位（0-based）。
_MARKER_COLUMN = 71

#: ``<数据集字母><6 位卡片号><序列字母><子记录类型><序号?>``。
_RE_RECORD_ID = re.compile(rf".{{{_MARKER_COLUMN}}}([A-Z])(\d{{6}})([A-Z])([A-Z0-9+*])(\d{{0,2}})")

#: 数据区右边界。第 67~70 列对多数子记录是标志区（孤立的 ``G``、
#: 续行标记 ``DB 1`` / ``SM 2``），混进化学式会把 ``G`` 当成元素污染元素过滤。
_DATA_END_DEFAULT = 67

#: 例外：类型 3 的晶胞体积恰好落在 63~68 列，按 67 截会把 ``365.97`` 砍成
#: ``365.``，体积全丢且建库哨兵恒为 0。全库 35,020 行逐行核对过。
_DATA_END_BY_CODE = {"3": 69}

#: 数值字段允许带尾随 ``E``（ICDD 用它标记「估计值」），如 ``2.680E``。
_RE_FLOAT_E = re.compile(r"^(\d*\.?\d+)E?$")

#: 类型 3 的列位。
_SG_SYM_SLICE = slice(0, 8)
_SG_FLAG_SLICE = slice(8, 11)
_SG_NUM_SLICE = slice(11, 14)
_SG_SETTING_INDEX = 14
_SG_TAIL_START = 15

#: 化学式里的候选元素符号（首字母大写 + 可选小写）。
_RE_ELEMENT_TOKEN = re.compile(r"[A-Z][a-z]?")

#: 氢的同位素写法归一（PDF-2 少数重水合物写 ``D2 O``）。
_ISOTOPE_TO_ELEMENT = {"D": "H", "T": "H"}

#: 元素白名单与 :mod:`cif_meta` 共用同一张 IUPAC 表，避免两套白名单分叉。
_ELEMENT_SET = frozenset(IUPAC_ELEMENTS)


@dataclass(slots=True)
class Pdf2Card:
    """一张 PDF-2 卡片的解析结果。

    Attributes
    ----------
    card_id:
        6 位卡片号（``"010837"`` 即 ICDD 的 ``01-0837``）。原文件里是整数，
        这里保留前导零，因为它同时被当作 ``.xinDB`` 的 ``entry_id`` 尾段。
    formula, formula_norm:
        ``X7`` 原始写法 / ``X8`` 规范化写法。前者为空时用后者兜底。
    name:
        ``X6`` 类别 ``P`` 的主名称（化合物名），覆盖率约 100%。
    mineral:
        ``X6`` 类别 ``M`` 的矿物名，覆盖率约 10.8% —— 老 Hanawalt 卡片本就没有。
    space_group, space_group_num, space_group_setting:
        ``X3`` 的空间群符号 / 国际编号 / 设置字母。编号为 0 表示未定，此时
        符号位放的是 Bravais 点阵字母（``P``/``F``/``R``…），本模块会清空它，
        否则 ``F`` 会被当成空间群 ``F``。
    cell:
        常规晶胞 ``(a, b, c, alpha, beta, gamma)``（``XE``），或 ``None``。
        注意**不是** ``XD`` 的约化晶胞：方解石的 E 胞是 4.983/4.983/17.019，
        D 胞是菱面体原胞（a=4.983, c=6.361, gamma=60）。
    cell_vol_reported:
        ``X3`` 末尾的报告体积，仅用于建库哨兵比对（986 条与 E 胞口径不同），
        **不入库**：库里的体积一律由六参数自算，保证与 a~gamma 自洽。
    pearson:
        ``XA`` 的 Pearson 符号（晶系 + 点阵），未定时可能是 ``mP?`` / ``t??``。
    radiation:
        ``XF`` 的辐射类型（``CuKa1`` / ``MoKa`` / ``FeKa1`` / ``DD``…）。
        **峰表是 d 值不是 2θ，所以与辐射无关**，这个字段只作溯源。
    quality:
        ``XG`` 的 ICDD 质量标记（``OO`` / ``V`` / ``SB`` / ``F`` / ``X``…）。
    reference:
        ``X9`` 的文献引用整串，形如 ``ANCHAM  10  457 1938 Hanawalt et al.``。
    peaks_d, peaks_i:
        ``XI`` 的 d（Å，降序）与 I/I0（0~100，ICDD 原尺度）。
    offset, size:
        该记录在原始文件中的字节偏移与长度，仅供溯源。
    """

    card_id: str = ""
    formula: str = ""
    formula_norm: str = ""
    name: str = ""
    mineral: str = ""
    space_group: str = ""
    space_group_num: int = 0
    space_group_setting: str = ""
    cell: tuple[float, float, float, float, float, float] | None = None
    cell_vol_reported: float | None = None
    pearson: str = ""
    radiation: str = ""
    quality: str = ""
    reference: str = ""
    peaks_d: list[float] = field(default_factory=list)
    peaks_i: list[float] = field(default_factory=list)
    offset: int = 0
    size: int = 0

    @property
    def n_peaks(self) -> int:
        """峰表行数。"""
        return len(self.peaks_d)

    @property
    def display_formula(self) -> str:
        """优先返回原始写法的化学式，缺失时用规范化写法。"""
        return self.formula or self.formula_norm

    @property
    def cell_volume(self) -> float | None:
        """由六参数自算的晶胞体积（Å³），无晶胞时为 ``None``。"""
        if self.cell is None:
            return None
        a, b, c, alpha, beta, gamma = self.cell
        return cell_volume(a, b, c, alpha, beta, gamma)


def cell_volume(
    a: float, b: float, c: float, alpha: float, beta: float, gamma: float
) -> float:
    """由六参数算晶胞体积（Å³）。

    用 ICDD 报告值反推是不可行的：986 条卡片的报告体积是菱面体原胞或超胞的
    体积，与六参数不同口径。自算才能保证 ``entries`` 里 a~gamma 与 volume
    始终自洽。
    """
    ca = math.cos(math.radians(alpha))
    cb = math.cos(math.radians(beta))
    cg = math.cos(math.radians(gamma))
    inside = 1.0 - ca * ca - cb * cb - cg * cg + 2.0 * ca * cb * cg
    if inside <= 0.0:
        return 0.0
    return a * b * c * math.sqrt(inside)


def _num(tok: str) -> float | None:
    """解析可能带 ``E`` 后缀的数值 token；非法返回 ``None``。"""
    m = _RE_FLOAT_E.match(tok)
    return float(m.group(1)) if m else None


def _data_area(line: str, code: str) -> str:
    """取一行子记录的数据区（右边界按类型选，见 :data:`_DATA_END_BY_CODE`）。"""
    return line[: _DATA_END_BY_CODE.get(code, _DATA_END_DEFAULT)].rstrip()


def _parse_peaks(line: str) -> list[tuple[float, float]]:
    """解析一行 ``XI``：3 对 (d, I)，跳过空位与越界值。"""
    out: list[tuple[float, float]] = []
    for offset in _PAIR_OFFSETS:
        d_txt = line[offset : offset + _D_FIELD_SIZE].strip()
        i_txt = line[offset + _D_FIELD_SIZE : offset + _D_FIELD_SIZE + _I_FIELD_SIZE].strip()
        if not d_txt:
            continue
        try:
            d_val = float(d_txt)
        except ValueError:
            continue
        if d_val <= 0.0 or d_val > _D_MAX_SANE:
            continue
        try:
            i_val = float(i_txt) if i_txt else 0.0
        except ValueError:
            i_val = 0.0
        # ICDD 的 I/I0 上限就是 100；越界只可能是列位串了。
        if i_val > 100.0:
            i_val = 100.0
        out.append((d_val, i_val))
    return out


def _parse_cell_fields(data: str, n_expected: int) -> tuple[float, ...] | None:
    """取一行里的前 *n_expected* 个数值作为 (a, b, c, α, β, γ)。

    类型 D / E 都是空格分隔的定宽数值，``split()`` 比硬编码列宽稳健。
    """
    toks = data.split()
    if len(toks) < n_expected:
        return None
    vals: list[float] = []
    for tok in toks[:n_expected]:
        v = _num(tok)
        if v is None:
            return None
        vals.append(v)
    return tuple(vals)  # type: ignore[return-value]


def elements_from_pdf2_formula(formula: str | None) -> list[str]:
    """从 PDF-2 的化学式写法里抽元素符号。

    为什么不直接用 :func:`xinmatch.core.db.cif_meta.elements_from_formula`
    ----------------------------------------------------------------
    那个函数带一层 base64 启发式（首字符不是字母、长度 ≥8 且无空格就判为
    编码块）—— 那条规则是为 COD 准备的，PDF-2 上有害：ICDD 的化学式惯例
    恰恰**会以数字或括号开头**

        ``3 Mg O !2 Si O2 !2 H2 O``（水合物的计量系数）
        ``( N H4 )2 W O4 !x H2 O``（括号分组 + 未知含水量）

    实测前 3000 张卡片里 15% 被这层启发式判成垃圾，元素列表为空的相在
    元素过滤检索里是**完全不可见**的。故这里只保留「符号必须落在 IUPAC
    白名单」这一层校验，跳过 base64 判定。

    Parameters
    ----------
    formula:
        ``X7`` / ``X8`` 的化学式写法。

    Returns
    -------
    list[str]
        去重升序的元素符号；写不出元素时返回空列表。
    """
    text = (formula or "").strip()
    if not text:
        return []
    found: set[str] = set()
    for symbol in _RE_ELEMENT_TOKEN.findall(text):
        if symbol in _ISOTOPE_TO_ELEMENT:
            symbol = _ISOTOPE_TO_ELEMENT[symbol]
        if symbol in _ELEMENT_SET:
            found.add(symbol)
    return sorted(found)


def _parse_symmetry(data: str, card: Pdf2Card) -> None:
    """解析类型 3：空间群 + 报告体积。

    ``cols 15~`` 里的 Z 与两个密度字段**有意不取**：约 20% 的卡片里 Z 来自
    原始文献、Dx 由 ICDD 另算，两者可以互相矛盾（刚玉报 Z=2 但 Dx=3.97 是按
    Z=6 的六方胞算的）。存一个自相矛盾的密度只会误导筛选。
    """
    symbol = data[_SG_SYM_SLICE].strip()
    num_txt = data[_SG_NUM_SLICE].strip()
    num = int(num_txt) if num_txt.isdigit() else 0

    if num > 0:
        card.space_group = symbol
        card.space_group_num = num
        card.space_group_setting = data[_SG_SETTING_INDEX].strip()
    else:
        # 编号 = 0 → 空间群未定，符号位此时是 Bravais 点阵字母，必须清空。
        card.space_group = ""
        card.space_group_num = 0
        card.space_group_setting = ""

    tail = data[_SG_TAIL_START:].split()
    if not tail:
        return
    # 报告体积是末字段且一律带小数点（94.67 / 2645.77）；Z 永远是纯整数，
    # 用小数点区分可避免把「只有 Z 没有体积」的行把 Z 当体积收下。
    if "." in tail[-1]:
        vol = _num(tail[-1])
        if vol is not None and vol > 0.0:
            card.cell_vol_reported = vol


def parse_record(lines: Sequence[str], offset: int) -> Pdf2Card | None:
    """解析一个完整记录（``X1`` 到 ``XK``）。

    Parameters
    ----------
    lines:
        该记录的全部 80 字符定宽行。
    offset:
        该记录在原始文件中的字节偏移（仅作溯源）。

    Returns
    -------
    Pdf2Card | None
        卡片号缺失时返回 ``None``。
    """
    card = Pdf2Card(offset=offset, size=len(lines) * RECORD_LINE_SIZE)

    for line in lines:
        m = _RE_RECORD_ID.match(line)
        if not m:
            continue
        card_id_digits = m.group(2)
        code = m.group(4)
        data = _data_area(line, code)

        if code == "3":
            _parse_symmetry(data, card)
        elif code == "6":
            name = data.strip()
            if name.startswith("$"):
                name = name[1:]
            kind = line[67:69].strip()
            if kind == "M":
                card.mineral = name
            elif kind in ("", "P") or not card.name:
                # ``P``（或空类别）是主名称，见到就覆盖；``C`` / ``Z`` 这类只在
                # 还没有主名称时兜底，避免俗名把主名称挤掉。
                card.name = name
        elif code == "7":
            card.formula = data.strip()
        elif code == "8":
            card.formula_norm = data.strip()
        elif code == "9":
            card.reference = data.strip()
        elif code == "A":
            toks = data.split()
            if toks:
                card.pearson = toks[0][:3]
        elif code == "D":
            # 约化晶胞：本模块**不保存**。PDF-2 的 D 胞与 E 胞在菱面体相上
            # 明显不同构（方解石 D=菱面体原胞、E=六方惯用胞），把两者并排放
            # 进同一个库只会让「a 到底指哪个胞」变成隐式假设。晶胞一律取 E。
            pass
        elif code == "E":
            card.cell = _parse_cell_fields(data, 6)  # type: ignore[assignment]
        elif code == "F":
            toks = data.strip().split()
            if toks:
                card.radiation = toks[0]
        elif code == "G":
            card.quality = data.strip()[:2]
        elif code == "I":
            for d_val, i_val in _parse_peaks(line):
                card.peaks_d.append(d_val)
                card.peaks_i.append(i_val)
        elif code == "K":
            break

        if code == "1" and card_id_digits:
            card.card_id = card_id_digits

    if not card.card_id:
        return None
    if not card.formula and card.formula_norm:
        card.formula = card.formula_norm
    return card


def iter_cards(
    raw_path: str | Path, *, limit: int | None = None, start: int = 0
) -> Iterator[Pdf2Card]:
    """逐记录迭代 ``pdf2.dat``。

    Parameters
    ----------
    raw_path:
        原始定宽文件（``C:\\pdf2\\pdf2.dat``，563 MiB）。
    limit:
        最多产出多少张卡片（切片验收用），``None`` 表示读到底。
    start:
        跳过前 *start* 张卡片（切片验收用）。

    Yields
    ------
    Pdf2Card
    """
    path = Path(raw_path)
    total_bytes = path.stat().st_size
    emitted = 0

    buf: list[str] = []
    rec_start = 0

    with path.open("rb") as fh:
        offset = 0
        while offset < total_bytes:
            chunk = fh.read(RECORD_LINE_SIZE * 8192)
            if not chunk:
                break
            n_lines = len(chunk) // RECORD_LINE_SIZE
            for i in range(n_lines):
                start_b = i * RECORD_LINE_SIZE
                line = chunk[start_b : start_b + RECORD_LINE_SIZE].decode(
                    "ascii", errors="replace"
                )
                m = _RE_RECORD_ID.match(line)
                if m and m.group(4) == "1" and buf:
                    card = parse_record(buf, rec_start)
                    if card is not None:
                        if emitted >= start:
                            yield card
                        emitted += 1
                        if limit is not None and emitted - start >= limit:
                            return
                    buf = []
                    rec_start = offset + start_b
                buf.append(line)
            offset += n_lines * RECORD_LINE_SIZE

    if buf:
        card = parse_record(buf, rec_start)
        if card is not None and emitted >= start:
            yield card


__all__ = [
    "Pdf2Card",
    "RECORD_LINE_SIZE",
    "cell_volume",
    "elements_from_pdf2_formula",
    "iter_cards",
    "parse_record",
]
