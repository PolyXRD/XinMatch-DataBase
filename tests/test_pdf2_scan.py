"""``pdf2_scan`` 定宽流解析的回归测试。

分两层：

* **合成行**（无条件跑）—— 把列位知识钉死。定宽解析最贵的错误是「字段读偏
  一格」，而这类错误在没有真值对照时完全看不出来，所以每个易错列位都配一条
  用例：类型 3 的体积落在 63~68 列、空间群号为 0 时符号位放的是 Bravais 点阵
  字母、``XI`` 每行三对 d-I、名称按 67~68 列分类。
* **真实文件**（``C:\\pdf2\\pdf2.dat`` 存在时才跑）—— 用 ICDD 已知卡片自证：
  方解石 / 石英 / 硅 / 刚玉 / 石墨的卡片号、空间群、晶胞、最强线。这几个值
  与 PolyXRD 侧的验收报告逐字一致，两边跑的是同一份数据。

没有原始文件时整组跳过，而不是假装通过。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from xinmatch.core.db.pdf2_scan import (
    Pdf2Card,
    cell_volume,
    elements_from_pdf2_formula,
    iter_cards,
    parse_record,
)

RAW = Path(r"C:\pdf2\pdf2.dat")

requires_raw = pytest.mark.skipif(
    not RAW.is_file(), reason=f"ICDD PDF-2 原始文件不在本机：{RAW}"
)


# ---------------------------------------------------------------------------
# 合成行工具
# ---------------------------------------------------------------------------
def _line(data: str, code: str, card: str = "010837") -> str:
    """拼一条 80 字符定宽行：数据区 + 第 71 列起的子记录标记。"""
    return data.ljust(71)[:71] + f"P{card}R{code}"


def _parse(*lines: str) -> Pdf2Card:
    card = parse_record(list(lines), 0)
    assert card is not None
    return card


def _sym_line(head: str, tail: str) -> str:
    """类型 3 行：``head`` 占 0~62 列，末字段右对齐到第 68 列。"""
    return head.ljust(63) + tail.rjust(6)


# ---------------------------------------------------------------------------
# 定宽列位
# ---------------------------------------------------------------------------
def test_peaks_three_pairs_per_line() -> None:
    """``XI`` 一行三对 (d, I)，列位 0 / 23 / 46，d 宽 7、I 宽 3。"""
    data = (
        "  4.255 20".ljust(23)
        + "  3.343100".ljust(23)
        + "  2.457 12".ljust(23)
    )
    card = _parse(_line("", "1"), _line(data, "I"))
    assert card.peaks_d == [4.255, 3.343, 2.457]
    assert card.peaks_i == [20.0, 100.0, 12.0]


def test_peak_line_skips_empty_slots_and_absurd_d() -> None:
    """空位与 d>100 的读数跳过，不产生假峰。"""
    data = "  0.000  0".ljust(23) + "         ".ljust(23) + "150.000 50".ljust(23)
    card = _parse(_line("", "1"), _line(data, "I"))
    assert card.peaks_d == []


def test_symmetry_keeps_symbol_only_when_number_is_set() -> None:
    """空间群号 >0 时收符号；=0 时 0~7 列是 Bravais 点阵字母，必须清空。"""
    ok = _parse(_line("", "1"), _line(_sym_line("R-3c    E  167A", "365.97"), "3"))
    assert (ok.space_group, ok.space_group_num) == ("R-3c", 167)

    unknown = _parse(_line("", "1"), _line(_sym_line("F            0", ""), "3"))
    assert unknown.space_group == ""
    assert unknown.space_group_num == 0


def test_type3_volume_survives_the_69_column_cut() -> None:
    """体积字段落在 63~68 列：按默认 67 列截断会把它砍成 ``365.`` 而丢掉。

    这是历史上真实踩过的坑（PolyXRD 的 ``tests/test_pdf2_space_group.py``
    抓出），单独立一条防止有人把 ``_DATA_END_BY_CODE`` 删掉。
    """
    data = _sym_line("R-3c    E  167A", "365.97")
    assert data.index("365.97") == 63
    card = _parse(_line("", "1"), _line(data, "3"))
    assert card.cell_vol_reported == pytest.approx(365.97)


def test_type3_zero_volume_without_decimal_is_ignored() -> None:
    """只有 Z 没有体积时，末字段是纯整数，不能被当成体积。"""
    card = _parse(_line("", "1"), _line(_sym_line("Fd-3m   E  227A", "8"), "3"))
    assert card.cell_vol_reported is None


def test_cell_comes_from_type_e_not_type_d() -> None:
    """晶胞取 ``XE``（常规胞），不取 ``XD``（约化胞）。"""
    # 方解石：E 是六方惯用胞，D 是菱面体原胞 —— 两者必须能区分开。
    regular = "4.983  4.983  17.019  90.0  90.0  120.0"
    reduced = "4.983  4.983  6.361   90.0  90.0  120.0"
    card = _parse(_line("", "1"), _line(reduced, "D"), _line(regular, "E"))
    assert card.cell == pytest.approx((4.983, 4.983, 17.019, 90.0, 90.0, 120.0))


def test_name_kind_column_separates_mineral_from_compound() -> None:
    """第 67~68 列是名称类别：``M`` = 矿物名，``P`` = 主名称。"""
    main = "Calcium Carbonate Oxide".ljust(67) + " P"
    mineral = "Calcite".ljust(67) + " M"
    card = _parse(
        _line("", "1"), _line(main, "6"), _line(mineral, "6")
    )
    assert card.name == "Calcium Carbonate Oxide"
    assert card.mineral == "Calcite"


def test_dollar_prefix_is_stripped_from_name() -> None:
    """``$`` 是 ICDD 的斜体/希腊字母标记，不属于名称正文。"""
    card = _parse(_line("", "1"), _line("$Quartz".ljust(67) + " P", "6"))
    assert card.name == "Quartz"


def test_formula_falls_back_to_normalised_form() -> None:
    card = _parse(_line("", "1"), _line("Ca C O3".ljust(67) + " ", "8"))
    assert card.formula == "Ca C O3"


def test_missing_card_number_returns_none() -> None:
    """第 71 列起没有合法标记（缺 6 位卡片号）时不计为一张卡片。"""
    assert parse_record([_line("", "1", card="")], 0) is None


# ---------------------------------------------------------------------------
# 化学式 → 元素
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    ("formula", "expected"),
    [
        ("Ca C O3", ["C", "Ca", "O"]),
        ("Al2 O3", ["Al", "O"]),
        # 水合物计量系数以数字开头 —— 走 COD 的 base64 启发式会被判成垃圾
        ("3 Mg O !2 Si O2 !2 H2 O", ["H", "Mg", "O", "Si"]),
        # 括号分组 + 未知含水量
        ("( N H4 )2 W O4 !x H2 O", ["H", "N", "O", "W"]),
        ("", []),
        (")x ]n", []),
    ],
)
def test_elements_from_pdf2_formula(formula: str, expected: list[str]) -> None:
    assert elements_from_pdf2_formula(formula) == expected


def test_cell_volume_matches_reported_for_calcite() -> None:
    """方解石六方胞体积对得上 ICDD 报告值（365.97 Å³）。"""
    vol = cell_volume(4.983, 4.983, 17.019, 90.0, 90.0, 120.0)
    assert vol == pytest.approx(365.97, abs=0.05)


# ---------------------------------------------------------------------------
# 真实文件：ICDD 已知卡片
# ---------------------------------------------------------------------------
@requires_raw
def test_known_icdd_cards() -> None:
    """五张已知卡片的卡片号 / 空间群 / 晶胞 / 最强线。

    这些值来自 PolyXRD 侧 ``docs/PDF2-2004重建验收报告-20261028.md`` §6.1
    的点检表，两侧跑同一份原始文件，所以结果必须逐字相同。
    """
    want = {"010837", "461045", "271402", "461212", "010646"}
    found: dict[str, Pdf2Card] = {}
    for card in iter_cards(RAW):
        if card.card_id in want:
            found[card.card_id] = card
            if len(found) == len(want):
                break
    assert set(found) == want

    calcite = found["010837"]
    assert (calcite.space_group, calcite.space_group_num) == ("R-3c", 167)
    assert calcite.cell == pytest.approx((4.983, 4.983, 17.019, 90.0, 90.0, 120.0))
    assert calcite.mineral == "Calcite"
    assert calcite.name == "Calcium Carbonate Oxide"

    quartz = found["461045"]
    assert (quartz.space_group, quartz.space_group_num) == ("P3221", 154)
    assert quartz.n_peaks == 58
    strongest = max(zip(quartz.peaks_i, quartz.peaks_d, strict=True))
    # 原表写的是 3.3435（7 字符字段），报告里四舍五入到 3.343。
    assert strongest[0] == pytest.approx(100.0)
    assert strongest[1] == pytest.approx(3.3435, abs=1e-4)

    silicon = found["271402"]
    assert (silicon.formula, silicon.space_group_num) == ("Si", 227)

    corundum = found["461212"]
    assert (corundum.formula, corundum.space_group_num) == ("Al2 O3", 167)

    graphite = found["010646"]
    assert (graphite.formula, graphite.space_group_num) == ("C", 194)


@requires_raw
def test_full_scan_counts() -> None:
    """全量扫描的计数与 ICDD 2004 的公认总数吻合。"""
    n = 0
    no_peaks = 0
    with_cell = 0
    for card in iter_cards(RAW):
        n += 1
        if not card.peaks_d:
            no_peaks += 1
        if card.cell is not None:
            with_cell += 1
    assert n == 163_835
    assert no_peaks == 1
    # 比建库统计（133,935）多 1：唯一那张「无峰」卡片其实带晶胞，建库时被跳过，
    # 所以写进库的是 133,935。两个数都对，差的就是这一张。
    assert with_cell == 133_936
