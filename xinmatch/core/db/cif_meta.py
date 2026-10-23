"""CIF metadata parsing and the canonical element bitmask (M1 / S1.3).

This module is the **single source of truth** for two things in the v3
``.xinDB`` format:

1. **Phase metadata** — :func:`parse_meta` turns a CIF file into a
   :class:`PhaseMeta` record (formula, space group, cell, density, elements,
   ...).  It is deliberately independent of the peak calculator so that the
   builder can stream metadata and reflections separately.

2. **The element bitmask** — :func:`element_mask` and
   :data:`element_bit_positions` define how a phase's element set is encoded
   as a 128-bit mask.  The writer stores that mask in the database and the
   reader uses the *same* table to translate a user's element restriction
   into a mask, so the two can never disagree.

Why the bitmask is 128 bit
--------------------------
The COD archive contains **97 distinct IUPAC elements** in the shipped
inorganic library — 118 exist in total.  A 64-bit mask can therefore not
represent every element that occurs in real data.  With a 64-bit table,
38.8 % of the shipped entries contain at least one element outside the
table (La, Pb, Bi, W, U, Ce, Nd, Ta, ... — all high-frequency heavy
elements), which would make the element filter systematically lose
recall.  The v3 format is not released yet, so the schema can simply carry
the mask in two 64-bit columns (``elem_mask_lo`` / ``elem_mask_hi``).

The bit positions are **frozen**
--------------------------------
The bit numbers are baked into the stored ``elem_mask_lo`` / ``elem_mask_hi``
columns.  If the table were ever recomputed from the data — even by one
position — every previously written database would be misread, with no
error raised: a phase would simply carry the bits of a different element.
The table below is therefore a hard-coded constant, it is versioned through
:data:`ELEMENT_MASK_SCHEME`, and the reader refuses a database whose
recorded scheme it does not recognise.

Element symbol whitelist
-------------------------
:func:`elements_from_formula` validates every parsed symbol against the
IUPAC whitelist.  This is not cosmetic: some COD entries carry a
*base64* blob in their ``_chemical_formula_sum`` field, and a naive
"uppercase letter + optional lowercase" regex turns such a blob into ~19
bogus "elements" (``A``, ``G``, ``Q``, ``Yw``, ...), which would then be
written into ``entries.elements`` and into the ``element_inv`` inverted
index.
"""

from __future__ import annotations

import hashlib
import math
import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path

import gemmi

__all__ = [
    "ELEMENT_MASK_SCHEME",
    "ELEMENT_BIT_TABLE_HASH",
    "IUPAC_ELEMENTS",
    "PhaseMeta",
    "element_bit_positions",
    "element_mask",
    "elements_from_formula",
    "elements_from_sites",
    "parse_meta",
]

#: Version tag of the frozen element bit table.  It is written to
#: ``meta.element_mask_scheme`` by the builder and checked by the reader,
#: which refuses a database built with a different table rather than
#: misreading its element masks.
ELEMENT_MASK_SCHEME = "iupac128-v1"

#: The 118 IUPAC element symbols, in atomic-number order.  This ordering is
#: what defines the bit numbers, so it must never change: see the module
#: docstring.
IUPAC_ELEMENTS: tuple[str, ...] = (
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)

#: ``element symbol -> bit number``.  Derived from :data:`IUPAC_ELEMENTS`
#: once at import time and then never touched — the frozen contract.
element_bit_positions: dict[str, int] = {
    symbol: index for index, symbol in enumerate(IUPAC_ELEMENTS)
}

#: Number of bits in the element mask.  128 bit covers all 118 IUPAC
#: elements with room to spare.
ELEMENT_MASK_BITS = 128

#: Digest of the bit table.  A reader can compare it with the digest
#: recorded in the database to detect a table that was edited in place
#: (which would silently invalidate every stored mask).
ELEMENT_BIT_TABLE_HASH: str = hashlib.sha256(
    "|".join(f"{symbol}:{index}" for symbol, index in element_bit_positions.items()).encode(
        "utf-8"
    )
).hexdigest()[:16]

#: Deuterium and tritium are hydrogen isotopes whose CIF formulas are written
#: as ``D`` / ``T``.  Neither is an IUPAC *element* symbol, so validating them
#: against the whitelist drops them -- and a pure deuteride (``D2 O``) then ends
#: up with an empty element list, which removes it from every element-filtered
#: search.  For X-ray scattering they are hydrogen, so they map to ``H``.
_HYDROGEN_ISOTOPES = {"D": "H", "T": "H"}

#: Fast membership test for the element whitelist.
_ELEMENT_SET = frozenset(IUPAC_ELEMENTS)

#: Regex for a chemical-formula token: an element symbol optionally
#: followed by a (possibly fractional) count.  The symbol is *validated*
#: against the whitelist afterwards, because a plain regex also matches the
#: capital letters of a base64 blob.
_FORMULA_TOKEN_RE = re.compile(r"([A-Z][a-z]?)(\d*\.?\d*)")

#: Formula strings that carry no chemistry: some COD entries store a
#: base64 payload in ``_chemical_formula_sum``.  They are recognised by
#: their character set and are rejected outright.
_BASE64_CHARS = set(
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/="
)


def _looks_like_base64(text: str) -> bool:
    """True when *text* is a base64 blob rather than a chemical formula.

    A chemical formula is made of element symbols, digits, spaces, dots and
    a few punctuation marks.  A base64 blob contains characters that never
    occur in a formula (``+``, ``/``, ``=``) or is a single long token with
    no separator, which is the giveaway for the COD entries that store an
    encoded attachment instead of a formula.

    Args:
        text: Candidate formula string.

    Returns:
        ``True`` when the string should be treated as garbage.
    """
    stripped = text.strip()
    if not stripped:
        return False
    if any(ch in stripped for ch in "+/="):
        return True
    # A single token longer than 12 characters with no space is not a
    # formula: even the longest real COD formula ("B6 La0.987 ...") is
    # shorter and contains spaces.
    if " " not in stripped and len(stripped) > 12:
        return True
    # More than half the characters being base64 letters that cannot start
    # a valid element symbol run is another giveaway.
    if len(stripped) >= 8 and not stripped[0].isalpha():
        return True
    return False


#: Leading element symbol of an ``_atom_site_type_symbol`` value.  CIF writes
#: oxidation states and disorder suffixes straight onto the symbol (``Fe3+``,
#: ``O2-``, ``K1+``), so the symbol is taken from the head of the token.
_SITE_SYMBOL_RE = re.compile(r"[A-Z][a-z]?")


def elements_from_sites(structure: object) -> list[str]:
    """Extract element symbols from the ``_atom_site`` list of a structure.

    This is the last-resort element source.  ``_chemical_formula_sum`` is the
    canonical field, but ~420 COD files omit every formula flavour while still
    carrying a complete atom-site list -- and a phase with no element list is
    invisible to every element-filtered search.  The site list is arguably more
    trustworthy than the formula anyway: it is what the peak calculation itself
    consumes.

    Args:
        structure: A ``gemmi.SmallStructure`` (or anything exposing ``sites``
            whose items carry a ``type_symbol``).

    Returns:
        Sorted IUPAC element symbols.  Empty when the structure carries no
        usable site.
    """
    found: set[str] = set()
    try:
        sites = structure.sites
    except Exception:
        return []
    for site in sites:
        match = _SITE_SYMBOL_RE.match(str(getattr(site, "type_symbol", "") or "").strip())
        if not match:
            continue
        symbol = _HYDROGEN_ISOTOPES.get(match.group(0), match.group(0))
        if symbol in _ELEMENT_SET:
            found.add(symbol)
    return sorted(found)


def elements_from_formula(formula: str | None) -> list[str]:
    """Extract the element symbols of a chemical formula.

    Every candidate symbol is validated against the IUPAC whitelist, so a
    formula that is really a base64 payload yields no elements instead of a
    dozen fictitious ones.

    Args:
        formula: The ``_chemical_formula_sum`` value, or ``None``.

    Returns:
        The sorted, de-duplicated element symbols (title case).  Empty when
        the formula is missing, unparseable, or not a formula at all.
    """
    if not formula:
        return []
    text = str(formula).strip()
    if not text or _looks_like_base64(text):
        return []
    found: set[str] = set()
    for symbol, _count in _FORMULA_TOKEN_RE.findall(text):
        symbol = _HYDROGEN_ISOTOPES.get(symbol, symbol)
        if symbol in _ELEMENT_SET:
            found.add(symbol)
    return sorted(found)


def element_mask(elements: Sequence[str]) -> tuple[int, int]:
    """Encode an element set as the 128-bit ``(low, high)`` mask.

    Args:
        elements: Element symbols of one phase.  Symbols are matched
            case-insensitively and must be IUPAC symbols; unknown ones
            contribute no bit.

    Returns:
        A ``(lo, hi)`` pair of unsigned 64-bit integers.  Bit ``n`` of the
        combined mask is set when the element whose bit number is ``n`` is
        present: ``n < 64`` lives in ``lo``, ``n >= 64`` in ``hi``.
    """
    lo = 0
    hi = 0
    for element in elements:
        bit = element_bit_positions.get(str(element).strip().title())
        if bit is None:
            continue
        if bit < 64:
            lo |= 1 << bit
        else:
            hi |= 1 << (bit - 64)
    return lo, hi


@dataclass
class PhaseMeta:
    """Scalar metadata of one phase, read from its CIF file.

    Attributes
    ----------
    entry_id:
        XinMatch entry id (``"96-XXX-XXXX"``), or ``""`` when unknown.
    formula:
        ``_chemical_formula_sum``, or ``None``.
    space_group:
        Hermann-Mauguin symbol of the space group, or ``None``.
    crystal_system:
        Crystal system derived from the space-group number, or ``None``.
    mineral_name:
        Mineral / common name, or ``None``.
    elements:
        Sorted IUPAC element symbols (see :func:`elements_from_formula`).
    cell:
        ``(a, b, c, alpha, beta, gamma)`` in A and degrees, or ``None``
        when the CIF carries no usable cell.
    volume:
        Unit-cell volume in A^3, or ``None``.
    density:
        Crystallographic density in g/cm^3, or ``None``.
    year, author, journal, quality:
        Bibliographic metadata, or ``None``.
    """

    entry_id: str = ""
    formula: str | None = None
    space_group: str | None = None
    crystal_system: str | None = None
    mineral_name: str | None = None
    elements: list[str] = field(default_factory=list)
    cell: tuple[float, float, float, float, float, float] | None = None
    volume: float | None = None
    density: float | None = None
    year: int | None = None
    author: str | None = None
    journal: str | None = None
    quality: str | None = None

    @property
    def element_mask(self) -> tuple[int, int]:
        """The 128-bit element mask of this phase."""
        return element_mask(self.elements)


#: ``_symmetry_cell_setting`` text to the crystal-system name the rest of the
#: library uses.  Only needed for phases whose space group cannot be resolved to
#: an international number -- magnetic structures and superspace groups record
#: their symmetry this way.  ``rhombohedral`` maps to ``trigonal`` because the
#: seven-system taxonomy folds the rhombohedral lattice into the trigonal one.
_CELL_SETTING_TO_SYSTEM: dict[str, str] = {
    "triclinic": "triclinic",
    "monoclinic": "monoclinic",
    "orthorhombic": "orthorhombic",
    "tetragonal": "tetragonal",
    "trigonal": "trigonal",
    "rhombohedral": "trigonal",
    "hexagonal": "hexagonal",
    "cubic": "cubic",
}


def _crystal_system(space_group_number: int | None) -> str | None:
    """Derive the crystal system from a space-group number.

    The seven crystal systems correspond to unambiguous ranges of the
    international space-group numbers 1-230.

    Args:
        space_group_number: The international number, or ``None``.

    Returns:
        The crystal-system name, or ``None`` when the number is unknown.
    """
    if space_group_number is None or not 1 <= space_group_number <= 230:
        return None
    if space_group_number <= 2:
        return "triclinic"
    if space_group_number <= 15:
        return "monoclinic"
    if space_group_number <= 74:
        return "orthorhombic"
    if space_group_number <= 142:
        return "tetragonal"
    if space_group_number <= 167:
        return "trigonal"
    if space_group_number <= 194:
        return "hexagonal"
    return "cubic"


def _first_block_value(block: gemmi.cif.Block, name: str) -> str | None:
    """Return the first value of a CIF data item, or ``None``.

    Args:
        block: The CIF block to read from.
        name: The data-item name, e.g. ``"_cell_length_a"``.

    Returns:
        The stripped value, or ``None`` when the item is absent or empty.
    """
    try:
        value = block.find_value(name)
    except Exception:
        return None
    if value is None:
        return None
    # gemmi returns the value as a plain string; the CIF placeholders "?"
    # (unknown) and "." (inapplicable) mean "no value".
    text = str(value).strip()
    if not text or text in {"?", "."}:
        return None
    # Keep the quotes out of the value.  gemmi preserves the CIF quoting, and
    # real COD files quote almost every text item (``_chemical_formula_sum``
    # and ``_symmetry_space_group_name_H-M`` are quoted in ~99 % of the files).
    # Leaving the quotes on poisons downstream parsing: the formula
    # "'C3 O6 Sr'" fails the element parser (and is even mistaken for a base64
    # blob), so every phase would silently lose its element set.
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        text = text[1:-1].strip()
        if not text:
            return None
    # CIF *text fields* are the third quoting form: a value that starts with
    # ``;`` on its own line and runs to the next line holding nothing but
    # ``;``.  gemmi hands the delimiters back verbatim, which broke ~1 % of the
    # COD phases -- ``";Al1.15 Ba0.04 ... Zn0.54\n;"`` starts with a semicolon,
    # so :func:`elements_from_formula` saw no element token and returned ``[]``.
    # Those phases then carried an empty element list, which silently removed
    # them from every element-filtered search.  A real formula never contains a
    # semicolon, so treating a leading/trailing one as a delimiter is safe.
    if text.startswith(";"):
        text = text[1:]
    if text.endswith(";"):
        text = text[:-1]
    # Collapse the internal line breaks a text field is allowed to contain;
    # without this the stored ``formula`` displays as garbage in the GUI.
    text = " ".join(text.split())
    if not text:
        return None
    return text


def _to_float(text: str | None) -> float | None:
    """Parse a CIF number, tolerating the standard uncertainty suffix.

    CIF stores a value as ``4.913(2)``; the parenthesised part is the
    standard uncertainty and must be stripped before parsing.

    Args:
        text: The raw CIF value, or ``None``.

    Returns:
        The parsed float, or ``None`` when the value is missing or not a
        number.
    """
    if not text:
        return None
    cleaned = text.split("(")[0].strip()
    try:
        value = float(cleaned)
    except ValueError:
        return None
    if not math.isfinite(value):
        return None
    return value


def _to_int(text: str | None) -> int | None:
    """Parse a CIF integer, tolerating the standard uncertainty suffix.

    Args:
        text: The raw CIF value, or ``None``.

    Returns:
        The parsed int, or ``None`` when the value is missing or not an
        integer.
    """
    if not text:
        return None
    cleaned = text.split("(")[0].strip()
    try:
        return int(float(cleaned))
    except ValueError:
        return None


def parse_meta(cif_path: str | Path, *, entry_id: str = "") -> PhaseMeta | None:
    """Read the scalar metadata of a CIF file.

    The file is parsed with :func:`gemmi.read_small_structure`, which is the
    same reader the peak calculator uses, so the two agree on the cell and
    the space group of a given file.

    Args:
        cif_path: Path of the CIF file.
        entry_id: XinMatch entry id to record on the result.

    Returns:
        The parsed metadata, or ``None`` when the file cannot be read or
        carries no cell (a phase without a cell has no reflections either).
    """
    path = Path(cif_path)
    try:
        ss = gemmi.read_small_structure(str(path))
    except Exception:
        # gemmi raises on malformed CIFs; a single bad file must never
        # abort a whole library build.
        return None
    if ss.cell is None:
        return None

    try:
        block = gemmi.cif.read(str(path)).sole_block()
    except Exception:
        block = None

    def item(name: str) -> str | None:
        if block is None:
            return None
        return _first_block_value(block, name)

    a = _to_float(item("_cell_length_a"))
    b = _to_float(item("_cell_length_b"))
    c = _to_float(item("_cell_length_c"))
    alpha = _to_float(item("_cell_angle_alpha"))
    beta = _to_float(item("_cell_angle_beta"))
    gamma = _to_float(item("_cell_angle_gamma"))
    cell: tuple[float, float, float, float, float, float] | None = None
    if None not in (a, b, c, alpha, beta, gamma):
        cell = (a, b, c, alpha, beta, gamma)  # type: ignore[arg-type]

    space_group = item("_symmetry_space_group_name_H-M") or item(
        "_space_group_name_H-M_alt"
    )
    if space_group:
        # The H-M symbol may carry a ":H" / ":R" setting suffix.
        space_group = space_group.split(":")[0].strip() or None

    sg_number: int | None = None
    try:
        sg_number = ss.spacegroup_number
    except Exception:
        sg_number = None
    if not sg_number:
        sg_number = _to_int(item("_space_group_it_number"))
    if not sg_number and space_group:
        # PCOD writes the Hermann-Mauguin symbol without the conventional
        # spaces (``"P1"``), which gemmi's small-structure reader does not
        # resolve -- it reports spacegroup_number 0, and the phase ends up with
        # no crystal system at all, i.e. invisible to the crystal-system
        # filter.  ``find_spacegroup_by_name`` does accept that spelling, so
        # fall back to it and adopt the canonical symbol it returns.
        try:
            resolved = gemmi.find_spacegroup_by_name(space_group)
        except Exception:  # pragma: no cover - gemmi API guard
            resolved = None
        if resolved is not None:
            sg_number = int(resolved.number)
            space_group = resolved.hm or space_group

    # Last resort for the crystal system: ``_symmetry_cell_setting`` states it
    # outright.  This is what saves the magnetic / superspace structures, whose
    # symmetry is written as ``_space_group_ssg_name  Imma(00\g)s00`` -- no IT
    # number, no resolvable H-M symbol, but the cell setting is spelled out.
    crystal_system = _crystal_system(sg_number)
    if crystal_system is None:
        crystal_system = _CELL_SETTING_TO_SYSTEM.get(
            (item("_symmetry_cell_setting") or "").strip().lower()
        )

    # ``_chemical_formula_sum`` is the canonical field, but ~460 COD files omit
    # it and only carry a structural or moiety formula.  Without a fallback those
    # phases end up with an empty element list and vanish from every
    # element-filtered search, so walk the standard CIF alternatives in order of
    # decreasing reliability.
    formula = item("_chemical_formula_sum")
    if elements_from_formula(formula) == []:
        for fallback in ("_chemical_formula_structural", "_chemical_formula_moiety"):
            candidate = item(fallback)
            if elements_from_formula(candidate):
                formula = candidate
                break
    mineral_name = item("_chemical_name_mineral") or item("_chemical_name_common")
    if mineral_name is not None and (not mineral_name or mineral_name.isdigit()):
        mineral_name = None

    year = _to_int(item("_journal_year"))
    density = _to_float(item("_cell_density"))

    volume: float | None = None
    try:
        if ss.cell.is_crystal():
            volume = float(ss.cell.volume)
    except Exception:
        volume = None
    if volume is None and cell is not None:
        # ``read_small_structure`` leaves ``ss.cell`` at its default (1, 1, 1,
        # 90, 90, 90) for a handful of COD files whose cell tags it cannot map,
        # so ``is_crystal()`` is False and the volume would be dropped even
        # though the six parameters were read fine from the text tags above.
        # Rebuild the cell from those parameters instead of losing the value.
        try:
            volume = float(gemmi.UnitCell(*cell).volume)
        except Exception:  # pragma: no cover - gemmi API guard
            volume = None

    return PhaseMeta(
        entry_id=entry_id,
        formula=formula,
        space_group=space_group,
        crystal_system=crystal_system,
        mineral_name=mineral_name,
        # Three-tier element resolution: the sum formula, then the structural /
        # moiety alternatives, and finally the atom-site list.  A phase without
        # an element set is unreachable through the element filter, so the
        # fallback matters more than its rarity suggests.
        elements=elements_from_formula(formula) or elements_from_sites(ss),
        cell=cell,
        volume=volume,
        density=density,
        year=year,
        author=item("_publ_author_name"),
        journal=item("_journal_name_full"),
        quality=item("_cod_database_code"),
    )
