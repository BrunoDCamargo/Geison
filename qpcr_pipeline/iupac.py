"""Dependency-free helpers for IUPAC nucleotide symbols."""

from collections.abc import Iterable


class IupacError(ValueError):
    """Raised when an invalid IUPAC nucleotide operation is requested."""


IUPAC_SUPPORT: dict[str, frozenset[str]] = {
    "A": frozenset(("A",)),
    "C": frozenset(("C",)),
    "G": frozenset(("G",)),
    "T": frozenset(("T",)),
    "R": frozenset(("A", "G")),
    "Y": frozenset(("C", "T")),
    "S": frozenset(("G", "C")),
    "W": frozenset(("A", "T")),
    "K": frozenset(("G", "T")),
    "M": frozenset(("A", "C")),
    "B": frozenset(("C", "G", "T")),
    "D": frozenset(("A", "G", "T")),
    "H": frozenset(("A", "C", "T")),
    "V": frozenset(("A", "C", "G")),
    "N": frozenset(("A", "C", "G", "T")),
}

SUPPORT_TO_IUPAC = {support: symbol for symbol, support in IUPAC_SUPPORT.items()}
IUPAC_COMPLEMENT: dict[str, str] = {
    "A": "T", "C": "G", "G": "C", "T": "A",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "V": "B",
    "D": "H", "H": "D", "N": "N",
}


def normalize_iupac(sequence: str, *, context: str) -> str:
    normalized = sequence.upper()
    for position, symbol in enumerate(normalized, 1):
        if symbol not in IUPAC_SUPPORT:
            raise IupacError(
                f"Invalid IUPAC symbol in {context} at position {position}: {symbol!r}."
            )
    return normalized


def iupac_support(symbol: str) -> frozenset[str]:
    normalized = symbol.upper()
    if len(normalized) != 1 or normalized not in IUPAC_SUPPORT:
        raise IupacError(f"Invalid IUPAC symbol {symbol!r}.")
    return IUPAC_SUPPORT[normalized]


def minimal_iupac_symbol(bases: Iterable[str]) -> str:
    support = frozenset(base.upper() for base in bases)
    try:
        return SUPPORT_TO_IUPAC[support]
    except KeyError as error:
        raise IupacError(
            "IUPAC support must be a non-empty subset of A, C, G, T."
        ) from error


def reverse_complement_iupac(sequence: str) -> str:
    normalized = normalize_iupac(sequence, context="sequence")
    return "".join(IUPAC_COMPLEMENT[symbol] for symbol in reversed(normalized))


def mismatch_positions(oligo: str, target: str) -> tuple[int, ...]:
    if len(oligo) != len(target):
        raise IupacError("Oligo and target lengths must match.")
    normalized_oligo = normalize_iupac(oligo, context="oligo")
    normalized_target = normalize_iupac(target, context="target")
    return tuple(
        index
        for index, (oligo_symbol, target_symbol) in enumerate(
            zip(normalized_oligo, normalized_target), 1
        )
        if not iupac_support(target_symbol) <= iupac_support(oligo_symbol)
    )


def sequence_degeneracy(sequence: str) -> int:
    result = 1
    for symbol in normalize_iupac(sequence, context="sequence"):
        result *= len(iupac_support(symbol))
    return result
