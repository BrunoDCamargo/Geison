"""Core domain models for target, discovery, and evaluation sequence sets."""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TargetSequenceSet:
    """Sequences that belong to one homologous qPCR analysis target."""

    sequence_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class DiscoverySet:
    """Diversity-preserving subset used by computationally expensive discovery steps."""

    sequence_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class EvaluationSet:
    """All QC-approved target sequences used to evaluate candidate assays."""

    sequence_ids: tuple[str, ...]
