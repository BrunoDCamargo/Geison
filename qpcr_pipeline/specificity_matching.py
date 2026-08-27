"""Deterministic exhaustive oligo matching for off-target specificity."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from qpcr_pipeline.config import SpecificityConfig, validate_specificity_config
from qpcr_pipeline.iupac import (
    IupacError,
    mismatch_positions,
    normalize_iupac,
    reverse_complement_iupac,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo


class SpecificityMatchingError(RuntimeError):
    """Raised when deterministic specificity matching cannot be completed."""


OligoRole = Literal["FORWARD", "PROBE", "REVERSE"]
Orientation = Literal["FORWARD", "REVERSE_COMPLEMENT"]
_ROLE_ORDER: tuple[OligoRole, ...] = ("FORWARD", "PROBE", "REVERSE")
_ORIENTATION_ORDER: tuple[Orientation, ...] = ("FORWARD", "REVERSE_COMPLEMENT")


@dataclass(frozen=True, slots=True)
class MatchHit:
    dataset_name: str
    assay_id: str
    sequence_id: str
    role: OligoRole
    orientation: Orientation
    source_start: int
    source_end: int
    oriented_start: int
    oriented_end: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class GeometryAmplicon:
    dataset_name: str
    assay_id: str
    sequence_id: str
    orientation: Orientation
    forward: MatchHit
    reverse: MatchHit
    probes: tuple[MatchHit, ...]
    source_start: int
    source_end: int
    amplicon_size: int
    primer_amplicon_plausible: bool
    detectable_off_target: bool


def _oligo_for_role(assay: AssayCandidate, role: OligoRole) -> DesignedOligo:
    try:
        return {
            "FORWARD": assay.forward_primer,
            "PROBE": assay.probe,
            "REVERSE": assay.reverse_primer,
        }[role]
    except (AttributeError, KeyError) as error:
        raise SpecificityMatchingError(
            f"Specificity role {role!r} is invalid for assay matching."
        ) from error


def _context(dataset_name: str, assay_id: str, sequence_id: str, role: OligoRole) -> str:
    return (
        f"dataset {dataset_name!r}, assay {assay_id!r}, "
        f"sequence {sequence_id!r}, role {role!r}"
    )


def _compatible_detail(
    oligo_sequence: str,
    target_sequence: str,
    role: OligoRole,
    config: SpecificityConfig,
) -> tuple[tuple[int, ...], bool, bool]:
    positions = mismatch_positions(oligo_sequence, target_sequence)
    three_prime = role != "PROBE" and any(
        position > len(oligo_sequence) - config.primer_3_prime_bases
        for position in positions
    )
    allowed = (
        config.max_probe_mismatches if role == "PROBE" else config.max_primer_mismatches
    )
    compatible = len(positions) <= allowed and not (
        role != "PROBE"
        and config.reject_primer_3_prime_mismatch
        and three_prime
    )
    return positions, three_prime, compatible


def enumerate_compatible_hits(
    dataset_name: str,
    record: LocalSequenceRecord,
    assay: AssayCandidate,
    role: OligoRole,
    config: SpecificityConfig,
) -> tuple[MatchHit, ...]:
    """Return every compatible site in both record orientations, without truncation."""
    validate_specificity_config(config)
    if not isinstance(dataset_name, str) or not dataset_name.strip():
        raise SpecificityMatchingError("Specificity dataset name must be non-blank.")
    if not isinstance(record, LocalSequenceRecord):
        raise SpecificityMatchingError(
            f"Specificity dataset {dataset_name!r} record must be a LocalSequenceRecord."
        )
    if not isinstance(assay, AssayCandidate):
        raise SpecificityMatchingError(
            f"Specificity dataset {dataset_name!r} assay must be an AssayCandidate."
        )
    if role not in _ROLE_ORDER:
        raise SpecificityMatchingError(f"Specificity oligo role {role!r} is invalid.")
    if not isinstance(record.sequence_id, str) or not record.sequence_id.strip():
        raise SpecificityMatchingError(
            f"Specificity dataset {dataset_name!r} has a record with a blank sequence ID."
        )

    context = _context(dataset_name, assay.assay_id, record.sequence_id, role)
    oligo = _oligo_for_role(assay, role)
    try:
        sequence = normalize_iupac(record.sequence, context=context)
        oligo_sequence = normalize_iupac(oligo.sequence, context=f"{context}, oligo")
    except (IupacError, TypeError) as error:
        raise SpecificityMatchingError(
            f"Invalid IUPAC input for {context}: {error}"
        ) from error

    oligo_length = len(oligo_sequence)
    record_length = len(sequence)
    if oligo_length == 0:
        raise SpecificityMatchingError(f"Specificity oligo is empty for {context}.")
    if record_length < oligo_length:
        return ()

    orientations: tuple[tuple[Orientation, str], ...] = (
        ("FORWARD", sequence),
        ("REVERSE_COMPLEMENT", reverse_complement_iupac(sequence)),
    )
    hits: list[MatchHit] = []
    for orientation, oriented_sequence in orientations:
        for oriented_start in range(1, record_length - oligo_length + 2):
            oriented_end = oriented_start + oligo_length - 1
            segment = oriented_sequence[oriented_start - 1 : oriented_end]
            target = (
                reverse_complement_iupac(segment)
                if role == "REVERSE"
                else segment
            )
            positions, three_prime, compatible = _compatible_detail(
                oligo_sequence, target, role, config
            )
            if not compatible:
                continue
            if orientation == "FORWARD":
                source_start, source_end = oriented_start, oriented_end
            else:
                source_start = record_length - oriented_end + 1
                source_end = record_length - oriented_start + 1
            hits.append(
                MatchHit(
                    dataset_name=dataset_name,
                    assay_id=assay.assay_id,
                    sequence_id=record.sequence_id,
                    role=role,
                    orientation=orientation,
                    source_start=source_start,
                    source_end=source_end,
                    oriented_start=oriented_start,
                    oriented_end=oriented_end,
                    mismatch_positions=positions,
                    mismatch_count=len(positions),
                    exact_match=not positions,
                    three_prime_mismatch=three_prime,
                    compatible=True,
                )
            )

    orientation_rank = {name: rank for rank, name in enumerate(_ORIENTATION_ORDER)}
    hits.sort(
        key=lambda hit: (
            orientation_rank[hit.orientation],
            hit.source_start,
            hit.source_end,
            hit.mismatch_count,
            hit.mismatch_positions,
        )
    )
    return tuple(hits)


def all_assay_hits(
    dataset_name: str,
    records: tuple[LocalSequenceRecord, ...],
    assays: tuple[AssayCandidate, ...],
    config: SpecificityConfig,
) -> tuple[MatchHit, ...]:
    """Return compatible hits in assay-major, record-major configured order."""
    validate_specificity_config(config)
    if not isinstance(records, tuple):
        raise SpecificityMatchingError("Specificity records must be a tuple.")
    if not isinstance(assays, tuple):
        raise SpecificityMatchingError("Specificity assays must be a tuple.")
    return tuple(
        hit
        for assay in assays
        for record in records
        for role in _ROLE_ORDER
        for hit in enumerate_compatible_hits(dataset_name, record, assay, role, config)
    )


def find_plausible_amplicons(
    hits: tuple[MatchHit, ...],
    config: SpecificityConfig,
) -> tuple[GeometryAmplicon, ...]:
    """Combine complete, untruncated hits into every plausible primer geometry."""
    validate_specificity_config(config)
    if not isinstance(hits, tuple):
        raise SpecificityMatchingError("Specificity hits must be a tuple.")

    groups: dict[tuple[str, str, str, Orientation], list[MatchHit]] = {}
    for hit in hits:
        if not isinstance(hit, MatchHit):
            raise SpecificityMatchingError("Specificity geometry requires MatchHit values.")
        key = (hit.dataset_name, hit.assay_id, hit.sequence_id, hit.orientation)
        groups.setdefault(key, []).append(hit)

    amplicons: list[GeometryAmplicon] = []
    for (dataset_name, assay_id, sequence_id, orientation), group_hits in groups.items():
        forward_hits = tuple(hit for hit in group_hits if hit.role == "FORWARD")
        reverse_hits = tuple(hit for hit in group_hits if hit.role == "REVERSE")
        probe_hits = tuple(hit for hit in group_hits if hit.role == "PROBE")
        for forward in forward_hits:
            for reverse in reverse_hits:
                if forward.oriented_end >= reverse.oriented_start:
                    continue
                amplicon_size = reverse.oriented_end - forward.oriented_start + 1
                if amplicon_size > config.max_amplicon_size:
                    continue
                probes = tuple(
                    probe
                    for probe in probe_hits
                    if forward.oriented_end
                    < probe.oriented_start
                    <= probe.oriented_end
                    < reverse.oriented_start
                )
                source_start = min(forward.source_start, reverse.source_start)
                source_end = max(forward.source_end, reverse.source_end)
                amplicons.append(
                    GeometryAmplicon(
                        dataset_name=dataset_name,
                        assay_id=assay_id,
                        sequence_id=sequence_id,
                        orientation=orientation,
                        forward=forward,
                        reverse=reverse,
                        probes=probes,
                        source_start=source_start,
                        source_end=source_end,
                        amplicon_size=amplicon_size,
                        primer_amplicon_plausible=True,
                        detectable_off_target=bool(probes),
                    )
                )
    return tuple(amplicons)
