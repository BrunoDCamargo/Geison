"""Bounded local oligo matching for qPCR assay inclusivity evaluation."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from qpcr_pipeline.config import InclusivityConfig
from qpcr_pipeline.iupac import (
    IupacError,
    mismatch_positions,
    normalize_iupac,
    reverse_complement_iupac,
)
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo


class InclusivityError(RuntimeError):
    """Raised when an inclusivity input cannot be evaluated safely."""


OligoRole = Literal["FORWARD", "PROBE", "REVERSE"]
Orientation = Literal["FORWARD", "REVERSE_COMPLEMENT"]


@dataclass(frozen=True, slots=True)
class OligoMatch:
    assay_id: str
    sequence_id: str
    role: OligoRole
    orientation: Orientation
    hit_rank: int
    source_start: int
    source_end: int
    expected_start: int
    expected_end: int
    displacement: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool
    selected: bool


@dataclass(frozen=True, slots=True)
class ProposedOligoCompatibility:
    role: OligoRole
    effective_sequence: str
    exact_match: bool
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    three_prime_mismatch: bool
    probe_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class AssayInclusivity:
    assay_id: str
    sequence_id: str
    orientation: Orientation | None
    geometry_found: bool
    source_amplicon_start: int | None
    source_amplicon_end: int | None
    amplicon_size: int | None
    forward_match: OligoMatch | None
    probe_match: OligoMatch | None
    reverse_match: OligoMatch | None
    original_compatible: bool
    proposed_forward: ProposedOligoCompatibility | None
    proposed_probe: ProposedOligoCompatibility | None
    proposed_reverse: ProposedOligoCompatibility | None
    proposed_compatible: bool


@dataclass(frozen=True, slots=True)
class OligoVariation:
    assay_id: str
    role: OligoRole
    oligo_position: int
    original_symbol: str
    original_support: tuple[str, ...]
    observed_symbol: str
    observed_support: tuple[str, ...]
    affected_sequence_ids: tuple[str, ...]
    affected_sequence_count: int
    affected_fraction: float | None
    primer_3_prime_position: bool


@dataclass(frozen=True, slots=True)
class DegeneracyProposal:
    assay_id: str
    role: OligoRole
    original_sequence: str
    proposed_sequence: str
    status: Literal["ACCEPTED", "UNCHANGED", "REJECTED"]
    reason: str
    original_degeneracy: int
    proposed_degeneracy: int
    changed_positions: tuple[int, ...]
    binding_site_count: int
    original_exact_count: int
    original_exact_fraction: float | None
    proposed_exact_count: int
    proposed_exact_fraction: float | None


@dataclass(frozen=True, slots=True)
class InclusivityResult:
    status: Literal["SKIPPED", "COMPLETE"]
    evaluation_sequence_ids: tuple[str, ...]
    oligo_matches: tuple[OligoMatch, ...]
    assay_results: tuple[AssayInclusivity, ...]
    variations: tuple[OligoVariation, ...]
    proposals: tuple[DegeneracyProposal, ...]
    oligo_matches_path: Path | None
    assay_inclusivity_path: Path | None
    oligo_variations_path: Path | None
    degeneracy_proposals_path: Path | None
    report_path: Path


@dataclass(frozen=True, slots=True)
class _Hit:
    public: OligoMatch
    oriented_start: int
    oriented_end: int
    target_in_synthesis_orientation: str


def _oligo_for_role(assay: AssayCandidate, role: OligoRole) -> DesignedOligo:
    return {
        "FORWARD": assay.forward_primer,
        "PROBE": assay.probe,
        "REVERSE": assay.reverse_primer,
    }[role]


def _enumerate_hits(
    *,
    assay_id: str,
    sequence_id: str,
    oriented_sequence: str,
    orientation: Orientation,
    role: OligoRole,
    oligo: DesignedOligo,
    config: InclusivityConfig,
) -> tuple[_Hit, ...]:
    """Return the bounded, deterministically ordered local sites for one oligo."""
    try:
        sequence = normalize_iupac(
            oriented_sequence,
            context=(
                f"record {sequence_id!r} for assay {assay_id!r} role {role!r}"
            ),
        )
        oligo_sequence = normalize_iupac(
            oligo.sequence,
            context=f"assay {assay_id!r} role {role!r} oligo",
        )
    except IupacError as error:
        raise InclusivityError(
            f"Invalid IUPAC input for assay {assay_id!r}, record {sequence_id!r}, "
            f"role {role!r}: {error}"
        ) from error

    oligo_length = len(oligo_sequence)
    record_length = len(sequence)
    if record_length < oligo_length:
        return ()

    window_start = max(1, oligo.reference_start - config.search_flank)
    window_end = min(record_length, oligo.reference_end + config.search_flank)
    first_start = window_start
    last_start = min(record_length - oligo_length + 1, window_end - oligo_length + 1)
    if first_start > last_start:
        return ()

    hits: list[tuple[_Hit, int, str]] = []
    for oriented_start in range(first_start, last_start + 1):
        oriented_end = oriented_start + oligo_length - 1
        target_segment = sequence[oriented_start - 1 : oriented_end]
        try:
            target_in_synthesis_orientation = (
                reverse_complement_iupac(target_segment)
                if role == "REVERSE"
                else target_segment
            )
            positions = mismatch_positions(oligo_sequence, target_in_synthesis_orientation)
        except IupacError as error:
            raise InclusivityError(
                f"Invalid IUPAC input for assay {assay_id!r}, record {sequence_id!r}, "
                f"role {role!r}: {error}"
            ) from error

        three_prime_mismatch = role != "PROBE" and any(
            position > oligo_length - config.primer_3_prime_bases
            for position in positions
        )
        compatible = (
            len(positions) <= (
                config.max_probe_mismatches
                if role == "PROBE"
                else config.max_primer_mismatches
            )
            and not (
                role != "PROBE"
                and config.reject_primer_3_prime_mismatch
                and three_prime_mismatch
            )
        )
        if orientation == "REVERSE_COMPLEMENT":
            source_start = record_length - oriented_end + 1
            source_end = record_length - oriented_start + 1
        else:
            source_start = oriented_start
            source_end = oriented_end

        public = OligoMatch(
            assay_id=assay_id,
            sequence_id=sequence_id,
            role=role,
            orientation=orientation,
            hit_rank=0,
            source_start=source_start,
            source_end=source_end,
            expected_start=oligo.reference_start,
            expected_end=oligo.reference_end,
            displacement=abs(oriented_start - oligo.reference_start),
            mismatch_positions=positions,
            mismatch_count=len(positions),
            exact_match=not positions,
            three_prime_mismatch=three_prime_mismatch,
            probe_mismatch=role == "PROBE" and bool(positions),
            compatible=compatible,
            selected=False,
        )
        hits.append(
            (
                _Hit(
                    public=public,
                    oriented_start=oriented_start,
                    oriented_end=oriented_end,
                    target_in_synthesis_orientation=target_in_synthesis_orientation,
                ),
                sum(
                    position > oligo_length - config.primer_3_prime_bases
                    for position in positions
                )
                if role != "PROBE"
                else 0,
                target_segment,
            )
        )

    hits.sort(
        key=lambda item: (
            item[0].public.mismatch_count,
            item[1],
            item[0].public.displacement,
            item[0].oriented_start,
            item[2],
        )
    )
    return tuple(
        replace(hit, public=replace(hit.public, hit_rank=rank))
        for rank, (hit, _, _) in enumerate(hits[: config.max_hits_per_oligo], 1)
    )
