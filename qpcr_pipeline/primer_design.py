"""Pure candidate-region selection for Primer3 assay design."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Literal

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)


class PrimerDesignError(RuntimeError):
    """Raised when primer-design inputs cannot produce traceable results."""


@dataclass(frozen=True, slots=True)
class CandidateRegion:
    region_id: str
    rank: int
    reference_start: int
    reference_end: int
    peak_start: int
    peak_end: int
    position_count: int
    usable_length: int
    usable_fraction: float
    mean_conservation: float
    minimum_conservation: float
    mean_coverage: float
    mean_gap_frequency: float
    mean_entropy_bits: float


@dataclass(frozen=True, slots=True)
class DesignedOligo:
    sequence: str
    reference_start: int
    reference_end: int
    length: int
    tm: float
    gc_percent: float
    penalty: float | None
    metrics: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class AssayCandidate:
    assay_id: str
    region_id: str
    primer3_index: int
    forward_primer: DesignedOligo
    probe: DesignedOligo
    reverse_primer: DesignedOligo
    product_size: int
    pair_penalty: float | None
    metrics: tuple[tuple[str, str], ...]


@dataclass(frozen=True, slots=True)
class PrimerDesignResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    candidates: tuple[CandidateRegion, ...]
    assays: tuple[AssayCandidate, ...]
    candidate_regions_path: Path | None
    assays_path: Path | None
    primer3_input_path: Path | None
    primer3_output_path: Path | None
    report_path: Path


def _select_candidate_regions(
    conservation: ConservationResult, config: PrimerDesignConfig
) -> tuple[CandidateRegion, ...]:
    """Expand conservation windows into fixed-length reference intervals."""
    _validate_conservation_input(conservation)
    if not conservation.positions:
        return ()

    reference_length = len(conservation.positions)
    generated: list[tuple[CandidateRegion, WindowConservation]] = []
    for window in conservation.windows:
        reference_start, reference_end = _expanded_interval(
            window.reference_start,
            window.reference_end,
            reference_length,
            config.candidate_region_length,
        )
        region = _candidate_region(
            region_id="",
            rank=0,
            reference_start=reference_start,
            reference_end=reference_end,
            peak_start=window.reference_start,
            peak_end=window.reference_end,
            positions=conservation.positions[reference_start - 1:reference_end],
            config=config,
        )
        if _is_eligible(region, config):
            generated.append((region, window))

    unique_regions: dict[
        tuple[int, int], tuple[CandidateRegion, WindowConservation]
    ] = {}
    for region, window in generated:
        interval = region.reference_start, region.reference_end
        existing = unique_regions.get(interval)
        if existing is None or _window_ranking_key(window) < _window_ranking_key(
            existing[1]
        ):
            unique_regions[interval] = region, window

    selected: list[CandidateRegion] = []
    for region, _ in sorted(
        unique_regions.values(), key=lambda item: _ranking_key(item[0])
    ):
        if all(
            _overlap_fraction(region, accepted) <= config.max_region_overlap_fraction
            for accepted in selected
        ):
            selected.append(region)
            if len(selected) == config.max_candidate_regions:
                break

    return tuple(
        replace(region, region_id=f"region-{rank:03d}", rank=rank)
        for rank, region in enumerate(selected, 1)
    )


def _expanded_interval(
    peak_start: int,
    peak_end: int,
    reference_length: int,
    requested_length: int,
) -> tuple[int, int]:
    if reference_length <= requested_length:
        return 1, reference_length

    midpoint = (peak_start + peak_end) // 2
    reference_start = midpoint - (requested_length - 1) // 2
    reference_end = reference_start + requested_length - 1
    if reference_start < 1:
        return 1, requested_length
    if reference_end > reference_length:
        return reference_length - requested_length + 1, reference_length
    return reference_start, reference_end


def _candidate_region(
    *,
    region_id: str,
    rank: int,
    reference_start: int,
    reference_end: int,
    peak_start: int,
    peak_end: int,
    positions: tuple[PositionConservation, ...],
    config: PrimerDesignConfig,
) -> CandidateRegion:
    position_count = len(positions)
    usable_length = sum(
        position.major_allele_frequency >= config.min_minimum_conservation
        and position.coverage >= config.min_mean_coverage
        and position.gap_frequency <= config.max_mean_gap_frequency
        and position.entropy_bits <= config.max_mean_entropy_bits
        for position in positions
    )
    return CandidateRegion(
        region_id=region_id,
        rank=rank,
        reference_start=reference_start,
        reference_end=reference_end,
        peak_start=peak_start,
        peak_end=peak_end,
        position_count=position_count,
        usable_length=usable_length,
        usable_fraction=usable_length / position_count,
        mean_conservation=math.fsum(
            position.major_allele_frequency for position in positions
        ) / position_count,
        minimum_conservation=min(
            position.major_allele_frequency for position in positions
        ),
        mean_coverage=(
            math.fsum(position.coverage for position in positions) / position_count
        ),
        mean_gap_frequency=math.fsum(
            position.gap_frequency for position in positions
        ) / position_count,
        mean_entropy_bits=math.fsum(
            position.entropy_bits for position in positions
        ) / position_count,
    )


def _is_eligible(region: CandidateRegion, config: PrimerDesignConfig) -> bool:
    return (
        region.mean_conservation >= config.min_mean_conservation
        and region.minimum_conservation >= config.min_minimum_conservation
        and region.mean_coverage >= config.min_mean_coverage
        and region.mean_gap_frequency <= config.max_mean_gap_frequency
        and region.mean_entropy_bits <= config.max_mean_entropy_bits
        and region.usable_fraction >= config.min_usable_fraction
    )


def _ranking_key(region: CandidateRegion) -> tuple[float | int, ...]:
    return (
        -region.mean_conservation,
        -region.minimum_conservation,
        -region.mean_coverage,
        region.mean_entropy_bits,
        region.mean_gap_frequency,
        -region.usable_length,
        region.reference_start,
        region.reference_end,
    )


def _window_ranking_key(window: WindowConservation) -> tuple[float | int, ...]:
    return (
        -window.mean_conservation,
        -window.minimum_conservation,
        -window.mean_coverage,
        window.mean_entropy_bits,
        window.mean_gap_frequency,
        window.reference_start,
        window.reference_end,
    )


def _overlap_fraction(first: CandidateRegion, second: CandidateRegion) -> float:
    overlap_length = max(
        0,
        min(first.reference_end, second.reference_end)
        - max(first.reference_start, second.reference_start)
        + 1,
    )
    shorter_length = min(
        first.reference_end - first.reference_start + 1,
        second.reference_end - second.reference_start + 1,
    )
    return overlap_length / shorter_length


def _validate_conservation_input(conservation: ConservationResult) -> None:
    if conservation.status != "COMPLETE":
        raise PrimerDesignError("Candidate selection requires COMPLETE conservation.")

    if (
        not isinstance(conservation.reference_id, str)
        or not conservation.reference_id.strip()
    ):
        raise PrimerDesignError("Complete conservation requires a reference ID.")

    if not conservation.positions:
        if conservation.windows or conservation.major_consensus or conservation.iupac_consensus:
            raise PrimerDesignError("Empty conservation cannot contain positions or windows.")
        return

    reference_length = len(conservation.positions)
    if (
        len(conservation.major_consensus) != reference_length
        or len(conservation.iupac_consensus) != reference_length
    ):
        raise PrimerDesignError(
            "Conservation consensus length must match reference positions."
        )

    for expected_position, position in enumerate(conservation.positions, 1):
        if position.reference_position != expected_position:
            raise PrimerDesignError(
                "Conservation reference positions must be contiguous from 1."
            )

    for window in conservation.windows:
        if not (
            1 <= window.reference_start <= window.reference_end <= reference_length
        ):
            raise PrimerDesignError(
                "Conservation windows must lie within reference bounds."
            )
