"""Shared deterministic target-region selection primitives."""

from __future__ import annotations

import math
from dataclasses import dataclass, replace

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)


class RegionSelectionError(RuntimeError):
    """Raised when conservation cannot produce traceable candidate regions."""


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


def candidate_region_from_window(
    conservation: ConservationResult,
    window: WindowConservation,
    config: PrimerDesignConfig,
) -> CandidateRegion:
    """Expand one conservation window and calculate target-side metrics."""
    reference_positions = _reference_positions(conservation)
    _validate_conservation_input(conservation, reference_positions)
    if not reference_positions:
        raise RegionSelectionError("Candidate selection requires reference positions.")
    if not (
        1
        <= window.reference_start
        <= window.reference_end
        <= len(reference_positions)
    ):
        raise RegionSelectionError(
            "Conservation windows must lie within reference bounds."
        )
    reference_start, reference_end = _expanded_interval(
        window.reference_start,
        window.reference_end,
        len(reference_positions),
        config.candidate_region_length,
    )
    return _candidate_region(
        region_id="",
        rank=0,
        reference_start=reference_start,
        reference_end=reference_end,
        peak_start=window.reference_start,
        peak_end=window.reference_end,
        positions=reference_positions[reference_start - 1 : reference_end],
        config=config,
    )


def is_target_eligible(region: CandidateRegion, config: PrimerDesignConfig) -> bool:
    """Return whether a complete design region satisfies the target-side policy."""
    return (
        region.mean_conservation >= config.min_mean_conservation
        and region.minimum_conservation >= config.min_minimum_conservation
        and region.mean_coverage >= config.min_mean_coverage
        and region.mean_gap_frequency <= config.max_mean_gap_frequency
        and region.mean_entropy_bits <= config.max_mean_entropy_bits
        and region.usable_fraction >= config.min_usable_fraction
    )


def is_window_target_eligible(
    conservation: ConservationResult,
    window: WindowConservation,
    config: PrimerDesignConfig,
) -> bool:
    """Return whether the conserved anchor window satisfies target-side policy.

    The window is the scientific anchor. ``candidate_region_length`` may expand
    around it to give Primer3 design space, but variability in those flanks must
    not retroactively make an otherwise eligible anchor ineligible.
    """
    reference_positions = _reference_positions(conservation)
    _validate_conservation_input(conservation, reference_positions)
    window_positions = reference_positions[
        window.reference_start - 1 : window.reference_end
    ]
    if not window_positions:
        return False
    usable_fraction = _usable_length(window_positions, config) / len(window_positions)
    return (
        window.mean_conservation >= config.min_mean_conservation
        and window.minimum_conservation >= config.min_minimum_conservation
        and window.mean_coverage >= config.min_mean_coverage
        and window.mean_gap_frequency <= config.max_mean_gap_frequency
        and window.mean_entropy_bits <= config.max_mean_entropy_bits
        and usable_fraction >= config.min_usable_fraction
    )


def overlap_fraction(first: CandidateRegion, second: CandidateRegion) -> float:
    """Return overlap length divided by the shorter region length."""
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


def select_conservation_candidate_regions(
    conservation: ConservationResult,
    config: PrimerDesignConfig,
) -> tuple[CandidateRegion, ...]:
    """Preserve the legacy conservation-only candidate selection semantics."""
    reference_positions = _reference_positions(conservation)
    _validate_conservation_input(conservation, reference_positions)
    if not reference_positions:
        return ()

    generated: list[tuple[CandidateRegion, WindowConservation]] = []
    for window in conservation.windows:
        region = candidate_region_from_window(conservation, window, config)
        if is_target_eligible(region, config):
            generated.append((region, window))

    unique_regions: dict[
        tuple[int, int], tuple[CandidateRegion, WindowConservation]
    ] = {}
    for region, window in generated:
        interval = region.reference_start, region.reference_end
        existing = unique_regions.get(interval)
        if existing is None or _window_ranking_key(
            window, reference_positions, config
        ) < _window_ranking_key(existing[1], reference_positions, config):
            unique_regions[interval] = region, window

    selected: list[CandidateRegion] = []
    for region, _ in sorted(
        unique_regions.values(), key=lambda item: _ranking_key(item[0])
    ):
        if all(
            overlap_fraction(region, accepted) <= config.max_region_overlap_fraction
            for accepted in selected
        ):
            selected.append(region)
            if len(selected) == config.max_candidate_regions:
                break

    return tuple(
        replace(region, region_id=f"region-{rank:03d}", rank=rank)
        for rank, region in enumerate(selected, 1)
    )


def _reference_positions(
    conservation: ConservationResult,
) -> tuple[PositionConservation, ...]:
    return tuple(
        position
        for position in conservation.positions
        if position.reference_position is not None
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
    usable_length = _usable_length(positions, config)
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
        )
        / position_count,
        minimum_conservation=min(
            position.major_allele_frequency for position in positions
        ),
        mean_coverage=math.fsum(position.coverage for position in positions)
        / position_count,
        mean_gap_frequency=math.fsum(
            position.gap_frequency for position in positions
        )
        / position_count,
        mean_entropy_bits=math.fsum(
            position.entropy_bits for position in positions
        )
        / position_count,
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


def _window_ranking_key(
    window: WindowConservation,
    reference_positions: tuple[PositionConservation, ...],
    config: PrimerDesignConfig,
) -> tuple[float | int, ...]:
    window_positions = reference_positions[
        window.reference_start - 1 : window.reference_end
    ]
    return (
        -window.mean_conservation,
        -window.minimum_conservation,
        -window.mean_coverage,
        window.mean_entropy_bits,
        window.mean_gap_frequency,
        -_usable_length(window_positions, config),
        window.reference_start,
        window.reference_end,
    )


def _usable_length(
    positions: tuple[PositionConservation, ...],
    config: PrimerDesignConfig,
) -> int:
    return sum(
        position.major_allele_frequency >= config.min_minimum_conservation
        and position.coverage >= config.min_mean_coverage
        and position.gap_frequency <= config.max_mean_gap_frequency
        and position.entropy_bits <= config.max_mean_entropy_bits
        for position in positions
    )


def _validate_conservation_input(
    conservation: ConservationResult,
    reference_positions: tuple[PositionConservation, ...],
) -> None:
    if conservation.status != "COMPLETE":
        raise RegionSelectionError(
            "Candidate selection requires COMPLETE conservation."
        )
    if (
        not isinstance(conservation.reference_id, str)
        or not conservation.reference_id.strip()
    ):
        raise RegionSelectionError("Complete conservation requires a reference ID.")

    if not conservation.positions:
        if (
            conservation.windows
            or conservation.major_consensus
            or conservation.iupac_consensus
        ):
            raise RegionSelectionError(
                "Empty conservation cannot contain positions or windows."
            )
        return

    reference_length = len(reference_positions)
    if (
        len(conservation.major_consensus) != reference_length
        or len(conservation.iupac_consensus) != reference_length
    ):
        raise RegionSelectionError(
            "Conservation consensus length must match reference positions."
        )

    for expected_position, position in enumerate(reference_positions, 1):
        if position.reference_position != expected_position:
            raise RegionSelectionError(
                "Conservation reference positions must be contiguous from 1."
            )

    for window in conservation.windows:
        if not (
            1
            <= window.reference_start
            <= window.reference_end
            <= reference_length
        ):
            raise RegionSelectionError(
                "Conservation windows must lie within reference bounds."
            )
