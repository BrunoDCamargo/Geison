"""Target-conservation constraints for Primer3 binding sites."""

from __future__ import annotations

from collections.abc import Mapping

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.region_selection import CandidateRegion


def binding_site_exclusion_intervals(
    conservation: object,
    candidate: CandidateRegion,
    config: PrimerDesignConfig,
) -> tuple[tuple[int, int], ...]:
    """Return 1-based inclusive intervals unsuitable for primer/probe binding.

    Hard exclusions use only the existing absolute binding-site floors:
    minimum target conservation, coverage, and gap frequency. Entropy remains a
    region/anchor quality metric and is not converted into a per-base hard mask.
    """
    positions = tuple(getattr(conservation, "positions", ()))
    if not positions:
        return ()

    by_reference = {
        position.reference_position: position
        for position in positions
        if getattr(position, "reference_position", None) is not None
    }
    excluded: list[int] = []
    for reference_position in range(
        candidate.reference_start,
        candidate.reference_end + 1,
    ):
        position = by_reference.get(reference_position)
        if position is None:
            raise ValueError(
                "Primer3 conservation constraints require contiguous reference positions."
            )
        if (
            position.major_allele_frequency < config.min_minimum_conservation
            or position.coverage < config.min_mean_coverage
            or position.gap_frequency > config.max_mean_gap_frequency
        ):
            excluded.append(reference_position)

    return _compress_positions(excluded)


def binding_site_exclusions_by_candidate(
    conservation: object,
    candidates: tuple[CandidateRegion, ...],
    config: PrimerDesignConfig,
) -> dict[str, tuple[tuple[int, int], ...]]:
    return {
        candidate.region_id: binding_site_exclusion_intervals(
            conservation,
            candidate,
            config,
        )
        for candidate in candidates
    }


def inject_primer3_exclusions(
    input_text: str,
    exclusions_by_candidate: Mapping[str, tuple[tuple[int, int], ...]],
) -> str:
    """Insert Primer3 primer/probe exclusion tags into Boulder-IO records."""
    if not exclusions_by_candidate:
        return input_text

    lines = input_text.splitlines()
    output: list[str] = []
    current_id: str | None = None
    for line in lines:
        output.append(line)
        if line.startswith("SEQUENCE_ID="):
            current_id = line.removeprefix("SEQUENCE_ID=")
            continue
        if line.startswith("SEQUENCE_INCLUDED_REGION=") and current_id is not None:
            intervals = exclusions_by_candidate.get(current_id, ())
            if intervals:
                serialized = " ".join(
                    f"{start - 1},{end - start + 1}"
                    for start, end in intervals
                )
                output.append(f"SEQUENCE_EXCLUDED_REGION={serialized}")
                output.append(f"SEQUENCE_INTERNAL_EXCLUDED_REGION={serialized}")
        if line == "=":
            current_id = None

    return "\n".join(output) + ("\n" if input_text.endswith("\n") else "")


def _compress_positions(positions: list[int]) -> tuple[tuple[int, int], ...]:
    if not positions:
        return ()

    intervals: list[tuple[int, int]] = []
    start = previous = positions[0]
    for position in positions[1:]:
        if position == previous + 1:
            previous = position
            continue
        intervals.append((start, previous))
        start = previous = position
    intervals.append((start, previous))
    return tuple(intervals)
