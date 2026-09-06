from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)
from qpcr_pipeline.region_selection import (
    candidate_region_from_window,
    is_target_eligible,
    is_window_target_eligible,
    overlap_fraction,
    select_conservation_candidate_regions,
)


def _position(reference_position: int, **changes) -> PositionConservation:
    return replace(
        PositionConservation(
            alignment_position=reference_position,
            reference_position=reference_position,
            reference_base="A",
            depth=10,
            coverage=1.0,
            frequency_a=1.0,
            frequency_c=0.0,
            frequency_g=0.0,
            frequency_t=0.0,
            gap_frequency=0.0,
            major_allele_frequency=1.0,
            entropy_bits=0.0,
            major_consensus="A",
            iupac_consensus="A",
        ),
        **changes,
    )


def _positions(count: int) -> tuple[PositionConservation, ...]:
    return tuple(_position(index) for index in range(1, count + 1))


def _window(start: int, end: int | None = None, **changes) -> WindowConservation:
    end = start if end is None else end
    return replace(
        WindowConservation(
            reference_start=start,
            reference_end=end,
            position_count=end - start + 1,
            mean_conservation=1.0,
            minimum_conservation=1.0,
            mean_coverage=1.0,
            mean_gap_frequency=0.0,
            mean_entropy_bits=0.0,
        ),
        **changes,
    )


def _conservation(
    positions: tuple[PositionConservation, ...],
    windows: tuple[WindowConservation, ...],
) -> ConservationResult:
    consensus = "".join(item.major_consensus for item in positions)
    return ConservationResult(
        status="COMPLETE",
        reference_id="ref",
        positions=positions,
        windows=windows,
        annotations=(),
        major_consensus=consensus,
        iupac_consensus=consensus,
        position_metrics_path=None,
        window_metrics_path=None,
        major_consensus_path=None,
        iupac_consensus_path=None,
        html_report_path=None,
        report_path=Path("unused.json"),
    )


def _config(**changes) -> PrimerDesignConfig:
    return replace(
        PrimerDesignConfig(
            candidate_region_length=200,
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
        ),
        **changes,
    )


def test_boundary_intervals_and_ids_match_legacy_selection():
    result = _conservation(
        _positions(500),
        (_window(1, 100), _window(400, 500)),
    )
    regions = select_conservation_candidate_regions(
        result,
        _config(candidate_region_length=300),
    )
    assert [(item.reference_start, item.reference_end) for item in regions] == [
        (1, 300),
        (201, 500),
    ]
    assert [item.region_id for item in regions] == ["region-001", "region-002"]
    assert [item.rank for item in regions] == [1, 2]


def test_duplicate_interval_keeps_best_originating_window():
    result = _conservation(
        _positions(500),
        (
            _window(50, mean_conservation=0.8),
            _window(100, mean_conservation=0.9),
            _window(450),
        ),
    )
    regions = select_conservation_candidate_regions(
        result,
        _config(candidate_region_length=300),
    )
    assert [
        (item.reference_start, item.reference_end, item.peak_start)
        for item in regions
    ] == [(1, 300, 100), (201, 500, 450)]


def test_candidate_region_metrics_and_eligibility_are_shared():
    positions = list(_positions(4))
    positions[1] = _position(
        2,
        major_allele_frequency=0.9,
        coverage=0.8,
        gap_frequency=0.1,
        entropy_bits=0.5,
    )
    result = _conservation(tuple(positions), (_window(1, 4),))
    region = candidate_region_from_window(result, result.windows[0], _config())
    assert region.position_count == 4
    assert region.mean_conservation == 0.975
    assert region.minimum_conservation == 0.9
    assert region.mean_coverage == 0.95
    assert region.mean_gap_frequency == 0.025
    assert region.mean_entropy_bits == 0.125
    assert is_target_eligible(region, _config()) is True
    assert is_target_eligible(
        region,
        _config(min_mean_conservation=0.99),
    ) is False


def test_window_anchor_can_be_eligible_when_expanded_design_space_is_not():
    positions = list(_positions(300))
    for index in list(range(0, 100)) + list(range(200, 300)):
        positions[index] = _position(
            index + 1,
            major_allele_frequency=0.60,
            entropy_bits=0.80,
        )
    anchor = _window(
        101,
        200,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )
    result = _conservation(tuple(positions), (anchor,))
    config = _config(
        candidate_region_length=300,
        min_mean_conservation=0.90,
        min_minimum_conservation=0.70,
        min_mean_coverage=0.90,
        max_mean_gap_frequency=0.05,
        max_mean_entropy_bits=0.50,
        min_usable_fraction=0.80,
    )

    design_space = candidate_region_from_window(result, anchor, config)

    assert is_target_eligible(design_space, config) is False
    assert is_window_target_eligible(result, anchor, config) is True


def test_overlap_fraction_uses_shorter_region_and_threshold_geometry():
    result = _conservation(_positions(500), (_window(100), _window(250)))
    first, second = select_conservation_candidate_regions(
        result,
        _config(max_region_overlap_fraction=0.25),
    )
    assert (first.reference_start, first.reference_end) == (1, 200)
    assert (second.reference_start, second.reference_end) == (151, 350)
    assert overlap_fraction(first, second) == 0.25


def test_excess_overlap_is_suppressed_and_ranks_are_compact():
    result = _conservation(
        _positions(600),
        (_window(100), _window(200), _window(500)),
    )
    regions = select_conservation_candidate_regions(
        result,
        _config(max_region_overlap_fraction=0.25),
    )
    assert [item.reference_start for item in regions] == [1, 401]
    assert [item.region_id for item in regions] == ["region-001", "region-002"]
    assert [item.rank for item in regions] == [1, 2]
