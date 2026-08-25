import unittest
from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)
from qpcr_pipeline.primer_design import PrimerDesignError, _select_candidate_regions


def _position(reference_position: int) -> PositionConservation:
    return PositionConservation(
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
    )


def _conservation(
    positions: tuple[PositionConservation, ...],
    windows: tuple[WindowConservation, ...],
) -> ConservationResult:
    consensus = "".join(position.major_consensus for position in positions)
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
        report_path=Path("unused-conservation-report.json"),
    )


def _positions(
    start: int,
    count: int,
    *,
    major_allele_frequency: float = 1.0,
    coverage: float = 1.0,
    gap_frequency: float = 0.0,
    entropy_bits: float = 0.0,
) -> tuple[PositionConservation, ...]:
    return tuple(
        replace(
            _position(reference_position),
            major_allele_frequency=major_allele_frequency,
            coverage=coverage,
            gap_frequency=gap_frequency,
            entropy_bits=entropy_bits,
        )
        for reference_position in range(start, start + count)
    )


def _window(
    peak: int,
    *,
    mean_conservation: float = 1.0,
    minimum_conservation: float = 1.0,
    mean_coverage: float = 1.0,
    mean_gap_frequency: float = 0.0,
    mean_entropy_bits: float = 0.0,
) -> WindowConservation:
    return WindowConservation(
        peak,
        peak,
        1,
        mean_conservation,
        minimum_conservation,
        mean_coverage,
        mean_gap_frequency,
        mean_entropy_bits,
    )


def _permissive_config(**changes: float | int) -> PrimerDesignConfig:
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


class CandidateRegionSelectionTests(unittest.TestCase):
    def test_boundary_intervals_shift_to_preserve_configured_length(self):
        conservation = _conservation(
            tuple(_position(index) for index in range(1, 501)),
            (
                WindowConservation(1, 100, 100, 1.0, 1.0, 1.0, 0.0, 0.0),
                WindowConservation(400, 500, 101, 1.0, 1.0, 1.0, 0.0, 0.0),
            ),
        )
        permissive_config = PrimerDesignConfig(
            candidate_region_length=300,
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
        )

        regions = _select_candidate_regions(conservation, permissive_config)

        self.assertEqual(
            [(item.reference_start, item.reference_end) for item in regions],
            [(1, 300), (201, 500)],
        )

    def test_odd_length_interval_is_centered_on_an_integer_peak_midpoint(self):
        regions = _select_candidate_regions(
            _conservation(_positions(1, 500), (_window(250),)),
            _permissive_config(candidate_region_length=201),
        )

        self.assertEqual(
            [(item.reference_start, item.reference_end) for item in regions],
            [(150, 350)],
        )

    def test_aggregates_literal_metrics_and_usable_positions(self):
        conservation = _conservation(
            (
                _position(1),
                replace(_position(2), major_allele_frequency=0.9, coverage=0.8, gap_frequency=0.1, entropy_bits=0.5),
                replace(_position(3), major_allele_frequency=0.9, coverage=0.9, entropy_bits=0.5),
                replace(_position(4), major_allele_frequency=0.9, coverage=0.9),
            ),
            (WindowConservation(1, 4, 4, 0.925, 0.9, 0.9, 0.025, 0.25),),
        )

        regions = _select_candidate_regions(
            conservation, PrimerDesignConfig(min_usable_fraction=0.75)
        )

        self.assertEqual(len(regions), 1)
        region = regions[0]
        self.assertEqual(region.position_count, 4)
        self.assertEqual(region.usable_length, 3)
        self.assertEqual(region.usable_fraction, 0.75)
        self.assertAlmostEqual(region.mean_conservation, 0.925)
        self.assertEqual(region.minimum_conservation, 0.9)
        self.assertAlmostEqual(region.mean_coverage, 0.9)
        self.assertAlmostEqual(region.mean_gap_frequency, 0.025)
        self.assertAlmostEqual(region.mean_entropy_bits, 0.25)

    def test_rejects_each_independent_eligibility_failure(self):
        positions = tuple(_position(index) for index in range(1, 5))
        windows = (WindowConservation(1, 4, 4, 1.0, 1.0, 1.0, 0.0, 0.0),)
        base_config = PrimerDesignConfig(
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
        )
        cases = (
            (
                "mean conservation",
                replace(base_config, min_mean_conservation=0.95),
                replace(positions[0], major_allele_frequency=0.79),
            ),
            (
                "minimum conservation",
                replace(base_config, min_minimum_conservation=0.75),
                replace(positions[0], major_allele_frequency=0.74),
            ),
            (
                "mean coverage",
                replace(base_config, min_mean_coverage=0.95),
                replace(positions[0], coverage=0.79),
            ),
            (
                "mean gap frequency",
                replace(base_config, max_mean_gap_frequency=0.05),
                replace(positions[0], gap_frequency=0.21),
            ),
            (
                "mean entropy",
                replace(base_config, max_mean_entropy_bits=0.49),
                replace(positions[0], entropy_bits=2.0),
            ),
            (
                "usable fraction",
                replace(
                    base_config,
                    min_mean_coverage=0.75,
                    min_usable_fraction=0.8,
                ),
                replace(positions[0], coverage=0.74),
            ),
        )

        for name, config, mutated_position in cases:
            with self.subTest(threshold=name):
                conservation = _conservation(
                    (mutated_position, *positions[1:]), windows
                )

                regions = _select_candidate_regions(conservation, config)

                self.assertEqual(regions, ())

    def test_deduplicates_intervals_using_the_best_ranked_peak(self):
        conservation = _conservation(
            _positions(1, 500),
            (
                _window(50, mean_conservation=0.8),
                _window(100, mean_conservation=0.9),
                _window(450),
            ),
        )

        regions = _select_candidate_regions(
            conservation, _permissive_config(candidate_region_length=300)
        )

        self.assertEqual(
            [
                (item.reference_start, item.reference_end, item.peak_start, item.peak_end)
                for item in regions
            ],
            [(1, 300, 100, 100), (201, 500, 450, 450)],
        )
        self.assertEqual([item.region_id for item in regions], ["region-001", "region-002"])
        self.assertEqual([item.rank for item in regions], [1, 2])

    def test_ranks_by_each_reachable_lexicographic_key(self):
        base_config = _permissive_config()
        cases = (
            (
                "mean conservation",
                _positions(1, 200, major_allele_frequency=0.9)
                + _positions(201, 200, major_allele_frequency=0.95),
                (_window(100), _window(300)),
                base_config,
                [300, 100],
            ),
            (
                "minimum conservation",
                _positions(1, 200, major_allele_frequency=0.95)
                + _positions(201, 199, major_allele_frequency=0.951)
                + _positions(400, 1, major_allele_frequency=0.751),
                (_window(300), _window(100)),
                base_config,
                [100, 300],
            ),
            (
                "mean coverage",
                _positions(1, 200, coverage=0.9)
                + _positions(201, 200, coverage=0.95),
                (_window(100), _window(300)),
                base_config,
                [300, 100],
            ),
            (
                "mean entropy",
                _positions(1, 200, entropy_bits=0.3)
                + _positions(201, 200, entropy_bits=0.2),
                (_window(100), _window(300)),
                base_config,
                [300, 100],
            ),
            (
                "mean gap frequency",
                _positions(1, 200, gap_frequency=0.2)
                + _positions(201, 200, gap_frequency=0.1),
                (_window(100), _window(300)),
                base_config,
                [300, 100],
            ),
            (
                "usable length",
                _positions(1, 100, coverage=0.6)
                + _positions(101, 100, gap_frequency=0.4)
                + _positions(201, 100)
                + _positions(301, 100, coverage=0.6, gap_frequency=0.4),
                (_window(100), _window(300)),
                _permissive_config(
                    min_mean_coverage=0.8,
                    max_mean_gap_frequency=0.2,
                ),
                [300, 100],
            ),
            (
                "reference start",
                _positions(1, 400),
                (_window(300), _window(100)),
                base_config,
                [100, 300],
            ),
        )

        for name, positions, windows, config, expected_peaks in cases:
            with self.subTest(ordering_key=name):
                regions = _select_candidate_regions(
                    _conservation(positions, windows), config
                )

                self.assertEqual([item.peak_start for item in regions], expected_peaks)

    def test_allows_equal_overlap_rejects_excess_and_reassigns_selected_ranks(self):
        equal_overlap = _select_candidate_regions(
            _conservation(_positions(1, 500), (_window(100), _window(250))),
            _permissive_config(max_region_overlap_fraction=0.25),
        )
        excess_overlap = _select_candidate_regions(
            _conservation(_positions(1, 500), (_window(100), _window(200))),
            _permissive_config(max_region_overlap_fraction=0.25),
        )
        rejected_before_later_region = _select_candidate_regions(
            _conservation(
                _positions(1, 600),
                (_window(100), _window(200), _window(500)),
            ),
            _permissive_config(max_region_overlap_fraction=0.25),
        )

        self.assertEqual(
            [(item.reference_start, item.reference_end) for item in equal_overlap],
            [(1, 200), (151, 350)],
        )
        self.assertEqual(
            [(item.reference_start, item.reference_end) for item in excess_overlap],
            [(1, 200)],
        )
        self.assertEqual(
            [item.reference_start for item in rejected_before_later_region], [1, 401]
        )
        self.assertEqual(
            [item.region_id for item in rejected_before_later_region],
            ["region-001", "region-002"],
        )
        self.assertEqual([item.rank for item in rejected_before_later_region], [1, 2])

    def test_stops_after_configured_candidate_cap(self):
        regions = _select_candidate_regions(
            _conservation(
                _positions(1, 600),
                (_window(100), _window(300), _window(500)),
            ),
            _permissive_config(max_candidate_regions=2, max_region_overlap_fraction=0.0),
        )

        self.assertEqual([item.peak_start for item in regions], [100, 300])

    def test_empty_complete_conservation_returns_no_regions(self):
        self.assertEqual(
            _select_candidate_regions(
                _conservation((), ()), _permissive_config()
            ),
            (),
        )

    def test_rejects_malformed_conservation_before_selection(self):
        base = _conservation(
            _positions(1, 4), (WindowConservation(1, 4, 4, 1.0, 1.0, 1.0, 0.0, 0.0),)
        )
        empty = _conservation((), ())
        malformed_cases = (
            ("incomplete status", replace(base, status="SKIPPED"), "COMPLETE"),
            ("missing reference", replace(base, reference_id=None), "reference ID"),
            ("blank reference", replace(base, reference_id=" "), "reference ID"),
            (
                "empty missing reference",
                replace(empty, reference_id=None),
                "reference ID",
            ),
            (
                "empty blank reference",
                replace(empty, reference_id=" "),
                "reference ID",
            ),
            (
                "coordinate gap",
                replace(
                    base,
                    positions=base.positions[:2]
                    + (replace(base.positions[2], reference_position=4),)
                    + (replace(base.positions[3], reference_position=5),),
                ),
                "contiguous",
            ),
            (
                "duplicate coordinate",
                replace(
                    base,
                    positions=base.positions[:2]
                    + (replace(base.positions[2], reference_position=2),)
                    + base.positions[3:],
                ),
                "contiguous",
            ),
            (
                "consensus mismatch",
                replace(base, major_consensus="AAA"),
                "consensus length",
            ),
            (
                "out of bounds window",
                replace(base, windows=(_window(5),)),
                "within reference bounds",
            ),
        )

        for name, conservation, message in malformed_cases:
            with self.subTest(input=name), self.assertRaisesRegex(
                PrimerDesignError, message
            ):
                _select_candidate_regions(conservation, _permissive_config())
