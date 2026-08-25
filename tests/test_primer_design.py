import json
import os
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)
from qpcr_pipeline.primer_design import (
    PrimerDesignError,
    _select_candidate_regions,
    design_primers,
)


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


_CANDIDATE_HEADER = (
    "region_id\trank\treference_start\treference_end\tpeak_start\tpeak_end\t"
    "position_count\tusable_length\tusable_fraction\tmean_conservation\t"
    "minimum_conservation\tmean_coverage\tmean_gap_frequency\tmean_entropy_bits\n"
)
_ASSAY_HEADER = (
    "assay_id\tregion_id\tprimer3_index\tforward_sequence\t"
    "forward_reference_start\tforward_reference_end\tforward_length\tforward_tm\t"
    "forward_gc_percent\tforward_penalty\tprobe_sequence\tprobe_reference_start\t"
    "probe_reference_end\tprobe_length\tprobe_tm\tprobe_gc_percent\tprobe_penalty\t"
    "reverse_sequence\treverse_reference_start\treverse_reference_end\treverse_length\t"
    "reverse_tm\treverse_gc_percent\treverse_penalty\tproduct_size\tpair_penalty\n"
)
_COMPLETE_PRIMER3_OUTPUT = (
    "SEQUENCE_ID=region-001\n"
    "PRIMER_WARNING=fixture warning\n"
    "PRIMER_LEFT_EXPLAIN=considered 1, ok 1\n"
    "PRIMER_INTERNAL_EXPLAIN=considered 1, ok 1\n"
    "PRIMER_RIGHT_EXPLAIN=considered 1, ok 1\n"
    "PRIMER_PAIR_EXPLAIN=considered 1, ok 1\n"
    "PRIMER_LEFT_NUM_RETURNED=1\n"
    "PRIMER_INTERNAL_NUM_RETURNED=1\n"
    "PRIMER_RIGHT_NUM_RETURNED=1\n"
    "PRIMER_PAIR_NUM_RETURNED=1\n"
    "PRIMER_LEFT_0=10,20\n"
    "PRIMER_LEFT_0_SEQUENCE=AAAAAAAAAAAAAAAAAAAA\n"
    "PRIMER_LEFT_0_TM=60.25\n"
    "PRIMER_LEFT_0_GC_PERCENT=42.5\n"
    "PRIMER_LEFT_0_PENALTY=0.5\n"
    "PRIMER_LEFT_0_SELF_ANY_TH=1.2\n"
    "PRIMER_INTERNAL_0=40,20\n"
    "PRIMER_INTERNAL_0_SEQUENCE=CCCCCCCCCCCCCCCCCCCC\n"
    "PRIMER_INTERNAL_0_TM=70.5\n"
    "PRIMER_INTERNAL_0_GC_PERCENT=55.0\n"
    "PRIMER_INTERNAL_0_PENALTY=0.25\n"
    "PRIMER_INTERNAL_0_HAIRPIN_TH=0.75\n"
    "PRIMER_RIGHT_0=109,20\n"
    "PRIMER_RIGHT_0_SEQUENCE=GGGGGGGGGGGGGGGGGGGG\n"
    "PRIMER_RIGHT_0_TM=61.0\n"
    "PRIMER_RIGHT_0_GC_PERCENT=47.5\n"
    "PRIMER_RIGHT_0_PENALTY=0.75\n"
    "PRIMER_RIGHT_0_SELF_END_TH=0.8\n"
    "PRIMER_PAIR_0_PRODUCT_SIZE=100\n"
    "PRIMER_PAIR_0_PENALTY=1.5\n"
    "PRIMER_PAIR_0_COMPL_ANY_TH=2.75\n"
    "=\n"
)
_ZERO_PAIR_PRIMER3_OUTPUT = (
    "SEQUENCE_ID=region-001\n"
    "PRIMER_WARNING=no viable pair\n"
    "PRIMER_LEFT_EXPLAIN=considered 1, low tm 1\n"
    "PRIMER_INTERNAL_EXPLAIN=considered 1, high tm 1\n"
    "PRIMER_RIGHT_EXPLAIN=considered 1, high tm 1\n"
    "PRIMER_PAIR_EXPLAIN=considered 0\n"
    "PRIMER_LEFT_NUM_RETURNED=0\n"
    "PRIMER_INTERNAL_NUM_RETURNED=0\n"
    "PRIMER_RIGHT_NUM_RETURNED=0\n"
    "PRIMER_PAIR_NUM_RETURNED=0\n"
    "=\n"
)
_PARSER_ERROR_OUTPUT = (
    "SEQUENCE_ID=region-001\n"
    "PRIMER_ERROR=literal parser failure\n"
    "=\n"
)


class _LiteralPrimer3Runner:
    def __init__(self, response: str):
        self.response = response
        self.inputs: list[str] = []

    def run(self, input_text: str) -> str:
        self.inputs.append(input_text)
        return self.response


class _RaisingPrimer3Runner:
    def __init__(self, error: Exception):
        self.error = error

    def run(self, input_text: str) -> str:
        del input_text
        raise self.error


class _UnreadableConservation:
    def __getattribute__(self, name: str):
        raise AssertionError(f"disabled design read conservation attribute {name!r}")


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


class PrimerDesignServiceTests(unittest.TestCase):
    def test_enabled_design_publishes_typed_auditable_artifacts(self):
        conservation = _conservation(
            _positions(1, 300),
            (WindowConservation(1, 300, 300, 1.0, 1.0, 1.0, 0.0, 0.0),),
        )
        config = _permissive_config(enabled=True, candidate_region_length=300)
        runner = _LiteralPrimer3Runner(_COMPLETE_PRIMER3_OUTPUT)

        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)

            result = design_primers(
                conservation, config, output_dir, runner=runner
            )

            artifact_dir = output_dir / "primer_design"
            expected_paths = {
                "candidate_regions": artifact_dir / "candidate_regions.tsv",
                "assays": artifact_dir / "assays.tsv",
                "primer3_input": artifact_dir / "primer3_input.txt",
                "primer3_output": artifact_dir / "primer3_output.txt",
                "report": artifact_dir / "primer_design_report.json",
            }
            self.assertEqual(result.status, "COMPLETE")
            self.assertEqual(result.reference_id, "ref")
            self.assertEqual(len(result.candidates), 1)
            self.assertEqual(len(result.assays), 1)
            self.assertEqual(result.candidate_regions_path, expected_paths["candidate_regions"])
            self.assertEqual(result.assays_path, expected_paths["assays"])
            self.assertEqual(result.primer3_input_path, expected_paths["primer3_input"])
            self.assertEqual(result.primer3_output_path, expected_paths["primer3_output"])
            self.assertEqual(result.report_path, expected_paths["report"])

            candidate_text = expected_paths["candidate_regions"].read_text(encoding="utf-8")
            assay_text = expected_paths["assays"].read_text(encoding="utf-8")
            self.assertEqual(
                candidate_text,
                _CANDIDATE_HEADER
                + "region-001\t1\t1\t300\t1\t300\t300\t300\t1.0\t1.0\t1.0\t1.0\t0.0\t0.0\n",
            )
            self.assertEqual(
                assay_text,
                _ASSAY_HEADER
                + "region-001-assay-001\tregion-001\t0\tAAAAAAAAAAAAAAAAAAAA\t"
                "11\t30\t20\t60.25\t42.5\t0.5\tCCCCCCCCCCCCCCCCCCCC\t"
                "41\t60\t20\t70.5\t55.0\t0.25\tGGGGGGGGGGGGGGGGGGGG\t"
                "91\t110\t20\t61.0\t47.5\t0.75\t100\t1.5\n",
            )
            self.assertEqual(
                expected_paths["primer3_input"].read_text(encoding="utf-8"),
                runner.inputs[0],
            )
            self.assertEqual(
                expected_paths["primer3_output"].read_text(encoding="utf-8"),
                _COMPLETE_PRIMER3_OUTPUT,
            )

            report_text = expected_paths["report"].read_text(encoding="utf-8")
            report = json.loads(report_text)
            self.assertEqual(report["schema_version"], 1)
            self.assertEqual(report["status"], "COMPLETE")
            self.assertIs(report["enabled"], True)
            self.assertEqual(report["reference_id"], "ref")
            self.assertEqual(report["counts"], {"candidates": 1, "assays": 1})
            self.assertIsInstance(report["configuration"]["assays_per_region"], int)
            self.assertIsInstance(report["configuration"]["min_mean_conservation"], float)
            self.assertEqual(report["configuration"]["min_mean_conservation"], 0.0)
            self.assertEqual(report["configuration"]["max_mean_entropy_bits"], 2.0)
            self.assertEqual(
                report["candidates"],
                [
                    {
                        "region_id": "region-001",
                        "rank": 1,
                        "reference_start": 1,
                        "reference_end": 300,
                        "peak_start": 1,
                        "peak_end": 300,
                        "position_count": 300,
                        "usable_length": 300,
                        "usable_fraction": 1.0,
                        "mean_conservation": 1.0,
                        "minimum_conservation": 1.0,
                        "mean_coverage": 1.0,
                        "mean_gap_frequency": 0.0,
                        "mean_entropy_bits": 0.0,
                    }
                ],
            )
            self.assertIsInstance(report["candidates"][0]["rank"], int)
            self.assertIsInstance(report["candidates"][0]["mean_conservation"], float)
            assay = report["assays"][0]
            self.assertIsInstance(assay["primer3_index"], int)
            self.assertIsInstance(assay["product_size"], int)
            self.assertIsInstance(assay["pair_penalty"], float)
            self.assertEqual(
                assay,
                {
                    "assay_id": "region-001-assay-001",
                    "region_id": "region-001",
                    "primer3_index": 0,
                    "forward_primer": {
                        "sequence": "AAAAAAAAAAAAAAAAAAAA",
                        "reference_start": 11,
                        "reference_end": 30,
                        "length": 20,
                        "tm": 60.25,
                        "gc_percent": 42.5,
                        "penalty": 0.5,
                        "metrics": {"PRIMER_LEFT_0_SELF_ANY_TH": "1.2"},
                    },
                    "probe": {
                        "sequence": "CCCCCCCCCCCCCCCCCCCC",
                        "reference_start": 41,
                        "reference_end": 60,
                        "length": 20,
                        "tm": 70.5,
                        "gc_percent": 55.0,
                        "penalty": 0.25,
                        "metrics": {"PRIMER_INTERNAL_0_HAIRPIN_TH": "0.75"},
                    },
                    "reverse_primer": {
                        "sequence": "GGGGGGGGGGGGGGGGGGGG",
                        "reference_start": 91,
                        "reference_end": 110,
                        "length": 20,
                        "tm": 61.0,
                        "gc_percent": 47.5,
                        "penalty": 0.75,
                        "metrics": {"PRIMER_RIGHT_0_SELF_END_TH": "0.8"},
                    },
                    "product_size": 100,
                    "pair_penalty": 1.5,
                    "metrics": {"PRIMER_PAIR_0_COMPL_ANY_TH": "2.75"},
                },
            )
            self.assertEqual(
                report["primer3_details"]["region-001"]["PRIMER_WARNING"],
                "fixture warning",
            )
            self.assertEqual(
                report["artifacts"],
                {
                    "candidate_regions": "primer_design/candidate_regions.tsv",
                    "assays": "primer_design/assays.tsv",
                    "primer3_input": "primer_design/primer3_input.txt",
                    "primer3_output": "primer_design/primer3_output.txt",
                },
            )
            self.assertNotIn("SEQUENCE_TEMPLATE", report_text)

    def test_disabled_design_removes_only_stale_artifacts_and_publishes_skipped(self):
        config = PrimerDesignConfig(enabled=False)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            artifact_dir = output_dir / "primer_design"
            artifact_dir.mkdir()
            stale_names = (
                "candidate_regions.tsv",
                "assays.tsv",
                "primer3_input.txt",
                "primer3_output.txt",
            )
            for name in stale_names:
                (artifact_dir / name).write_text("stale", encoding="utf-8")
            report_path = artifact_dir / "primer_design_report.json"
            report_path.write_text('{"status":"COMPLETE"}\n', encoding="utf-8")
            sibling = artifact_dir / "keep-me.txt"
            sibling.write_text("preserved", encoding="utf-8")

            with patch(
                "qpcr_pipeline.primer3.SubprocessPrimer3Runner",
                side_effect=AssertionError("runner constructed"),
            ):
                result = design_primers(_UnreadableConservation(), config, output_dir)

            self.assertEqual(result.status, "SKIPPED")
            self.assertIsNone(result.reference_id)
            self.assertEqual(result.candidates, ())
            self.assertEqual(result.assays, ())
            self.assertIsNone(result.candidate_regions_path)
            self.assertIsNone(result.assays_path)
            self.assertIsNone(result.primer3_input_path)
            self.assertIsNone(result.primer3_output_path)
            self.assertEqual(result.report_path, report_path)
            for name in stale_names:
                self.assertFalse((artifact_dir / name).exists(), name)
            self.assertEqual(sibling.read_text(encoding="utf-8"), "preserved")
            self.assertEqual(
                json.loads(report_path.read_text(encoding="utf-8")),
                {
                    "schema_version": 1,
                    "status": "SKIPPED",
                    "enabled": False,
                    "configuration": {
                        "enabled": False,
                        "max_candidate_regions": 10,
                        "assays_per_region": 5,
                        "candidate_region_length": 300,
                        "max_region_overlap_fraction": 0.5,
                        "min_mean_conservation": 0.9,
                        "min_minimum_conservation": 0.7,
                        "min_mean_coverage": 0.9,
                        "max_mean_gap_frequency": 0.05,
                        "max_mean_entropy_bits": 0.5,
                        "min_usable_fraction": 0.8,
                        "product_size_min": 70,
                        "product_size_max": 200,
                        "primer": {
                            "min_size": 18,
                            "opt_size": 20,
                            "max_size": 25,
                            "min_tm": 58.0,
                            "opt_tm": 60.0,
                            "max_tm": 62.0,
                            "min_gc_percent": 40.0,
                            "max_gc_percent": 60.0,
                        },
                        "probe": {
                            "min_size": 18,
                            "opt_size": 25,
                            "max_size": 30,
                            "min_tm": 68.0,
                            "opt_tm": 70.0,
                            "max_tm": 72.0,
                            "min_gc_percent": 30.0,
                            "max_gc_percent": 80.0,
                        },
                    },
                    "reference_id": None,
                    "counts": {"candidates": 0, "assays": 0},
                    "candidates": [],
                    "assays": [],
                    "primer3_details": {},
                    "artifacts": {
                        "candidate_regions": None,
                        "assays": None,
                        "primer3_input": None,
                        "primer3_output": None,
                    },
                },
            )

    def test_enabled_design_without_candidates_publishes_header_only_tsvs(self):
        conservation = _conservation((), ())
        config = _permissive_config(enabled=True)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            artifact_dir = output_dir / "primer_design"
            artifact_dir.mkdir()
            (artifact_dir / "primer3_input.txt").write_text("stale", encoding="utf-8")
            (artifact_dir / "primer3_output.txt").write_text("stale", encoding="utf-8")

            with patch(
                "qpcr_pipeline.primer3.SubprocessPrimer3Runner",
                side_effect=AssertionError("runner constructed"),
            ):
                result = design_primers(conservation, config, output_dir)

            self.assertEqual(result.status, "COMPLETE")
            self.assertEqual(result.reference_id, "ref")
            self.assertEqual(result.candidates, ())
            self.assertEqual(result.assays, ())
            self.assertEqual(
                result.candidate_regions_path, artifact_dir / "candidate_regions.tsv"
            )
            self.assertEqual(result.assays_path, artifact_dir / "assays.tsv")
            self.assertIsNone(result.primer3_input_path)
            self.assertIsNone(result.primer3_output_path)
            self.assertEqual(
                result.candidate_regions_path.read_text(encoding="utf-8"),
                _CANDIDATE_HEADER,
            )
            self.assertEqual(
                result.assays_path.read_text(encoding="utf-8"), _ASSAY_HEADER
            )
            self.assertFalse((artifact_dir / "primer3_input.txt").exists())
            self.assertFalse((artifact_dir / "primer3_output.txt").exists())
            report = json.loads(result.report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "COMPLETE")
            self.assertEqual(report["counts"], {"candidates": 0, "assays": 0})
            self.assertEqual(
                report["artifacts"],
                {
                    "candidate_regions": "primer_design/candidate_regions.tsv",
                    "assays": "primer_design/assays.tsv",
                    "primer3_input": None,
                    "primer3_output": None,
                },
            )

    def test_candidates_with_zero_pairs_publish_raw_exchange_and_explanations(self):
        conservation = _conservation(
            _positions(1, 300),
            (WindowConservation(1, 300, 300, 1.0, 1.0, 1.0, 0.0, 0.0),),
        )
        config = _permissive_config(enabled=True, candidate_region_length=300)
        runner = _LiteralPrimer3Runner(_ZERO_PAIR_PRIMER3_OUTPUT)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)

            result = design_primers(
                conservation, config, output_dir, runner=runner
            )

            self.assertEqual(len(result.candidates), 1)
            self.assertEqual(result.assays, ())
            self.assertEqual(result.assays_path.read_text(encoding="utf-8"), _ASSAY_HEADER)
            self.assertEqual(
                result.primer3_input_path.read_text(encoding="utf-8"),
                runner.inputs[0],
            )
            self.assertEqual(
                result.primer3_output_path.read_text(encoding="utf-8"),
                _ZERO_PAIR_PRIMER3_OUTPUT,
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["counts"], {"candidates": 1, "assays": 0})
            self.assertEqual(
                report["primer3_details"]["region-001"],
                {
                    "PRIMER_WARNING": "no viable pair",
                    "PRIMER_LEFT_EXPLAIN": "considered 1, low tm 1",
                    "PRIMER_INTERNAL_EXPLAIN": "considered 1, high tm 1",
                    "PRIMER_RIGHT_EXPLAIN": "considered 1, high tm 1",
                    "PRIMER_PAIR_EXPLAIN": "considered 0",
                },
            )

    def test_runner_and_parser_errors_invalidate_prior_complete_report(self):
        conservation = _conservation(
            _positions(1, 300),
            (WindowConservation(1, 300, 300, 1.0, 1.0, 1.0, 0.0, 0.0),),
        )
        config = _permissive_config(enabled=True, candidate_region_length=300)
        cases = (
            (
                "runner error",
                _RaisingPrimer3Runner(PrimerDesignError("bounded runner failure")),
                "bounded runner failure",
            ),
            (
                "parser error",
                _LiteralPrimer3Runner(_PARSER_ERROR_OUTPUT),
                "literal parser failure",
            ),
        )
        for name, runner, message in cases:
            with self.subTest(failure=name), tempfile.TemporaryDirectory() as temporary_directory:
                output_dir = Path(temporary_directory)
                artifact_dir = output_dir / "primer_design"
                artifact_dir.mkdir()
                report_path = artifact_dir / "primer_design_report.json"
                report_path.write_text('{"status":"COMPLETE"}\n', encoding="utf-8")

                with self.assertRaisesRegex(PrimerDesignError, message):
                    design_primers(conservation, config, output_dir, runner=runner)

                self.assertFalse(report_path.exists())
                self.assertFalse((artifact_dir / "primer3_input.txt").exists())
                self.assertFalse((artifact_dir / "primer3_output.txt").exists())

    def test_each_enabled_replacement_is_atomic_report_last_and_hardlink_safe(self):
        conservation = _conservation(
            _positions(1, 300),
            (WindowConservation(1, 300, 300, 1.0, 1.0, 1.0, 0.0, 0.0),),
        )
        config = _permissive_config(enabled=True, candidate_region_length=300)
        publication_order = (
            "candidate_regions.tsv",
            "assays.tsv",
            "primer3_input.txt",
            "primer3_output.txt",
            "primer_design_report.json",
        )
        original_replace = Path.replace

        for failed_name in publication_order:
            with self.subTest(destination=failed_name), tempfile.TemporaryDirectory() as temporary_directory:
                output_dir = Path(temporary_directory)
                artifact_dir = output_dir / "primer_design"
                artifact_dir.mkdir()
                destination = artifact_dir / failed_name
                frozen_source = output_dir / f"frozen-{failed_name}"
                frozen_text = (
                    '{"status":"COMPLETE"}\n'
                    if failed_name == "primer_design_report.json"
                    else f"frozen {failed_name}\n"
                )
                frozen_source.write_text(frozen_text, encoding="utf-8")
                os.link(frozen_source, destination)
                report_path = artifact_dir / "primer_design_report.json"
                if failed_name != "primer_design_report.json":
                    report_path.write_text(
                        '{"status":"COMPLETE"}\n', encoding="utf-8"
                    )
                sibling = artifact_dir / "keep-me.txt"
                sibling.write_text("preserved", encoding="utf-8")
                attempted_destinations: list[str] = []

                def fail_selected_replace(source: Path, target: Path) -> Path:
                    target_path = Path(target)
                    attempted_destinations.append(target_path.name)
                    if target_path.name == failed_name:
                        raise OSError(f"replace failed for {failed_name}")
                    return original_replace(source, target_path)

                with patch.object(
                    Path,
                    "replace",
                    autospec=True,
                    side_effect=fail_selected_replace,
                ), self.assertRaisesRegex(OSError, f"replace failed for {failed_name}"):
                    design_primers(
                        conservation,
                        config,
                        output_dir,
                        runner=_LiteralPrimer3Runner(_COMPLETE_PRIMER3_OUTPUT),
                    )

                failed_index = publication_order.index(failed_name)
                self.assertEqual(
                    attempted_destinations, list(publication_order[: failed_index + 1])
                )
                self.assertFalse(report_path.exists())
                self.assertEqual(frozen_source.read_text(encoding="utf-8"), frozen_text)
                self.assertEqual(sibling.read_text(encoding="utf-8"), "preserved")
                self.assertEqual(list(artifact_dir.glob(".*.tmp")), [])

    def test_disabled_cleanup_precedes_atomic_report_publication(self):
        config = PrimerDesignConfig(enabled=False)
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            artifact_dir = output_dir / "primer_design"
            artifact_dir.mkdir()
            stale_paths = tuple(
                artifact_dir / name
                for name in (
                    "candidate_regions.tsv",
                    "assays.tsv",
                    "primer3_input.txt",
                    "primer3_output.txt",
                )
            )
            for path in stale_paths:
                path.write_text("stale", encoding="utf-8")
            frozen_report = output_dir / "frozen-report.json"
            frozen_report.write_text('{"status":"COMPLETE"}\n', encoding="utf-8")
            report_path = artifact_dir / "primer_design_report.json"
            os.link(frozen_report, report_path)
            sibling = artifact_dir / "keep-me.txt"
            sibling.write_text("preserved", encoding="utf-8")
            attempted_destinations: list[str] = []

            def fail_report_replace(source: Path, target: Path) -> Path:
                del source
                attempted_destinations.append(Path(target).name)
                raise OSError("disabled report replace failed")

            with patch.object(
                Path,
                "replace",
                autospec=True,
                side_effect=fail_report_replace,
            ), self.assertRaisesRegex(OSError, "disabled report replace failed"):
                design_primers(_UnreadableConservation(), config, output_dir)

            self.assertEqual(attempted_destinations, ["primer_design_report.json"])
            self.assertFalse(report_path.exists())
            self.assertTrue(all(not path.exists() for path in stale_paths))
            self.assertEqual(
                frozen_report.read_text(encoding="utf-8"),
                '{"status":"COMPLETE"}\n',
            )
            self.assertEqual(sibling.read_text(encoding="utf-8"), "preserved")
            self.assertEqual(list(artifact_dir.glob(".*.tmp")), [])
