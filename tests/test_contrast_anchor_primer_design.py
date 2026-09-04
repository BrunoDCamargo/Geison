from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.primer_design import CandidateRegion, design_primers


class LiteralRunner:
    def __init__(self, response: str) -> None:
        self.response = response
        self.inputs: list[str] = []

    def run(self, input_text: str) -> str:
        self.inputs.append(input_text)
        return self.response


def candidate() -> CandidateRegion:
    return CandidateRegion(
        region_id="contrast-region-001",
        rank=1,
        reference_start=501,
        reference_end=800,
        peak_start=601,
        peak_end=700,
        position_count=300,
        usable_length=300,
        usable_fraction=1.0,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def anchored_output() -> str:
    return (
        "SEQUENCE_ID=contrast-region-001\n"
        "PRIMER_LEFT_NUM_RETURNED=1\n"
        "PRIMER_INTERNAL_NUM_RETURNED=1\n"
        "PRIMER_RIGHT_NUM_RETURNED=1\n"
        "PRIMER_PAIR_NUM_RETURNED=1\n"
        "PRIMER_LEFT_0=579,20\n"
        "PRIMER_LEFT_0_SEQUENCE=AAAAAAAAAAAAAAAAAAAA\n"
        "PRIMER_LEFT_0_TM=60.0\n"
        "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
        "PRIMER_INTERNAL_0=619,25\n"
        "PRIMER_INTERNAL_0_SEQUENCE=CCCCCCCCCCCCCCCCCCCCCCCCC\n"
        "PRIMER_INTERNAL_0_TM=70.0\n"
        "PRIMER_INTERNAL_0_GC_PERCENT=50.0\n"
        "PRIMER_RIGHT_0=719,20\n"
        "PRIMER_RIGHT_0_SEQUENCE=GGGGGGGGGGGGGGGGGGGG\n"
        "PRIMER_RIGHT_0_TM=60.0\n"
        "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
        "PRIMER_PAIR_0_PRODUCT_SIZE=141\n"
        "=\n"
    )


class ContrastAnchorPrimerDesignTests(unittest.TestCase):
    def test_complete_contrastive_result_activates_primer3_anchor_mode(self) -> None:
        region = candidate()
        contrastive = SimpleNamespace(
            status="COMPLETE",
            candidates=(SimpleNamespace(region=region),),
        )
        conservation = SimpleNamespace(
            reference_id="ref",
            major_consensus="A" * 1200,
        )
        runner = LiteralRunner(anchored_output())

        with tempfile.TemporaryDirectory() as temporary:
            result = design_primers(
                conservation,
                PrimerDesignConfig(enabled=True, assays_per_region=1),
                Path(temporary),
                contrastive=contrastive,
                runner=runner,
            )

        self.assertEqual(result.candidate_source, "CONTRASTIVE_CONSERVATION")
        self.assertEqual(len(result.assays), 1)
        self.assertIn("SEQUENCE_INCLUDED_REGION=500,300\n", runner.inputs[0])
        self.assertIn("SEQUENCE_TARGET=600,100\n", runner.inputs[0])


if __name__ == "__main__":
    unittest.main()
