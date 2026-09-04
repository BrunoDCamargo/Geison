from __future__ import annotations

import unittest

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.primer3 import build_primer3_input, parse_primer3_output
from qpcr_pipeline.primer_design import CandidateRegion, PrimerDesignError


CONSENSUS = "A" * 1200


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


def single_assay_output(
    *,
    forward: tuple[int, int],
    probe: tuple[int, int],
    reverse: tuple[int, int],
    product_size: int,
) -> str:
    forward_position, forward_length = forward
    probe_position, probe_length = probe
    reverse_position, reverse_length = reverse
    return (
        "SEQUENCE_ID=contrast-region-001\n"
        "PRIMER_LEFT_NUM_RETURNED=1\n"
        "PRIMER_INTERNAL_NUM_RETURNED=1\n"
        "PRIMER_RIGHT_NUM_RETURNED=1\n"
        "PRIMER_PAIR_NUM_RETURNED=1\n"
        f"PRIMER_LEFT_0={forward_position},{forward_length}\n"
        f"PRIMER_LEFT_0_SEQUENCE={'A' * forward_length}\n"
        "PRIMER_LEFT_0_TM=60.0\n"
        "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
        f"PRIMER_INTERNAL_0={probe_position},{probe_length}\n"
        f"PRIMER_INTERNAL_0_SEQUENCE={'C' * probe_length}\n"
        "PRIMER_INTERNAL_0_TM=70.0\n"
        "PRIMER_INTERNAL_0_GC_PERCENT=50.0\n"
        f"PRIMER_RIGHT_0={reverse_position},{reverse_length}\n"
        f"PRIMER_RIGHT_0_SEQUENCE={'G' * reverse_length}\n"
        "PRIMER_RIGHT_0_TM=60.0\n"
        "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
        f"PRIMER_PAIR_0_PRODUCT_SIZE={product_size}\n"
        "=\n"
    )


class ContrastAnchorPrimer3Tests(unittest.TestCase):
    def test_contrastive_input_serializes_peak_as_sequence_target(self) -> None:
        text = build_primer3_input(
            CONSENSUS,
            (candidate(),),
            PrimerDesignConfig(),
            require_contrast_anchor=True,
        )

        self.assertIn("SEQUENCE_INCLUDED_REGION=500,300\n", text)
        self.assertIn("SEQUENCE_TARGET=600,100\n", text)

    def test_conservation_only_input_does_not_emit_sequence_target(self) -> None:
        text = build_primer3_input(
            CONSENSUS,
            (candidate(),),
            PrimerDesignConfig(),
            require_contrast_anchor=False,
        )

        self.assertNotIn("SEQUENCE_TARGET=", text)

    def test_rejects_amplicon_entirely_downstream_of_contrast_anchor(self) -> None:
        output = single_assay_output(
            forward=(713, 20),
            probe=(743, 25),
            reverse=(799, 20),
            product_size=87,
        )

        with self.assertRaisesRegex(PrimerDesignError, "contrast anchor"):
            parse_primer3_output(
                output,
                (candidate(),),
                CONSENSUS,
                require_contrast_anchor=True,
            )

    def test_accepts_amplicon_that_contains_complete_contrast_anchor(self) -> None:
        output = single_assay_output(
            forward=(579, 20),
            probe=(619, 25),
            reverse=(719, 20),
            product_size=141,
        )

        assays, _ = parse_primer3_output(
            output,
            (candidate(),),
            CONSENSUS,
            require_contrast_anchor=True,
        )

        self.assertEqual(len(assays), 1)
        self.assertEqual(assays[0].assay_id, "contrast-region-001-assay-001")


if __name__ == "__main__":
    unittest.main()
