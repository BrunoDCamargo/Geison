import os
from pathlib import Path
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest.mock import patch

from qpcr_pipeline.config import OligoConstraints, PrimerDesignConfig
from qpcr_pipeline.primer3 import (
    SubprocessPrimer3Runner,
    build_primer3_input,
    parse_primer3_output,
)
from qpcr_pipeline.primer_design import CandidateRegion, PrimerDesignError


CONSENSUS_300 = (
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
)


def candidate(
    region_id: str = "region-001",
    reference_start: int = 1,
    reference_end: int = 300,
) -> CandidateRegion:
    return CandidateRegion(
        region_id=region_id,
        rank=1,
        reference_start=reference_start,
        reference_end=reference_end,
        peak_start=101,
        peak_end=200,
        position_count=reference_end - reference_start + 1,
        usable_length=reference_end - reference_start + 1,
        usable_fraction=1.0,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


class BuildPrimer3InputTests(unittest.TestCase):
    def test_serializes_complete_record_in_fixed_order(self) -> None:
        text = build_primer3_input(
            CONSENSUS_300,
            (candidate(),),
            PrimerDesignConfig(),
        )

        self.assertEqual(
            text,
            "SEQUENCE_ID=region-001\n"
            "SEQUENCE_TEMPLATE=ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n"
            "SEQUENCE_INCLUDED_REGION=0,300\n"
            "PRIMER_TASK=generic\n"
            "PRIMER_FIRST_BASE_INDEX=0\n"
            "PRIMER_PICK_LEFT_PRIMER=1\n"
            "PRIMER_PICK_INTERNAL_OLIGO=1\n"
            "PRIMER_PICK_RIGHT_PRIMER=1\n"
            "PRIMER_NUM_RETURN=5\n"
            "PRIMER_PRODUCT_SIZE_RANGE=70-200\n"
            "PRIMER_EXPLAIN_FLAG=1\n"
            "PRIMER_MIN_SIZE=18\n"
            "PRIMER_OPT_SIZE=20\n"
            "PRIMER_MAX_SIZE=25\n"
            "PRIMER_MIN_TM=58.0\n"
            "PRIMER_OPT_TM=60.0\n"
            "PRIMER_MAX_TM=62.0\n"
            "PRIMER_MIN_GC=40.0\n"
            "PRIMER_MAX_GC=60.0\n"
            "PRIMER_INTERNAL_MIN_SIZE=18\n"
            "PRIMER_INTERNAL_OPT_SIZE=25\n"
            "PRIMER_INTERNAL_MAX_SIZE=30\n"
            "PRIMER_INTERNAL_MIN_TM=68.0\n"
            "PRIMER_INTERNAL_OPT_TM=70.0\n"
            "PRIMER_INTERNAL_MAX_TM=72.0\n"
            "PRIMER_INTERNAL_MIN_GC=30.0\n"
            "PRIMER_INTERNAL_MAX_GC=80.0\n"
            "=\n",
        )

    def test_rejects_candidate_ids_with_boulder_delimiters(self) -> None:
        for unsafe_id in ("region=001", "region\n001", "region\r001"):
            with self.subTest(unsafe_id=repr(unsafe_id)):
                with self.assertRaisesRegex(
                    PrimerDesignError, "candidate region ID contains unsafe"
                ):
                    build_primer3_input(
                        CONSENSUS_300,
                        (candidate(region_id=unsafe_id),),
                        PrimerDesignConfig(),
                    )

    def test_rejects_empty_or_noncanonical_consensus(self) -> None:
        for unsafe_consensus in ("", "ACGTN", "ACGT=ACGT", "ACGT\nACGT", "acgt"):
            with self.subTest(unsafe_consensus=repr(unsafe_consensus)):
                with self.assertRaisesRegex(
                    PrimerDesignError,
                    "consensus must contain only uppercase canonical DNA bases",
                ):
                    build_primer3_input(
                        unsafe_consensus,
                        (candidate(),),
                        PrimerDesignConfig(),
                    )

    def test_serializes_all_configured_constraints_and_included_region(self) -> None:
        config = PrimerDesignConfig(
            assays_per_region=3,
            product_size_min=90,
            product_size_max=180,
            primer=OligoConstraints(
                19, 21, 24, 57.5, 59.5, 61.5, 35.5, 65.5
            ),
            probe=OligoConstraints(
                20, 24, 28, 67.5, 69.5, 71.5, 32.5, 77.5
            ),
        )

        text = build_primer3_input(
            CONSENSUS_300,
            (candidate(reference_start=51, reference_end=250),),
            config,
        )

        self.assertEqual(
            text.splitlines()[2:],
            [
                "SEQUENCE_INCLUDED_REGION=50,200",
                "PRIMER_TASK=generic",
                "PRIMER_FIRST_BASE_INDEX=0",
                "PRIMER_PICK_LEFT_PRIMER=1",
                "PRIMER_PICK_INTERNAL_OLIGO=1",
                "PRIMER_PICK_RIGHT_PRIMER=1",
                "PRIMER_NUM_RETURN=3",
                "PRIMER_PRODUCT_SIZE_RANGE=90-180",
                "PRIMER_EXPLAIN_FLAG=1",
                "PRIMER_MIN_SIZE=19",
                "PRIMER_OPT_SIZE=21",
                "PRIMER_MAX_SIZE=24",
                "PRIMER_MIN_TM=57.5",
                "PRIMER_OPT_TM=59.5",
                "PRIMER_MAX_TM=61.5",
                "PRIMER_MIN_GC=35.5",
                "PRIMER_MAX_GC=65.5",
                "PRIMER_INTERNAL_MIN_SIZE=20",
                "PRIMER_INTERNAL_OPT_SIZE=24",
                "PRIMER_INTERNAL_MAX_SIZE=28",
                "PRIMER_INTERNAL_MIN_TM=67.5",
                "PRIMER_INTERNAL_OPT_TM=69.5",
                "PRIMER_INTERNAL_MAX_TM=71.5",
                "PRIMER_INTERNAL_MIN_GC=32.5",
                "PRIMER_INTERNAL_MAX_GC=77.5",
                "=",
            ],
        )


class ParsePrimer3OutputTests(unittest.TestCase):
    def test_parses_complete_assay_with_exact_coordinates_and_metrics(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.25\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_LEFT_0_PENALTY=0.125\n"
            "PRIMER_LEFT_0_SELF_ANY_TH=0.7\n"
            "PRIMER_LEFT_0_END_STABILITY=3.1\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.5\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_INTERNAL_0_PENALTY=0.25\n"
            "PRIMER_INTERNAL_0_HAIRPIN_TH=0.2\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=59.75\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_RIGHT_0_PENALTY=0.375\n"
            "PRIMER_RIGHT_0_SELF_END_TH=0.4\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "PRIMER_PAIR_0_PENALTY=0.75\n"
            "PRIMER_PAIR_0_TEMPLATE_MISPRIMING_TH=1.2\n"
            "PRIMER_PAIR_0_COMPL_ANY_TH=0.0\n"
            "=\n"
        )

        assays, details = parse_primer3_output(
            output,
            (candidate(),),
            CONSENSUS_300,
        )

        self.assertEqual(details, {"region-001": {}})
        self.assertEqual(len(assays), 1)
        assay = assays[0]
        self.assertEqual(assay.assay_id, "region-001-assay-001")
        self.assertEqual(assay.region_id, "region-001")
        self.assertEqual(assay.primer3_index, 0)
        self.assertEqual(assay.product_size, 80)
        self.assertEqual(assay.pair_penalty, 0.75)
        self.assertEqual(
            assay.metrics,
            (
                ("PRIMER_PAIR_0_COMPL_ANY_TH", "0.0"),
                ("PRIMER_PAIR_0_TEMPLATE_MISPRIMING_TH", "1.2"),
            ),
        )

        self.assertEqual(assay.forward_primer.sequence, "ACGTACGTACGTACGTACGT")
        self.assertEqual(
            (assay.forward_primer.reference_start, assay.forward_primer.reference_end),
            (11, 30),
        )
        self.assertEqual(assay.forward_primer.length, 20)
        self.assertEqual(assay.forward_primer.tm, 60.25)
        self.assertEqual(assay.forward_primer.gc_percent, 50.0)
        self.assertEqual(assay.forward_primer.penalty, 0.125)
        self.assertEqual(
            assay.forward_primer.metrics,
            (
                ("PRIMER_LEFT_0_END_STABILITY", "3.1"),
                ("PRIMER_LEFT_0_SELF_ANY_TH", "0.7"),
            ),
        )

        self.assertEqual(assay.probe.sequence, "ACGTACGTACGTACGTACGTACGTA")
        self.assertEqual(
            (assay.probe.reference_start, assay.probe.reference_end),
            (36, 60),
        )
        self.assertEqual(assay.probe.length, 25)
        self.assertEqual(assay.probe.tm, 70.5)
        self.assertEqual(assay.probe.gc_percent, 48.0)
        self.assertEqual(assay.probe.penalty, 0.25)
        self.assertEqual(
            assay.probe.metrics,
            (("PRIMER_INTERNAL_0_HAIRPIN_TH", "0.2"),),
        )

        self.assertEqual(assay.reverse_primer.sequence, "TGCATGCATGCATGCATGCA")
        self.assertEqual(
            (assay.reverse_primer.reference_start, assay.reverse_primer.reference_end),
            (71, 90),
        )
        self.assertEqual(assay.reverse_primer.length, 20)
        self.assertEqual(assay.reverse_primer.tm, 59.75)
        self.assertEqual(assay.reverse_primer.gc_percent, 50.0)
        self.assertEqual(assay.reverse_primer.penalty, 0.375)
        self.assertEqual(
            assay.reverse_primer.metrics,
            (("PRIMER_RIGHT_0_SELF_END_TH", "0.4"),),
        )

    def test_matches_out_of_order_records_by_sequence_id(self) -> None:
        output = (
            "SEQUENCE_ID=region-002\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        assays, details = parse_primer3_output(
            output,
            (candidate(), candidate(region_id="region-002")),
            CONSENSUS_300,
        )

        self.assertEqual(assays, ())
        self.assertEqual(details, {"region-001": {}, "region-002": {}})

    def test_parses_multiple_assays_in_primer3_index_order(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=2\n"
            "PRIMER_INTERNAL_NUM_RETURNED=2\n"
            "PRIMER_RIGHT_NUM_RETURNED=2\n"
            "PRIMER_PAIR_NUM_RETURNED=2\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "PRIMER_LEFT_1=100,20\n"
            "PRIMER_LEFT_1_SEQUENCE=AAAACCCCGGGGTTTTAAAA\n"
            "PRIMER_LEFT_1_TM=59.0\n"
            "PRIMER_LEFT_1_GC_PERCENT=40.0\n"
            "PRIMER_INTERNAL_1=130,20\n"
            "PRIMER_INTERNAL_1_SEQUENCE=CCCCAAAATTTTGGGGCCCC\n"
            "PRIMER_INTERNAL_1_TM=69.0\n"
            "PRIMER_INTERNAL_1_GC_PERCENT=60.0\n"
            "PRIMER_RIGHT_1=169,20\n"
            "PRIMER_RIGHT_1_SEQUENCE=TTTTGGGGCCCCAAAATTTT\n"
            "PRIMER_RIGHT_1_TM=61.0\n"
            "PRIMER_RIGHT_1_GC_PERCENT=40.0\n"
            "PRIMER_PAIR_1_PRODUCT_SIZE=70\n"
            "=\n"
        )

        assays, _ = parse_primer3_output(output, (candidate(),), CONSENSUS_300)

        self.assertEqual(
            [assay.assay_id for assay in assays],
            ["region-001-assay-001", "region-001-assay-002"],
        )
        self.assertEqual([assay.primer3_index for assay in assays], [0, 1])
        self.assertEqual([assay.product_size for assay in assays], [80, 70])

    def test_accepts_zero_pairs_and_preserves_explanations_and_warning(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_EXPLAIN=considered 12, low tm 12, ok 0\n"
            "PRIMER_INTERNAL_EXPLAIN=considered 7, high tm 7, ok 0\n"
            "PRIMER_RIGHT_EXPLAIN=considered 12, high hairpin stability 12, ok 0\n"
            "PRIMER_PAIR_EXPLAIN=considered 0, ok 0\n"
            "PRIMER_WARNING=design space was constrained\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        assays, details = parse_primer3_output(output, (candidate(),), CONSENSUS_300)

        self.assertEqual(assays, ())
        self.assertEqual(
            details,
            {
                "region-001": {
                    "PRIMER_LEFT_EXPLAIN": "considered 12, low tm 12, ok 0",
                    "PRIMER_INTERNAL_EXPLAIN": "considered 7, high tm 7, ok 0",
                    "PRIMER_RIGHT_EXPLAIN": (
                        "considered 12, high hairpin stability 12, ok 0"
                    ),
                    "PRIMER_PAIR_EXPLAIN": "considered 0, ok 0",
                    "PRIMER_WARNING": "design space was constrained",
                }
            },
        )

    def test_rejects_global_primer3_error_tag(self) -> None:
        output = "PRIMER_ERROR=Unrecognized tag PRIMER_BAD_OPTION\n=\n"

        with self.assertRaisesRegex(
            PrimerDesignError,
            "Primer3 reported an error: Unrecognized tag PRIMER_BAD_OPTION",
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_duplicate_output_sequence_ids(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "duplicate SEQUENCE_ID 'region-001'"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_missing_requested_sequence_id(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "missing requested SEQUENCE_ID 'region-002'"
        ):
            parse_primer3_output(
                output,
                (candidate(), candidate(region_id="region-002")),
                CONSENSUS_300,
            )

    def test_rejects_unknown_output_sequence_id(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
            "SEQUENCE_ID=region-999\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "unknown SEQUENCE_ID 'region-999'"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_partial_assay_with_missing_member_field(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError,
            "missing required tag 'PRIMER_INTERNAL_0_SEQUENCE'",
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_malformed_numeric_values(self) -> None:
        valid_output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_LEFT_0_PENALTY=0.1\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "PRIMER_PAIR_0_PENALTY=0.2\n"
            "=\n"
        )
        mutations = (
            ("pair count", "PRIMER_PAIR_NUM_RETURNED=1", "PRIMER_PAIR_NUM_RETURNED=one"),
            ("coordinate", "PRIMER_LEFT_0=10,20", "PRIMER_LEFT_0=ten,20"),
            ("Tm", "PRIMER_LEFT_0_TM=60.0", "PRIMER_LEFT_0_TM=hot"),
            ("GC", "PRIMER_LEFT_0_GC_PERCENT=50.0", "PRIMER_LEFT_0_GC_PERCENT=half"),
            ("oligo penalty", "PRIMER_LEFT_0_PENALTY=0.1", "PRIMER_LEFT_0_PENALTY=low"),
            ("product size", "PRIMER_PAIR_0_PRODUCT_SIZE=80", "PRIMER_PAIR_0_PRODUCT_SIZE=eighty"),
            ("pair penalty", "PRIMER_PAIR_0_PENALTY=0.2", "PRIMER_PAIR_0_PENALTY=low"),
        )

        for label, old, new in mutations:
            with self.subTest(label=label):
                with self.assertRaisesRegex(
                    PrimerDesignError, "invalid numeric value"
                ):
                    parse_primer3_output(
                        valid_output.replace(old, new),
                        (candidate(),),
                        CONSENSUS_300,
                    )

    def test_rejects_nonpositive_oligo_length(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,0\n"
            "PRIMER_LEFT_0_SEQUENCE=A\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_LEFT_0.*positive length"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_oligo_sequence_length_mismatch(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACG\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_LEFT_0_SEQUENCE.*does not match.*20"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_oligo_coordinates_outside_candidate_region(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=80,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=129,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=120\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_LEFT_0.*outside candidate 'region-001'"
        ):
            parse_primer3_output(
                output,
                (candidate(reference_start=51, reference_end=250),),
                CONSENSUS_300,
            )

    def test_rejects_product_size_disagreeing_with_primer_coordinates(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=79\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "product size 79.*coordinate-derived size 80"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_output_record_missing_sequence_id(self) -> None:
        output = (
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=0\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "missing required tag 'SEQUENCE_ID'"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_malformed_boulder_record_structure(self) -> None:
        malformed_outputs = (
            ("empty", ""),
            (
                "unterminated",
                "SEQUENCE_ID=region-001\nPRIMER_PAIR_NUM_RETURNED=0\n",
            ),
            (
                "line without delimiter",
                "SEQUENCE_ID=region-001\nBROKEN LINE\n=\n",
            ),
            (
                "duplicate tag",
                "SEQUENCE_ID=region-001\n"
                "SEQUENCE_ID=region-001\n"
                "PRIMER_PAIR_NUM_RETURNED=0\n"
                "=\n",
            ),
        )

        for label, output in malformed_outputs:
            with self.subTest(label=label):
                with self.assertRaisesRegex(
                    PrimerDesignError, "malformed Boulder-IO"
                ):
                    parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_pair_count_exceeding_member_count(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError,
            "PRIMER_PAIR_NUM_RETURNED.*exceeds PRIMER_INTERNAL_NUM_RETURNED",
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_negative_returned_count(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=0\n"
            "PRIMER_INTERNAL_NUM_RETURNED=0\n"
            "PRIMER_RIGHT_NUM_RETURNED=0\n"
            "PRIMER_PAIR_NUM_RETURNED=-1\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_PAIR_NUM_RETURNED.*nonnegative"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_nonfinite_numeric_value(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=nan\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_LEFT_0_TM.*finite numeric value"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_coordinate_value_without_position_length_pair(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=89,20\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATGCATGCATGCA\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=80\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "PRIMER_LEFT_0.*position,length pair"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)

    def test_rejects_nonpositive_coordinate_derived_product_size(self) -> None:
        output = (
            "SEQUENCE_ID=region-001\n"
            "PRIMER_LEFT_NUM_RETURNED=1\n"
            "PRIMER_INTERNAL_NUM_RETURNED=1\n"
            "PRIMER_RIGHT_NUM_RETURNED=1\n"
            "PRIMER_PAIR_NUM_RETURNED=1\n"
            "PRIMER_LEFT_0=10,20\n"
            "PRIMER_LEFT_0_SEQUENCE=ACGTACGTACGTACGTACGT\n"
            "PRIMER_LEFT_0_TM=60.0\n"
            "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
            "PRIMER_INTERNAL_0=35,25\n"
            "PRIMER_INTERNAL_0_SEQUENCE=ACGTACGTACGTACGTACGTACGTA\n"
            "PRIMER_INTERNAL_0_TM=70.0\n"
            "PRIMER_INTERNAL_0_GC_PERCENT=48.0\n"
            "PRIMER_RIGHT_0=9,10\n"
            "PRIMER_RIGHT_0_SEQUENCE=TGCATGCATG\n"
            "PRIMER_RIGHT_0_TM=60.0\n"
            "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
            "PRIMER_PAIR_0_PRODUCT_SIZE=0\n"
            "=\n"
        )

        with self.assertRaisesRegex(
            PrimerDesignError, "product size must be positive"
        ):
            parse_primer3_output(output, (candidate(),), CONSENSUS_300)


class SubprocessPrimer3RunnerTests(unittest.TestCase):
    def test_overlapping_runs_do_not_mutate_process_thread_exception_hook(self) -> None:
        first_entered = threading.Event()
        second_entered = threading.Event()
        release_first = threading.Event()
        release_second = threading.Event()
        call_lock = threading.Lock()
        call_count = 0
        outputs: list[str] = []
        errors: list[BaseException] = []
        original_hook = threading.excepthook
        self.addCleanup(setattr, threading, "excepthook", original_hook)

        def coordinated_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess:
            nonlocal call_count
            with call_lock:
                call_count += 1
                call_number = call_count
            if call_number == 1:
                first_entered.set()
                release_first.wait(timeout=5.0)
            else:
                second_entered.set()
                release_second.wait(timeout=5.0)
            return self._successful_completed_process(kwargs["input"])

        def invoke() -> None:
            try:
                outputs.append(
                    SubprocessPrimer3Runner(sys.executable).run(
                        "SEQUENCE_ID=region-001\n=\n"
                    )
                )
            except BaseException as error:
                errors.append(error)

        with patch("qpcr_pipeline.primer3.subprocess.run", coordinated_run):
            first = threading.Thread(target=invoke)
            second = threading.Thread(target=invoke)
            first.start()
            self.assertTrue(first_entered.wait(timeout=5.0))
            second.start()
            self.assertTrue(second_entered.wait(timeout=5.0))
            release_first.set()
            first.join(timeout=5.0)
            self.assertFalse(first.is_alive())
            release_second.set()
            second.join(timeout=5.0)
            self.assertFalse(second.is_alive())

        self.assertEqual(errors, [])
        self.assertEqual(outputs, ["=\n", "=\n"])
        self.assertIs(threading.excepthook, original_hook)

    def test_unrelated_thread_unicode_error_is_not_attributed_to_runner(self) -> None:
        runner_entered = threading.Event()
        release_runner = threading.Event()
        runner_outputs: list[str] = []
        runner_errors: list[BaseException] = []
        unrelated_errors: list[BaseException] = []
        original_hook = threading.excepthook

        def record_unrelated_error(args: threading.ExceptHookArgs) -> None:
            unrelated_errors.append(args.exc_value)

        threading.excepthook = record_unrelated_error
        self.addCleanup(setattr, threading, "excepthook", original_hook)

        def blocking_run(*args: object, **kwargs: object) -> subprocess.CompletedProcess:
            runner_entered.set()
            release_runner.wait(timeout=5.0)
            return self._successful_completed_process(kwargs["input"])

        def invoke_runner() -> None:
            try:
                runner_outputs.append(
                    SubprocessPrimer3Runner(sys.executable).run(
                        "SEQUENCE_ID=region-001\n=\n"
                    )
                )
            except BaseException as error:
                runner_errors.append(error)

        def raise_unrelated_error() -> None:
            raise UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid byte")

        with patch("qpcr_pipeline.primer3.subprocess.run", blocking_run):
            runner = threading.Thread(target=invoke_runner)
            runner.start()
            self.assertTrue(runner_entered.wait(timeout=5.0))
            unrelated = threading.Thread(target=raise_unrelated_error)
            unrelated.start()
            unrelated.join(timeout=5.0)
            self.assertFalse(unrelated.is_alive())
            release_runner.set()
            runner.join(timeout=5.0)
            self.assertFalse(runner.is_alive())

        self.assertEqual(runner_errors, [])
        self.assertEqual(runner_outputs, ["=\n"])
        self.assertEqual(len(unrelated_errors), 1)
        self.assertIsInstance(unrelated_errors[0], UnicodeDecodeError)
        self.assertIs(threading.excepthook, record_unrelated_error)

    def test_uses_fixed_arguments_and_exchanges_stdin_stdout(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    'if not "%~1"=="--strict_tags" exit /b 31\n'
                    'if not "%~2=%~3"=="--io_version=4" exit /b 32\n'
                    'if not "%~4"=="" exit /b 33\n'
                    "echo(%cmdcmdline% | findstr /L /C:\"--io_version=4\" >nul "
                    "|| exit /b 34\n"
                    "powershell.exe -NoProfile -NonInteractive -Command "
                    '"$text = [Console]::In.ReadToEnd(); '
                    '[Console]::Out.Write($text)"\n'
                ),
                posix=(
                    "#!/bin/sh\n"
                    '[ "$1" = "--strict_tags" ] || exit 31\n'
                    '[ "$2" = "--io_version=4" ] || exit 32\n'
                    '[ "$#" -eq 2 ] || exit 33\n'
                    "cat\n"
                ),
            )
            input_text = "SEQUENCE_ID=region-001\n=\n"

            output = SubprocessPrimer3Runner(str(executable)).run(input_text)

        self.assertEqual(output, input_text)

    def test_missing_executable_raises_domain_error_without_input(self) -> None:
        input_text = f"SEQUENCE_TEMPLATE={CONSENSUS_300}\n=\n"

        with self.assertRaisesRegex(
            PrimerDesignError,
            "Primer3 executable 'definitely-missing-primer3-core-for-tests' was not found",
        ) as raised:
            SubprocessPrimer3Runner(
                "definitely-missing-primer3-core-for-tests"
            ).run(input_text)

        self.assertNotIn(CONSENSUS_300, str(raised.exception))

    def test_nonzero_exit_bounds_stderr_excerpt_to_2000_characters(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    "powershell.exe -NoProfile -NonInteractive -Command "
                    '"[Console]::Error.Write((\'X\' * 2500))"\n'
                    "exit /b 7\n"
                ),
                posix=(
                    "#!/bin/sh\n"
                    "i=0\n"
                    "while [ \"$i\" -lt 2500 ]; do printf X >&2; i=$((i + 1)); done\n"
                    "exit 7\n"
                ),
            )
            input_text = f"SEQUENCE_TEMPLATE={CONSENSUS_300}\n=\n"

            with self.assertRaisesRegex(
                PrimerDesignError, "Primer3 exited with status 7"
            ) as raised:
                SubprocessPrimer3Runner(str(executable)).run(input_text)

        message = str(raised.exception)
        self.assertIn("X" * 2000, message)
        self.assertNotIn("X" * 2001, message)
        self.assertNotIn(CONSENSUS_300, message)

    def test_empty_stdout_raises_domain_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    "powershell.exe -NoProfile -NonInteractive -Command "
                    '"[void][Console]::In.ReadToEnd()"\n'
                ),
                posix="#!/bin/sh\ncat >/dev/null\n",
            )

            with self.assertRaisesRegex(
                PrimerDesignError, "Primer3 produced empty stdout"
            ):
                SubprocessPrimer3Runner(str(executable)).run(
                    "SEQUENCE_ID=region-001\n=\n"
                )

    def test_invalid_utf8_stdout_raises_domain_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    "powershell.exe -NoProfile -NonInteractive -Command "
                    '"[void][Console]::In.ReadToEnd(); '
                    "$stdout = [Console]::OpenStandardOutput(); "
                    '$stdout.WriteByte(255)"\n'
                ),
                posix="#!/bin/sh\ncat >/dev/null\nprintf '\\377'\n",
            )

            with self.assertRaisesRegex(
                PrimerDesignError, "Primer3 output was not valid UTF-8"
            ):
                SubprocessPrimer3Runner(str(executable)).run(
                    "SEQUENCE_ID=region-001\n=\n"
                )

    def test_unexpected_execution_error_is_normalized_without_input(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    "echo =\n"
                ),
                posix="#!/bin/sh\nprintf '=\\n'\n",
            )
            input_text = f"SEQUENCE_TEMPLATE={CONSENSUS_300}\ud800\n=\n"

            with self.assertRaisesRegex(
                PrimerDesignError,
                r"Primer3 execution failed \(UnicodeEncodeError\)",
            ) as raised:
                SubprocessPrimer3Runner(str(executable)).run(input_text)

        self.assertNotIn(CONSENSUS_300, str(raised.exception))

    def test_nonzero_exit_redacts_consensus_from_stderr(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            executable = self._write_executable(
                Path(temporary_directory),
                windows=(
                    "@echo off\n"
                    "powershell.exe -NoProfile -NonInteractive -Command "
                    '"$text = [Console]::In.ReadToEnd(); '
                    '[Console]::Error.Write($text)"\n'
                    "exit /b 9\n"
                ),
                posix="#!/bin/sh\ncat >&2\nexit 9\n",
            )
            input_text = f"SEQUENCE_TEMPLATE={CONSENSUS_300}\n=\n"

            with self.assertRaisesRegex(
                PrimerDesignError, "Primer3 exited with status 9"
            ) as raised:
                SubprocessPrimer3Runner(str(executable)).run(input_text)

        self.assertNotIn(CONSENSUS_300, str(raised.exception))
        self.assertIn("<consensus omitted>", str(raised.exception))

    @staticmethod
    def _successful_completed_process(
        input_value: object,
    ) -> subprocess.CompletedProcess:
        if isinstance(input_value, bytes):
            return subprocess.CompletedProcess((), 0, stdout=b"=\n", stderr=b"")
        return subprocess.CompletedProcess((), 0, stdout="=\n", stderr="")

    @staticmethod
    def _write_executable(
        directory: Path,
        *,
        windows: str,
        posix: str,
    ) -> Path:
        executable = directory / (
            "fake_primer3.cmd" if os.name == "nt" else "fake_primer3"
        )
        executable.write_text(
            windows if os.name == "nt" else posix,
            encoding="utf-8",
            newline="\n",
        )
        if os.name != "nt":
            executable.chmod(0o755)
        return executable


if __name__ == "__main__":
    unittest.main()
