import shutil
import unittest

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.primer3 import (
    SubprocessPrimer3Runner,
    build_primer3_input,
    parse_primer3_output,
)
from qpcr_pipeline.primer_design import CandidateRegion


SYNTHETIC_CONSENSUS = (
    "GAACTAACCGGAACCGCTAGGGTTCTCAGTTTGGGGCGCTGTCCTCCCCGAGGTCTACATGAGG"
    "TGGGTGCTGCTGCGTTGTTAGACTATGTGCAGCGCATGGGCCACTGCGCAAGCGTTCTGAACAAA"
    "CCTCTAACCGATTTCCAGCCCCGGCTATGCACACCCGAGATCAATAGTAGAGTGTATATTAATTC"
    "TTGGATTAACACGTCCATTCGCCCCGCCGATCAGAAACTGAGGCGTACCGTCTCTCCAAGCCCTG"
    "GATGAATTAATGTTGCAATCCGGTGTTACCTAGTAGCGACTGGCTGCACCCCGGCTCTACCGGTC"
    "TTAGTACACCTACATCGGTTGGAGACCGGAATGTAT"
)
CANDIDATE = CandidateRegion(
    region_id="region-001",
    rank=1,
    reference_start=31,
    reference_end=330,
    peak_start=131,
    peak_end=230,
    position_count=300,
    usable_length=300,
    usable_fraction=1.0,
    mean_conservation=1.0,
    minimum_conservation=1.0,
    mean_coverage=1.0,
    mean_gap_frequency=0.0,
    mean_entropy_bits=0.0,
)


@unittest.skipUnless(
    shutil.which("primer3_core"),
    "primer3_core is not installed or is absent from PATH",
)
class Primer3CoreIntegrationTests(unittest.TestCase):
    def test_real_binary_returns_complete_assays_with_consistent_coordinates(self):
        config = PrimerDesignConfig(enabled=True, assays_per_region=5)
        candidates = (CANDIDATE,)
        input_text = build_primer3_input(
            SYNTHETIC_CONSENSUS,
            candidates,
            config,
        )

        output_text = SubprocessPrimer3Runner().run(input_text)
        assays, details = parse_primer3_output(
            output_text,
            candidates,
            SYNTHETIC_CONSENSUS,
        )

        self.assertEqual(set(details), {CANDIDATE.region_id})
        self.assertGreaterEqual(len(assays), 1)
        for assay in assays:
            self.assertEqual(assay.region_id, CANDIDATE.region_id)
            self.assertTrue(assay.forward_primer.sequence)
            self.assertTrue(assay.probe.sequence)
            self.assertTrue(assay.reverse_primer.sequence)
            for oligo in (
                assay.forward_primer,
                assay.probe,
                assay.reverse_primer,
            ):
                self.assertGreaterEqual(
                    oligo.reference_start,
                    CANDIDATE.reference_start,
                )
                self.assertLessEqual(oligo.reference_end, CANDIDATE.reference_end)
            self.assertEqual(
                assay.product_size,
                assay.reverse_primer.reference_end
                - assay.forward_primer.reference_start
                + 1,
            )


if __name__ == "__main__":
    unittest.main()
