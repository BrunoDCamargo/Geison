import shutil
import unittest

from qpcr_pipeline.config import OligoConstraints, PrimerDesignConfig
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
INTEGRATION_CONFIG = PrimerDesignConfig(
    enabled=True,
    assays_per_region=5,
    product_size_min=70,
    product_size_max=220,
    primer=OligoConstraints(
        min_size=16,
        opt_size=20,
        max_size=28,
        min_tm=50.0,
        opt_tm=60.0,
        max_tm=72.0,
        min_gc_percent=20.0,
        max_gc_percent=80.0,
    ),
    probe=OligoConstraints(
        min_size=16,
        opt_size=24,
        max_size=32,
        min_tm=55.0,
        opt_tm=65.0,
        max_tm=80.0,
        min_gc_percent=20.0,
        max_gc_percent=80.0,
    ),
)


@unittest.skipUnless(
    shutil.which("primer3_core"),
    "primer3_core is not installed or is absent from PATH",
)
class Primer3CoreIntegrationTests(unittest.TestCase):
    def test_real_binary_returns_complete_assays_with_consistent_coordinates(self):
        candidates = (CANDIDATE,)
        input_text = build_primer3_input(
            SYNTHETIC_CONSENSUS,
            candidates,
            INTEGRATION_CONFIG,
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
