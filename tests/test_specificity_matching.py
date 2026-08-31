import unittest

from qpcr_pipeline.config import SpecificityConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo
from qpcr_pipeline.specificity_matching import (
    SpecificityMatchingError,
    all_assay_hits,
    enumerate_compatible_hits,
)


class SpecificityMatchingTests(unittest.TestCase):
    @staticmethod
    def _oligo(sequence: str, start: int) -> DesignedOligo:
        return DesignedOligo(
            sequence=sequence,
            reference_start=start,
            reference_end=start + len(sequence) - 1,
            length=len(sequence),
            tm=60.0,
            gc_percent=50.0,
            penalty=None,
            metrics=(),
        )

    def _assay(
        self,
        assay_id: str = "a1",
        *,
        forward: str = "ACGT",
        probe: str = "TTAA",
        reverse: str = "AGTC",
    ) -> AssayCandidate:
        f = self._oligo(forward, 3)
        p = self._oligo(probe, 9)
        r = self._oligo(reverse, 15)
        return AssayCandidate(
            assay_id=assay_id,
            region_id="r1",
            primer3_index=0,
            forward_primer=f,
            probe=p,
            reverse_primer=r,
            product_size=r.reference_end - f.reference_start + 1,
            pair_penalty=None,
            metrics=(),
        )

    def test_finds_exact_forward_hit_and_source_coordinates(self):
        hits = enumerate_compatible_hits(
            "human",
            LocalSequenceRecord("off-1", "TTACGTACTT"),
            self._assay(),
            "FORWARD",
            SpecificityConfig(max_primer_mismatches=0),
        )
        exact = [h for h in hits if h.orientation == "FORWARD" and h.source_start == 3]
        self.assertEqual(len(exact), 1)
        self.assertEqual(exact[0].source_end, 6)
        self.assertEqual(exact[0].mismatch_positions, ())
        self.assertTrue(exact[0].exact_match)
        self.assertTrue(exact[0].compatible)

    def test_internal_mismatch_can_pass_but_three_prime_mismatch_is_rejected(self):
        assay = self._assay(forward="ACGT")
        config = SpecificityConfig(
            max_primer_mismatches=1,
            primer_3_prime_bases=2,
            reject_primer_3_prime_mismatch=True,
        )
        internal = enumerate_compatible_hits(
            "d", LocalSequenceRecord("internal", "ATGT"), assay, "FORWARD", config
        )
        three_prime = enumerate_compatible_hits(
            "d", LocalSequenceRecord("three-prime", "ACGA"), assay, "FORWARD", config
        )
        self.assertTrue(any(h.mismatch_positions == (2,) for h in internal))
        self.assertFalse(any(h.orientation == "FORWARD" for h in three_prime))

    def test_conservative_iupac_semantics_match_target_support_subset(self):
        assay = self._assay(forward="ARGT")
        config = SpecificityConfig(max_primer_mismatches=0)
        canonical = enumerate_compatible_hits(
            "d", LocalSequenceRecord("canonical", "AGGT"), assay, "FORWARD", config
        )
        ambiguous = enumerate_compatible_hits(
            "d", LocalSequenceRecord("ambiguous", "ANGT"), assay, "FORWARD", config
        )
        self.assertTrue(any(h.orientation == "FORWARD" for h in canonical))
        self.assertFalse(any(h.orientation == "FORWARD" for h in ambiguous))

    def test_reverse_complement_orientation_maps_coordinates_back_to_source(self):
        assay = self._assay(forward="ACGT")
        # Source reverse complement is ACGT; the physical site spans the full record.
        hits = enumerate_compatible_hits(
            "d", LocalSequenceRecord("rev", "ACGT"), assay, "FORWARD",
            SpecificityConfig(max_primer_mismatches=0),
        )
        reverse_hits = [h for h in hits if h.orientation == "REVERSE_COMPLEMENT"]
        self.assertTrue(reverse_hits)
        self.assertEqual((reverse_hits[0].source_start, reverse_hits[0].source_end), (1, 4))

    def test_assay_and_record_order_is_preserved_not_lexically_resorted(self):
        assays = (self._assay("z-assay"), self._assay("a-assay"))
        records = (
            LocalSequenceRecord("z-record", "ACGT"),
            LocalSequenceRecord("a-record", "ACGT"),
        )
        hits = all_assay_hits(
            "z-dataset", records, assays,
            SpecificityConfig(max_primer_mismatches=0, max_probe_mismatches=0),
        )
        first_seen = []
        for hit in hits:
            key = (hit.assay_id, hit.sequence_id)
            if key not in first_seen:
                first_seen.append(key)
        self.assertEqual(
            first_seen,
            [
                ("z-assay", "z-record"),
                ("z-assay", "a-record"),
                ("a-assay", "z-record"),
                ("a-assay", "a-record"),
            ],
        )

    def test_invalid_iupac_error_is_contextual_without_dumping_full_sequence(self):
        bad = "ACGTXACGT"
        with self.assertRaisesRegex(
            SpecificityMatchingError,
            "dataset 'human'.*assay 'a1'.*sequence 'bad'.*role 'FORWARD'",
        ) as caught:
            enumerate_compatible_hits(
                "human", LocalSequenceRecord("bad", bad), self._assay(), "FORWARD",
                SpecificityConfig(),
            )
        self.assertNotIn(bad, str(caught.exception))


if __name__ == "__main__":
    unittest.main()
