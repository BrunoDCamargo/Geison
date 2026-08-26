from dataclasses import fields
import unittest

from qpcr_pipeline.config import InclusivityConfig
from qpcr_pipeline.inclusivity import (
    AssayInclusivity,
    DegeneracyProposal,
    InclusivityError,
    InclusivityResult,
    OligoMatch,
    OligoVariation,
    ProposedOligoCompatibility,
    _Hit,
    _best_fallback,
    _enumerate_hits,
    _evaluate_original,
    _oligo_for_role,
    _select_binding,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import EvaluationSet
from qpcr_pipeline.primer_design import AssayCandidate, DesignedOligo


class InclusivitySearchTests(unittest.TestCase):
    def _oligo(self, sequence: str, start: int, end: int | None = None) -> DesignedOligo:
        end = end if end is not None else start + len(sequence) - 1
        return DesignedOligo(
            sequence=sequence,
            reference_start=start,
            reference_end=end,
            length=len(sequence),
            tm=60.0,
            gc_percent=50.0,
            penalty=None,
            metrics=(),
        )

    def _assay(self) -> AssayCandidate:
        return AssayCandidate(
            assay_id="a1",
            region_id="r1",
            primer3_index=0,
            forward_primer=self._oligo("ACGT", 3, 6),
            probe=self._oligo("TTAA", 9, 12),
            reverse_primer=self._oligo("AGTC", 15, 18),
            product_size=16,
            pair_penalty=None,
            metrics=(),
        )

    def _hits(
        self,
        sequence: str,
        *,
        oligo: DesignedOligo,
        role: str = "FORWARD",
        orientation: str = "FORWARD",
        config: InclusivityConfig | None = None,
    ):
        return _enumerate_hits(
            assay_id="a1",
            sequence_id="s1",
            oriented_sequence=sequence,
            orientation=orientation,
            role=role,
            oligo=oligo,
            config=config or InclusivityConfig(search_flank=0, max_hits_per_oligo=20),
        )

    def test_public_models_keep_only_public_result_fields(self):
        expected_fields = {
            OligoMatch: (
                "assay_id", "sequence_id", "role", "orientation", "hit_rank",
                "source_start", "source_end", "expected_start", "expected_end",
                "displacement", "mismatch_positions", "mismatch_count", "exact_match",
                "three_prime_mismatch", "probe_mismatch", "compatible", "selected",
            ),
            ProposedOligoCompatibility: (
                "role", "effective_sequence", "exact_match", "mismatch_positions",
                "mismatch_count", "three_prime_mismatch", "probe_mismatch", "compatible",
            ),
            AssayInclusivity: (
                "assay_id", "sequence_id", "orientation", "geometry_found",
                "source_amplicon_start", "source_amplicon_end", "amplicon_size",
                "forward_match", "probe_match", "reverse_match", "original_compatible",
                "proposed_forward", "proposed_probe", "proposed_reverse", "proposed_compatible",
            ),
            OligoVariation: (
                "assay_id", "role", "oligo_position", "original_symbol", "original_support",
                "observed_symbol", "observed_support", "affected_sequence_ids",
                "affected_sequence_count", "affected_fraction", "primer_3_prime_position",
            ),
            DegeneracyProposal: (
                "assay_id", "role", "original_sequence", "proposed_sequence", "status",
                "reason", "original_degeneracy", "proposed_degeneracy", "changed_positions",
                "binding_site_count", "original_exact_count", "original_exact_fraction",
                "proposed_exact_count", "proposed_exact_fraction",
            ),
            InclusivityResult: (
                "status", "evaluation_sequence_ids", "oligo_matches", "assay_results",
                "variations", "proposals", "oligo_matches_path", "assay_inclusivity_path",
                "oligo_variations_path", "degeneracy_proposals_path", "report_path",
            ),
        }
        for model, names in expected_fields.items():
            self.assertEqual(tuple(field.name for field in fields(model)), names)
            self.assertEqual(model.__slots__, names)

    def test_looks_up_assay_oligos_by_role(self):
        assay = self._assay()
        self.assertIs(_oligo_for_role(assay, "FORWARD"), assay.forward_primer)
        self.assertIs(_oligo_for_role(assay, "PROBE"), assay.probe)
        self.assertIs(_oligo_for_role(assay, "REVERSE"), assay.reverse_primer)

    def test_enumerates_exact_forward_hit_with_source_coordinates(self):
        hits = self._hits("TTACGTACTT", oligo=self._oligo("ACGT", 3, 6))
        hit = hits[0].public
        self.assertEqual((hit.source_start, hit.source_end), (3, 6))
        self.assertEqual(hit.mismatch_positions, ())
        self.assertTrue(hit.exact_match)
        self.assertTrue(hit.compatible)
        self.assertFalse(hit.selected)

    def test_preserves_a_first_base_mismatch_as_5_prime_and_compatible(self):
        hit = self._hits(
            "TTTACCGGTTAATT",
            oligo=self._oligo("AACCGGTTAA", 3, 12),
        )[0].public
        self.assertEqual(hit.mismatch_positions, (1,))
        self.assertFalse(hit.three_prime_mismatch)
        self.assertFalse(hit.probe_mismatch)
        self.assertTrue(hit.compatible)

    def test_keeps_an_internal_two_mismatch_primer_compatible(self):
        hit = self._hits(
            "ATTCGGTTAA",
            oligo=self._oligo("AACCGGTTAA", 1, 10),
        )[0].public
        self.assertEqual(hit.mismatch_positions, (2, 3))
        self.assertFalse(hit.three_prime_mismatch)
        self.assertTrue(hit.compatible)

    def test_marks_each_position_in_final_five_as_incompatible_3_prime_primer_mismatch(self):
        oligo = self._oligo("AACCGGTTAA", 1, 10)
        for position in range(6, 11):
            target = list(oligo.sequence)
            target[position - 1] = "C" if target[position - 1] != "C" else "A"
            hit = self._hits("".join(target), oligo=oligo)[0].public
            self.assertEqual(hit.mismatch_positions, (position,))
            self.assertTrue(hit.three_prime_mismatch)
            self.assertFalse(hit.compatible)

    def test_classifies_probe_mismatch_without_primer_3_prime_penalty(self):
        hit = self._hits(
            "ACTT",
            oligo=self._oligo("ACGT", 1, 4),
            role="PROBE",
        )[0].public
        self.assertEqual(hit.mismatch_positions, (3,))
        self.assertFalse(hit.three_prime_mismatch)
        self.assertTrue(hit.probe_mismatch)
        self.assertTrue(hit.compatible)

    def test_uses_conservative_iupac_target_subset_matching(self):
        exact = self._hits("AR", oligo=self._oligo("AN", 1, 2))[0].public
        mismatch = self._hits("AN", oligo=self._oligo("AR", 1, 2))[0].public
        self.assertTrue(exact.exact_match)
        self.assertEqual(mismatch.mismatch_positions, (2,))

    def test_compares_reverse_primer_in_its_synthesis_orientation(self):
        hit = self._hits(
            "TTGACTTT",
            oligo=self._oligo("AGTC", 3, 6),
            role="REVERSE",
        )[0]
        self.assertTrue(hit.public.exact_match)
        self.assertEqual(hit.target_in_synthesis_orientation, "AGTC")

    def test_converts_reverse_complement_oriented_interval_to_source_coordinates(self):
        hit = self._hits(
            "TTACGTACTTTT",
            oligo=self._oligo("ACGT", 3, 6),
            orientation="REVERSE_COMPLEMENT",
        )[0]
        self.assertEqual((hit.public.source_start, hit.public.source_end), (7, 10))
        self.assertEqual(hit.public.orientation, "REVERSE_COMPLEMENT")

    def test_clips_expanded_windows_to_record_boundaries(self):
        hit = self._hits(
            "ACGT",
            oligo=self._oligo("ACGT", 2, 5),
            config=InclusivityConfig(search_flank=10, max_hits_per_oligo=3),
        )[0].public
        self.assertEqual((hit.source_start, hit.source_end), (1, 4))
        self.assertEqual(hit.displacement, 1)

    def test_returns_no_hit_when_record_is_shorter_than_oligo(self):
        self.assertEqual(self._hits("ACG", oligo=self._oligo("ACGT", 1, 4)), ())

    def test_ranks_internal_mismatches_before_3_prime_mismatches_then_applies_cap(self):
        hits = self._hits(
            "AAAAATCTAAAAA",
            oligo=self._oligo("AAAAAA", 1, 6),
            config=InclusivityConfig(search_flank=7, max_hits_per_oligo=2),
        )
        self.assertEqual([hit.public.source_start for hit in hits], [8, 1])
        self.assertEqual([hit.public.hit_rank for hit in hits], [1, 2])
        self.assertEqual([hit.public.three_prime_mismatch for hit in hits], [False, True])

    def test_repeated_equally_scoring_hits_use_oriented_start_then_target_text(self):
        hits = self._hits(
            "ACGTACGTACGT",
            oligo=self._oligo("ACGT", 5, 8),
            config=InclusivityConfig(search_flank=8, max_hits_per_oligo=2),
        )
        self.assertEqual([hit.public.source_start for hit in hits], [5, 1])
        self.assertEqual([hit.public.hit_rank for hit in hits], [1, 2])

    def test_invalid_iupac_input_reports_assay_record_and_role_without_sequence(self):
        with self.assertRaisesRegex(InclusivityError, r"a1.*s1.*FORWARD.*position 3"):
            self._hits("ACX", oligo=self._oligo("ACG", 1, 3))


class InclusivityGeometryTests(unittest.TestCase):
    def _oligo(self, sequence: str, start: int, end: int | None = None) -> DesignedOligo:
        end = end if end is not None else start + len(sequence) - 1
        return DesignedOligo(
            sequence=sequence,
            reference_start=start,
            reference_end=end,
            length=len(sequence),
            tm=60.0,
            gc_percent=50.0,
            penalty=None,
            metrics=(),
        )

    def _assay(self, assay_id: str = "a1", *, product_size: int = 16) -> AssayCandidate:
        return AssayCandidate(
            assay_id=assay_id,
            region_id="r1",
            primer3_index=0,
            forward_primer=self._oligo("ACGT", 3, 6),
            probe=self._oligo("TTAA", 9, 12),
            reverse_primer=self._oligo("AGTC", 15, 18),
            product_size=product_size,
            pair_penalty=None,
            metrics=(),
        )

    def _record(self, sequence_id: str = "s1", sequence: str = "TTACGTCCTTAAGGGACTAA") -> LocalSequenceRecord:
        return LocalSequenceRecord(sequence_id, sequence)

    def test_selects_complete_forward_geometry_with_strict_order_and_size(self):
        selected = _select_binding(
            self._record(), self._assay(), InclusivityConfig(search_flank=2)
        )
        self.assertTrue(selected.geometry_found)
        self.assertEqual(selected.orientation, "FORWARD")
        self.assertEqual((selected.source_amplicon_start, selected.source_amplicon_end), (3, 18))
        self.assertEqual(selected.amplicon_size, 16)
        self.assertTrue(selected.original_compatible)
        self.assertTrue(all(hit.public.selected for hit in (selected.forward, selected.probe, selected.reverse)))

    def test_rejects_triplets_when_probe_overlaps_a_primer(self):
        record = self._record(sequence="TTACGTTAAGGGGGGACTAA")
        selected = _select_binding(
            record, self._assay(), InclusivityConfig(search_flank=3, max_hits_per_oligo=1)
        )
        self.assertFalse(selected.geometry_found)
        self.assertIsNone(selected.orientation)
        self.assertFalse(selected.original_compatible)

    def test_rejects_triplets_outside_the_amplicon_size_delta(self):
        selected = _select_binding(
            self._record(sequence="TTACGTCCTTAAGGGGGGGACTAA"),
            self._assay(),
            InclusivityConfig(search_flank=4, max_hits_per_oligo=1, max_amplicon_size_delta=1),
        )
        self.assertFalse(selected.geometry_found)
        self.assertIsNone(selected.amplicon_size)
        self.assertFalse(selected.original_compatible)

    def test_selects_reverse_complement_orientation_and_ascending_source_interval(self):
        selected = _select_binding(
            self._record(sequence="TTAGTCCCTTAAGGACGTAA"),
            self._assay(),
            InclusivityConfig(search_flank=2),
        )
        self.assertTrue(selected.geometry_found)
        self.assertEqual(selected.orientation, "REVERSE_COMPLEMENT")
        self.assertEqual((selected.source_amplicon_start, selected.source_amplicon_end), (3, 18))
        self.assertEqual(selected.amplicon_size, 16)

    def test_breaks_an_exact_candidate_tie_with_forward_orientation(self):
        sequence = list("C" * 100)
        for start, target in (
            (3, "ACGT"),
            (9, "TTAA"),
            (15, "GACT"),
            (83, "AGTC"),
            (89, "TTAA"),
            (95, "ACGT"),
        ):
            sequence[start - 1:start - 1 + len(target)] = target
        selected = _select_binding(
            self._record(sequence="".join(sequence)),
            self._assay(),
            InclusivityConfig(search_flank=0),
        )
        self.assertEqual(selected.orientation, "FORWARD")
        self.assertEqual(selected.forward.public.source_start, 3)

    def test_breaks_reverse_fallback_ties_by_raw_oriented_target_segment(self):
        def hit(target_in_synthesis_orientation: str) -> _Hit:
            return _Hit(
                public=OligoMatch(
                    assay_id="a1",
                    sequence_id="s1",
                    role="REVERSE",
                    orientation="FORWARD",
                    hit_rank=1,
                    source_start=1,
                    source_end=3,
                    expected_start=1,
                    expected_end=3,
                    displacement=0,
                    mismatch_positions=(),
                    mismatch_count=0,
                    exact_match=True,
                    three_prime_mismatch=False,
                    probe_mismatch=False,
                    compatible=True,
                    selected=False,
                ),
                oriented_start=1,
                oriented_end=3,
                oriented_target_segment={"AAA": "TTT", "CAA": "TTG"}[target_in_synthesis_orientation],
                target_in_synthesis_orientation=target_in_synthesis_orientation,
            )

        synthesis_first = hit("AAA")
        raw_first = hit("CAA")
        selected = _best_fallback(
            (synthesis_first, raw_first), "REVERSE", InclusivityConfig()
        )
        self.assertIs(selected, raw_first)
        self.assertEqual(selected.oriented_target_segment, "TTG")

    def test_evaluates_every_evaluation_record_once_in_assay_major_order(self):
        records = (self._record("s1"), self._record("s2"))
        selected = _evaluate_original(
            records,
            EvaluationSet(("s1", "s2")),
            (self._assay("a1"), self._assay("a2")),
            InclusivityConfig(search_flank=2),
        )
        self.assertEqual(
            [(item.assay.assay_id, item.sequence_id) for item in selected],
            [("a1", "s1"), ("a1", "s2"), ("a2", "s1"), ("a2", "s2")],
        )

    def test_rejects_invalid_evaluation_inputs_without_leaking_sequences(self):
        cases = (
            ((self._record("s1"),), EvaluationSet(("",)), r"Evaluation Set"),
            ((self._record("s1"), self._record("s1")), EvaluationSet(("s1",)), r"record"),
            ((self._record("s2"),), EvaluationSet(("s1",)), r"record"),
            ((self._record("s2"), self._record("s1")), EvaluationSet(("s1", "s2")), r"order"),
        )
        for records, evaluation_set, message in cases:
            with self.subTest(records=records, evaluation_set=evaluation_set):
                with self.assertRaisesRegex(InclusivityError, message) as raised:
                    _evaluate_original(
                        records,
                        evaluation_set,
                        (self._assay(),),
                        InclusivityConfig(search_flank=2),
                    )
                self.assertNotIn("TTACGTCCTTAAGGGACTAA", str(raised.exception))

    def test_accepts_an_empty_evaluation_set(self):
        self.assertEqual(
            _evaluate_original((), EvaluationSet(()), (self._assay(),), InclusivityConfig(search_flank=2)),
            (),
        )

    def test_rejects_invalid_assay_iupac_with_role_position_context_without_sequence(self):
        invalid = self._assay()
        invalid = AssayCandidate(
            assay_id=invalid.assay_id,
            region_id=invalid.region_id,
            primer3_index=invalid.primer3_index,
            forward_primer=self._oligo("ACXT", 3, 6),
            probe=invalid.probe,
            reverse_primer=invalid.reverse_primer,
            product_size=invalid.product_size,
            pair_penalty=invalid.pair_penalty,
            metrics=invalid.metrics,
        )
        with self.assertRaisesRegex(InclusivityError, r"a1.*FORWARD.*position 3") as raised:
            _evaluate_original(
                (self._record(),), EvaluationSet(("s1",)), (invalid,), InclusivityConfig(search_flank=2)
            )
        self.assertNotIn("ACXT", str(raised.exception))


if __name__ == "__main__":
    unittest.main()
