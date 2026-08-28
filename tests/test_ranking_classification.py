from dataclasses import replace
import unittest

from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.inclusivity import OligoMatch, OligoVariation
from qpcr_pipeline.ranking import RankingError, classify_assays
from qpcr_pipeline.specificity import HitRetentionSummary
from ranking_fixtures import (
    make_amplicon,
    make_assay,
    make_hit,
    make_inclusivity_result,
    make_primer_result,
    make_proposal,
    make_region,
    make_specificity_result,
)


class RankingClassificationTests(unittest.TestCase):
    def _classify_fraction(self, compatible_count: int, config=None):
        sequence_ids = tuple(f"s{i}" for i in range(10))
        primer = make_primer_result()
        compatibility = {
            ("a1", sequence_id): index < compatible_count
            for index, sequence_id in enumerate(sequence_ids)
        }
        inclusivity = make_inclusivity_result(
            primer, compatibility=compatibility, sequence_ids=sequence_ids
        )
        specificity = make_specificity_result(primer)
        return classify_assays(
            primer, inclusivity, specificity, config or RankingConfig(enabled=True)
        )[0]

    def test_default_inclusivity_boundaries_classify_before_score(self):
        passed = self._classify_fraction(10)
        review = self._classify_fraction(9)
        high_risk = self._classify_fraction(8)

        self.assertEqual(passed.classification, "IN SILICO PASS")
        self.assertEqual(passed.reasons, ())
        self.assertEqual(review.classification, "REVIEW")
        self.assertIn("INCLUSIVITY_BELOW_PASS", {r.code for r in review.reasons})
        self.assertEqual(high_risk.classification, "HIGH_RISK")
        self.assertIn(
            "INCLUSIVITY_BELOW_MINIMUM", {r.code for r in high_risk.reasons}
        )

    def test_configurable_inclusivity_thresholds_respect_boundary(self):
        config = RankingConfig(
            enabled=True,
            min_inclusivity_for_pass=0.95,
            min_inclusivity_before_high_risk=0.80,
        )
        review = self._classify_fraction(8, config)
        high_risk = self._classify_fraction(7, config)
        self.assertEqual(review.classification, "REVIEW")
        self.assertEqual(high_risk.classification, "HIGH_RISK")

    def test_detectable_off_target_forces_high_risk(self):
        primer = make_primer_result()
        result = classify_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(
                primer, amplicons=(make_amplicon(detectable=True),)
            ),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "HIGH_RISK")
        self.assertIn("DETECTABLE_OFF_TARGET", {r.code for r in result.reasons})

    def test_primer_amplicon_without_probe_forces_review(self):
        primer = make_primer_result()
        result = classify_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(
                primer, amplicons=(make_amplicon(detectable=False),)
            ),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "REVIEW")
        self.assertIn(
            "PLAUSIBLE_OFF_TARGET_AMPLICON", {r.code for r in result.reasons}
        )

    def test_isolated_hits_are_advisory_only(self):
        primer = make_primer_result()
        specificity = make_specificity_result(
            primer, hit_totals={("off", "a1", "FORWARD"): 3}
        )
        result = classify_assays(
            primer,
            make_inclusivity_result(primer),
            specificity,
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "IN SILICO PASS")
        self.assertEqual({r.code for r in result.reasons}, {"ISOLATED_OFF_TARGET_HITS"})

    def test_skipped_evidence_forces_review_with_deduplicated_codes(self):
        primer = make_primer_result()
        result = classify_assays(
            primer,
            make_inclusivity_result(primer, status="SKIPPED"),
            make_specificity_result(primer, status="SKIPPED"),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "REVIEW")
        self.assertEqual(
            {reason.code for reason in result.reasons},
            {
                "EVIDENCE_INCOMPLETE",
                "INCLUSIVITY_EVIDENCE_MISSING",
                "SPECIFICITY_EVIDENCE_MISSING",
            },
        )
        self.assertEqual(
            sum(reason.code == "EVIDENCE_INCOMPLETE" for reason in result.reasons), 1
        )

    def test_available_high_risk_evidence_wins_over_missing_source(self):
        primer = make_primer_result()
        result = classify_assays(
            primer,
            make_inclusivity_result(primer, status="SKIPPED"),
            make_specificity_result(
                primer, amplicons=(make_amplicon(detectable=True),)
            ),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "HIGH_RISK")
        self.assertIn("EVIDENCE_INCOMPLETE", {r.code for r in result.reasons})
        self.assertIn("DETECTABLE_OFF_TARGET", {r.code for r in result.reasons})

    def test_degeneracy_proposals_are_aggregated_advisories(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(
            primer,
            proposals=(
                make_proposal(role="FORWARD", status="ACCEPTED"),
                make_proposal(role="PROBE", status="REJECTED"),
            ),
        )
        result = classify_assays(
            primer,
            inclusivity,
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(result.classification, "IN SILICO PASS")
        by_code = {reason.code: reason for reason in result.reasons}
        self.assertEqual(
            set(by_code), {"IUPAC_PROPOSAL_ACCEPTED", "IUPAC_PROPOSAL_REJECTED"}
        )
        self.assertIn("FORWARD", dict(by_code["IUPAC_PROPOSAL_ACCEPTED"].evidence)["roles"])
        self.assertIn("PROBE", dict(by_code["IUPAC_PROPOSAL_REJECTED"].evidence)["roles"])

    def test_reason_order_is_deterministic_when_input_order_changes(self):
        primer = make_primer_result()
        proposals = (
            make_proposal(role="FORWARD", status="ACCEPTED"),
            make_proposal(role="PROBE", status="REJECTED"),
        )
        first = classify_assays(
            primer,
            make_inclusivity_result(primer, proposals=proposals),
            make_specificity_result(
                primer, hit_totals={("off", "a1", "FORWARD"): 1}
            ),
            RankingConfig(enabled=True),
        )[0]
        second = classify_assays(
            primer,
            make_inclusivity_result(primer, proposals=tuple(reversed(proposals))),
            make_specificity_result(
                primer, hit_totals={("off", "a1", "FORWARD"): 1}
            ),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(
            tuple((r.severity, r.source, r.code) for r in first.reasons),
            tuple((r.severity, r.source, r.code) for r in second.reasons),
        )


class RankingEvidenceIntegrityTests(unittest.TestCase):
    def setUp(self):
        self.config = RankingConfig(enabled=True)

    def _assert_invalid(self, primer, inclusivity, specificity):
        with self.assertRaises(RankingError):
            classify_assays(primer, inclusivity, specificity, self.config)

    def test_rejects_duplicate_assay_id_and_missing_region(self):
        duplicate = make_primer_result(
            assays=(make_assay(), make_assay()), candidates=(make_region(),)
        )
        self._assert_invalid(
            duplicate,
            make_inclusivity_result(duplicate),
            make_specificity_result(duplicate),
        )

        missing_region = make_primer_result(
            assays=(make_assay(region_id="missing"),), candidates=(make_region(),)
        )
        self._assert_invalid(
            missing_region,
            make_inclusivity_result(missing_region),
            make_specificity_result(missing_region),
        )

    def test_rejects_missing_and_duplicate_inclusivity_matrix_rows(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer, sequence_ids=("s1", "s2"))
        missing = replace(inclusivity, assay_results=inclusivity.assay_results[:-1])
        duplicate = replace(
            inclusivity,
            assay_results=inclusivity.assay_results + (inclusivity.assay_results[0],),
        )
        for invalid in (missing, duplicate):
            with self.subTest(row_count=len(invalid.assay_results)):
                self._assert_invalid(primer, invalid, make_specificity_result(primer))

    def test_rejects_unknown_assays_in_inclusivity_detail_rows(self):
        primer = make_primer_result()
        base = make_inclusivity_result(primer)
        unknown_match = OligoMatch(
            assay_id="unknown",
            sequence_id="s1",
            role="FORWARD",
            orientation="FORWARD",
            hit_rank=1,
            source_start=1,
            source_end=4,
            expected_start=1,
            expected_end=4,
            displacement=0,
            mismatch_positions=(),
            mismatch_count=0,
            exact_match=True,
            three_prime_mismatch=False,
            probe_mismatch=False,
            compatible=True,
            selected=True,
        )
        unknown_variation = OligoVariation(
            assay_id="unknown",
            role="FORWARD",
            oligo_position=1,
            original_symbol="A",
            original_support=("A",),
            observed_symbol="G",
            observed_support=("G",),
            affected_sequence_ids=("s1",),
            affected_sequence_count=1,
            affected_fraction=1.0,
            primer_3_prime_position=False,
        )
        invalid_results = (
            replace(base, oligo_matches=(unknown_match,)),
            replace(base, variations=(unknown_variation,)),
            replace(base, proposals=(make_proposal(assay_id="unknown"),)),
        )
        for invalid in invalid_results:
            with self.subTest(invalid=invalid):
                self._assert_invalid(primer, invalid, make_specificity_result(primer))

    def test_rejects_duplicate_degeneracy_proposal_for_assay_role(self):
        primer = make_primer_result()
        proposal = make_proposal()
        inclusivity = make_inclusivity_result(
            primer, proposals=(proposal, replace(proposal, proposed_sequence="AYGT"))
        )
        self._assert_invalid(primer, inclusivity, make_specificity_result(primer))

    def test_rejects_specificity_assay_count_and_retention_matrix_errors(self):
        primer = make_primer_result()
        base = make_specificity_result(primer)
        invalid_results = (
            replace(base, assay_count=2),
            replace(base, retention=base.retention[:-1]),
            replace(base, retention=base.retention + (base.retention[0],)),
            replace(
                base,
                retention=(
                    replace(base.retention[0], assay_id="unknown"),
                    *base.retention[1:],
                ),
            ),
        )
        inclusivity = make_inclusivity_result(primer)
        for invalid in invalid_results:
            with self.subTest(invalid=invalid):
                self._assert_invalid(primer, inclusivity, invalid)

    def test_rejects_invalid_retention_counts(self):
        primer = make_primer_result()
        base = make_specificity_result(primer)
        invalid_rows = (
            replace(base.retention[0], total_hit_count=-1),
            replace(base.retention[0], retained_hit_count=-1),
            replace(base.retention[0], total_hit_count=1, retained_hit_count=2),
        )
        for row in invalid_rows:
            specificity = replace(base, retention=(row, *base.retention[1:]))
            with self.subTest(row=row):
                self._assert_invalid(
                    primer, make_inclusivity_result(primer), specificity
                )

    def test_rejects_unknown_assays_in_specificity_hits_and_amplicons(self):
        primer = make_primer_result()
        base = make_specificity_result(primer)
        invalid_results = (
            replace(base, hits=(make_hit(assay_id="unknown"),)),
            replace(base, amplicons=(make_amplicon(assay_id="unknown"),)),
        )
        for invalid in invalid_results:
            with self.subTest(invalid=invalid):
                self._assert_invalid(
                    primer, make_inclusivity_result(primer), invalid
                )

    def test_rejects_duplicate_dataset_names(self):
        primer = make_primer_result()
        base = make_specificity_result(primer)
        specificity = replace(base, dataset_names=("off", "off"))
        self._assert_invalid(
            primer, make_inclusivity_result(primer), specificity
        )

    def test_rejects_non_complete_primer_design(self):
        primer = make_primer_result(assays=(), candidates=(), status="SKIPPED")
        self._assert_invalid(
            primer,
            make_inclusivity_result(primer, status="SKIPPED"),
            make_specificity_result(primer, status="SKIPPED"),
        )


if __name__ == "__main__":
    unittest.main()
