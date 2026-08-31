from dataclasses import replace
import unittest

from qpcr_pipeline.config import RankingConfig, RankingWeights
from qpcr_pipeline.ranking import RankingError, rank_assays
from ranking_fixtures import (
    make_amplicon,
    make_assay,
    make_inclusivity_result,
    make_primer_result,
    make_proposal,
    make_region,
    make_specificity_result,
)


class RankingScoringTests(unittest.TestCase):
    def _fraction_inputs(self, compatible_count: int, *, pair_penalty=1.0, proposals=()):
        sequence_ids = tuple(f"s{i}" for i in range(10))
        primer = make_primer_result(
            assays=(make_assay(pair_penalty=pair_penalty),),
            candidates=(make_region(),),
        )
        compatibility = {
            ("a1", sequence_id): index < compatible_count
            for index, sequence_id in enumerate(sequence_ids)
        }
        inclusivity = make_inclusivity_result(
            primer,
            compatibility=compatibility,
            sequence_ids=sequence_ids,
            proposals=proposals,
        )
        return primer, inclusivity

    def test_named_components_and_weighted_score_use_absolute_formulas(self):
        proposal = make_proposal(
            role="FORWARD",
            status="ACCEPTED",
            original_degeneracy=1,
            proposed_degeneracy=2,
        )
        primer, inclusivity = self._fraction_inputs(9, proposals=(proposal,))
        specificity = make_specificity_result(
            primer, hit_totals={("off", "a1", "FORWARD"): 2}
        )

        item = rank_assays(
            primer, inclusivity, specificity, RankingConfig(enabled=True)
        )[0]

        self.assertAlmostEqual(item.components.inclusivity, 0.9)
        self.assertAlmostEqual(item.components.specificity, 0.96)
        self.assertAlmostEqual(item.components.conservation, 1.0)
        self.assertAlmostEqual(item.components.primer3_quality, 0.5)
        self.assertAlmostEqual(
            item.components.robustness, (0.5 + 1.0 + 1.0) / 3.0
        )
        self.assertAlmostEqual(item.final_score, 88.83333333333333)
        self.assertEqual(item.score_status, "COMPLETE")
        self.assertEqual(item.classification, "REVIEW")

    def test_specificity_component_uses_full_retention_totals_and_risk_geometry(self):
        primer, inclusivity = self._fraction_inputs(10)
        no_hits = make_specificity_result(primer)
        one_hit = make_specificity_result(
            primer, hit_totals={("off", "a1", "FORWARD"): 1}
        )
        twenty_five = make_specificity_result(
            primer, hit_totals={("off", "a1", "FORWARD"): 25}
        )
        twenty_five = replace(
            twenty_five,
            retention=tuple(
                replace(row, retained_hit_count=1, truncated=True)
                if row.role == "FORWARD"
                else row
                for row in twenty_five.retention
            ),
        )
        plausible = make_specificity_result(
            primer, amplicons=(make_amplicon(detectable=False),)
        )
        detectable = make_specificity_result(
            primer, amplicons=(make_amplicon(detectable=True),)
        )

        cases = (
            (no_hits, 1.00),
            (one_hit, 0.98),
            (twenty_five, 0.80),
            (plausible, 0.40),
            (detectable, 0.00),
        )
        for specificity, expected in cases:
            with self.subTest(expected=expected):
                item = rank_assays(
                    primer,
                    inclusivity,
                    specificity,
                    RankingConfig(enabled=True),
                )[0]
                self.assertAlmostEqual(item.components.specificity, expected)

    def test_missing_pair_penalty_makes_score_incomplete_and_forces_review(self):
        primer, inclusivity = self._fraction_inputs(10, pair_penalty=None)
        item = rank_assays(
            primer,
            inclusivity,
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )[0]
        self.assertEqual(item.score_status, "INCOMPLETE")
        self.assertIsNone(item.final_score)
        self.assertIsNone(item.components.primer3_quality)
        self.assertEqual(item.classification, "REVIEW")
        self.assertIn("EVIDENCE_INCOMPLETE", {reason.code for reason in item.reasons})

    def test_zero_weight_does_not_hide_missing_component(self):
        primer, inclusivity = self._fraction_inputs(10, pair_penalty=None)
        config = RankingConfig(
            enabled=True,
            weights=RankingWeights(
                inclusivity=0.40,
                specificity=0.30,
                conservation=0.20,
                primer3_quality=0.00,
                robustness=0.10,
            ),
        )
        item = rank_assays(
            primer, inclusivity, make_specificity_result(primer), config
        )[0]
        self.assertEqual(item.score_status, "INCOMPLETE")
        self.assertIsNone(item.final_score)

    def test_empty_evaluation_set_is_incomplete_without_division_by_zero(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer, sequence_ids=())
        item = rank_assays(
            primer,
            inclusivity,
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )[0]
        self.assertIsNone(item.components.inclusivity)
        self.assertIsNone(item.components.robustness)
        self.assertEqual(item.score_status, "INCOMPLETE")
        self.assertIsNone(item.final_score)
        self.assertEqual(item.classification, "REVIEW")

    def test_invalid_present_scoring_metrics_fail_explicitly(self):
        base_assay = make_assay()
        invalid_regions = (
            replace(make_region(), mean_conservation=1.1),
            replace(make_region(), mean_entropy_bits=-0.1),
        )
        for region in invalid_regions:
            primer = make_primer_result(
                assays=(base_assay,), candidates=(region,)
            )
            with self.subTest(region=region):
                with self.assertRaises(RankingError):
                    rank_assays(
                        primer,
                        make_inclusivity_result(primer),
                        make_specificity_result(primer),
                        RankingConfig(enabled=True),
                    )

        primer = make_primer_result(
            assays=(make_assay(pair_penalty=-1.0),), candidates=(make_region(),)
        )
        with self.assertRaises(RankingError):
            rank_assays(
                primer,
                make_inclusivity_result(primer),
                make_specificity_result(primer),
                RankingConfig(enabled=True),
            )

        primer = make_primer_result()
        invalid_proposal = make_proposal(
            original_degeneracy=4, proposed_degeneracy=2
        )
        with self.assertRaises(RankingError):
            rank_assays(
                primer,
                make_inclusivity_result(primer, proposals=(invalid_proposal,)),
                make_specificity_result(primer),
                RankingConfig(enabled=True),
            )


if __name__ == "__main__":
    unittest.main()
