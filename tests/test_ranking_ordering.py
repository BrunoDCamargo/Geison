import unittest

from qpcr_pipeline.config import RankingConfig, RankingWeights
from qpcr_pipeline.ranking import rank_assays
from ranking_fixtures import (
    make_assay,
    make_inclusivity_result,
    make_primer_result,
    make_region,
    make_specificity_result,
)


class RankingOrderingTests(unittest.TestCase):
    @staticmethod
    def _compatibility(sequence_ids, counts):
        result = {}
        for assay_id, compatible_count in counts.items():
            for index, sequence_id in enumerate(sequence_ids):
                result[(assay_id, sequence_id)] = index < compatible_count
        return result

    def test_class_priority_beats_aggregate_score(self):
        assays = (
            make_assay("pass", "r-pass", primer3_index=2, pair_penalty=20.0),
            make_assay("review", "r-review", primer3_index=1, pair_penalty=10.0),
            make_assay("risk", "r-risk", primer3_index=0, pair_penalty=0.0),
        )
        candidates = (
            make_region("r-pass"),
            make_region("r-review"),
            make_region("r-risk"),
        )
        primer = make_primer_result(assays=assays, candidates=candidates)
        sequence_ids = tuple(f"s{i}" for i in range(10))
        inclusivity = make_inclusivity_result(
            primer,
            compatibility=self._compatibility(
                sequence_ids, {"pass": 10, "review": 9, "risk": 8}
            ),
            sequence_ids=sequence_ids,
        )
        ranked = rank_assays(
            primer,
            inclusivity,
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )
        self.assertEqual(
            [item.classification for item in ranked],
            ["IN SILICO PASS", "REVIEW", "HIGH_RISK"],
        )
        self.assertEqual([item.assay_id for item in ranked], ["pass", "review", "risk"])
        self.assertEqual([item.rank for item in ranked], [1, 2, 3])

    def test_complete_score_sorts_before_incomplete_within_class(self):
        assays = (
            make_assay("complete", "r-complete", pair_penalty=1.0),
            make_assay("incomplete", "r-incomplete", pair_penalty=None),
        )
        primer = make_primer_result(
            assays=assays,
            candidates=(make_region("r-complete"), make_region("r-incomplete")),
        )
        sequence_ids = tuple(f"s{i}" for i in range(10))
        inclusivity = make_inclusivity_result(
            primer,
            compatibility=self._compatibility(
                sequence_ids, {"complete": 9, "incomplete": 9}
            ),
            sequence_ids=sequence_ids,
        )
        ranked = rank_assays(
            primer,
            inclusivity,
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )
        self.assertEqual([item.assay_id for item in ranked], ["complete", "incomplete"])
        self.assertEqual(
            [item.score_status for item in ranked], ["COMPLETE", "INCOMPLETE"]
        )

    def test_final_score_descending_orders_same_class(self):
        assays = (
            make_assay("lower", "r-lower", pair_penalty=5.0),
            make_assay("higher", "r-higher", pair_penalty=0.0),
        )
        primer = make_primer_result(
            assays=assays,
            candidates=(make_region("r-lower"), make_region("r-higher")),
        )
        ranked = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )
        self.assertEqual([item.assay_id for item in ranked], ["higher", "lower"])
        self.assertGreater(ranked[0].final_score, ranked[1].final_score)

    def test_inclusivity_breaks_equal_score_when_its_weight_is_zero(self):
        assays = (
            make_assay("eight", "r-eight"),
            make_assay("nine", "r-nine"),
        )
        primer = make_primer_result(
            assays=assays,
            candidates=(make_region("r-eight"), make_region("r-nine")),
        )
        sequence_ids = tuple(f"s{i}" for i in range(10))
        inclusivity = make_inclusivity_result(
            primer,
            compatibility=self._compatibility(
                sequence_ids, {"eight": 8, "nine": 9}
            ),
            sequence_ids=sequence_ids,
        )
        config = RankingConfig(
            enabled=True,
            min_inclusivity_for_pass=1.0,
            min_inclusivity_before_high_risk=0.0,
            weights=RankingWeights(
                inclusivity=0.0,
                specificity=0.40,
                conservation=0.30,
                primer3_quality=0.15,
                robustness=0.15,
            ),
        )
        ranked = rank_assays(
            primer, inclusivity, make_specificity_result(primer), config
        )
        self.assertAlmostEqual(ranked[0].final_score, ranked[1].final_score)
        self.assertEqual([item.assay_id for item in ranked], ["nine", "eight"])

    def test_pair_penalty_breaks_equal_score_when_primer3_weight_is_zero(self):
        assays = (
            make_assay("penalty-high", "r-high", pair_penalty=4.0),
            make_assay("penalty-low", "r-low", pair_penalty=0.5),
        )
        primer = make_primer_result(
            assays=assays,
            candidates=(make_region("r-high"), make_region("r-low")),
        )
        config = RankingConfig(
            enabled=True,
            weights=RankingWeights(
                inclusivity=0.40,
                specificity=0.30,
                conservation=0.20,
                primer3_quality=0.0,
                robustness=0.10,
            ),
        )
        ranked = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(primer),
            config,
        )
        self.assertAlmostEqual(ranked[0].final_score, ranked[1].final_score)
        self.assertEqual(
            [item.assay_id for item in ranked], ["penalty-low", "penalty-high"]
        )

    def test_primer3_index_then_assay_id_are_stable_final_ties(self):
        assays = (
            make_assay("b", "r-b", primer3_index=1),
            make_assay("a", "r-a", primer3_index=1),
            make_assay("c", "r-c", primer3_index=0),
        )
        primer = make_primer_result(
            assays=assays,
            candidates=(make_region("r-b"), make_region("r-a"), make_region("r-c")),
        )
        ranked = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )
        self.assertEqual([item.assay_id for item in ranked], ["c", "a", "b"])
        self.assertEqual([item.rank for item in ranked], [1, 2, 3])
        self.assertEqual({item.assay_id for item in ranked}, {"a", "b", "c"})

    def test_score_is_absolute_not_relative_to_other_assays(self):
        assay = make_assay("stable", "r-stable", pair_penalty=1.0)
        single_primer = make_primer_result(
            assays=(assay,), candidates=(make_region("r-stable"),)
        )
        single = rank_assays(
            single_primer,
            make_inclusivity_result(single_primer),
            make_specificity_result(single_primer),
            RankingConfig(enabled=True),
        )[0]

        multi_primer = make_primer_result(
            assays=(assay, make_assay("other", "r-other", pair_penalty=0.0)),
            candidates=(make_region("r-stable"), make_region("r-other")),
        )
        multi = rank_assays(
            multi_primer,
            make_inclusivity_result(multi_primer),
            make_specificity_result(multi_primer),
            RankingConfig(enabled=True),
        )
        stable = next(item for item in multi if item.assay_id == "stable")
        self.assertAlmostEqual(single.final_score, stable.final_score)


if __name__ == "__main__":
    unittest.main()
