import unittest

from qpcr_pipeline.assay_report_html import render_assay_report_html
from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.ranking import rank_assays
from ranking_fixtures import (
    make_amplicon,
    make_inclusivity_result,
    make_primer_result,
    make_proposal,
    make_specificity_result,
)


class AssayReportHtmlTests(unittest.TestCase):
    def test_report_consolidates_design_evidence_classification_and_score(self):
        primer = make_primer_result()
        proposal = make_proposal(role="FORWARD", status="ACCEPTED")
        inclusivity = make_inclusivity_result(primer, proposals=(proposal,))
        specificity = make_specificity_result(
            primer,
            dataset_names=("human",),
            hit_totals={("human", "a1", "FORWARD"): 2},
        )
        ranked = rank_assays(
            primer, inclusivity, specificity, RankingConfig(enabled=True)
        )

        html = render_assay_report_html(
            target_name="target",
            primer_design=primer,
            inclusivity=inclusivity,
            specificity=specificity,
            assays=ranked,
        )

        for expected in (
            "a1",
            "IN SILICO PASS",
            "Rank",
            "Final score",
            "inclusivity",
            "specificity",
            "conservation",
            "primer3_quality",
            "robustness",
            "ACGT",
            "TTAA",
            "AGTC",
            "Product size",
            "Pair penalty",
            "Mean conservation",
            "Minimum conservation",
            "Mean coverage",
            "Mean gap frequency",
            "Mean entropy",
            "Original compatible",
            "1 / 1",
            "ARGT",
            "ACCEPTED",
            "human",
            "Compatible hits",
        ):
            with self.subTest(expected=expected):
                self.assertIn(expected, html)

    def test_report_escapes_dynamic_strings_and_has_no_remote_resources_or_script(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        ranked = rank_assays(
            primer, inclusivity, specificity, RankingConfig(enabled=True)
        )
        html = render_assay_report_html(
            target_name='<img src=x onerror="boom">',
            primer_design=primer,
            inclusivity=inclusivity,
            specificity=specificity,
            assays=ranked,
        )
        self.assertNotIn('<img src=x onerror="boom">', html)
        self.assertIn("&lt;img", html)
        lowered = html.lower()
        self.assertNotIn("http://", lowered)
        self.assertNotIn("https://", lowered)
        self.assertNotIn("<script", lowered)

    def test_zero_assays_has_visible_empty_state(self):
        primer = make_primer_result(assays=(), candidates=())
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        html = render_assay_report_html(
            target_name="target",
            primer_design=primer,
            inclusivity=inclusivity,
            specificity=specificity,
            assays=(),
        )
        self.assertIn("No assay candidates", html)


if __name__ == "__main__":
    unittest.main()
