import json
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.ranking import RankingError, evaluate_ranking
from ranking_fixtures import (
    make_assay,
    make_inclusivity_result,
    make_primer_result,
    make_region,
    make_specificity_result,
)


class RankingArtifactTests(unittest.TestCase):
    def test_disabled_removes_stale_ranking_tsv_but_preserves_root_conservation_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir)
            ranking_dir = output / "ranking"
            ranking_dir.mkdir()
            (ranking_dir / "assay_ranking.tsv").write_text("stale\n", encoding="utf-8")
            root_report = output / "report.html"
            root_report.write_text("conservation report", encoding="utf-8")

            result = evaluate_ranking(
                object(),
                object(),
                object(),
                RankingConfig(enabled=False),
                output,
                target_name="target",
            )

            self.assertEqual(result.status, "SKIPPED")
            self.assertIsNone(result.ranking_tsv_path)
            self.assertIsNone(result.html_report_path)
            self.assertFalse((ranking_dir / "assay_ranking.tsv").exists())
            self.assertEqual(root_report.read_text(encoding="utf-8"), "conservation report")
            report = json.loads(result.ranking_report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["status"], "SKIPPED")
            self.assertEqual(
                {path.name for path in ranking_dir.iterdir()}, {"ranking_report.json"}
            )

    def test_enabled_failure_clears_all_stale_ranking_owned_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir)
            ranking_dir = output / "ranking"
            ranking_dir.mkdir()
            for path in (
                ranking_dir / "assay_ranking.tsv",
                ranking_dir / "ranking_report.json",
                output / "report.html",
            ):
                path.write_text("stale", encoding="utf-8")

            primer = make_primer_result(assays=(), candidates=(), status="SKIPPED")
            with self.assertRaises(RankingError):
                evaluate_ranking(
                    primer,
                    make_inclusivity_result(primer, status="SKIPPED"),
                    make_specificity_result(primer, status="SKIPPED"),
                    RankingConfig(enabled=True),
                    output,
                    target_name="target",
                )

            self.assertFalse((ranking_dir / "assay_ranking.tsv").exists())
            self.assertFalse((ranking_dir / "ranking_report.json").exists())
            self.assertFalse((output / "report.html").exists())

    def test_enabled_success_publishes_complete_auditable_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir)
            primer = make_primer_result()
            inclusivity = make_inclusivity_result(primer)
            specificity = make_specificity_result(primer)

            result = evaluate_ranking(
                primer,
                inclusivity,
                specificity,
                RankingConfig(enabled=True),
                output,
                target_name="target",
            )

            self.assertEqual(result.status, "COMPLETE")
            self.assertEqual(result.ranking_tsv_path, output / "ranking" / "assay_ranking.tsv")
            self.assertEqual(result.ranking_report_path, output / "ranking" / "ranking_report.json")
            self.assertEqual(result.html_report_path, output / "report.html")
            self.assertTrue(result.ranking_tsv_path.exists())
            self.assertTrue(result.ranking_report_path.exists())
            self.assertTrue(result.html_report_path.exists())

            header = result.ranking_tsv_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
            for column in (
                "rank",
                "assay_id",
                "region_id",
                "classification",
                "score_status",
                "final_score",
                "inclusivity",
                "specificity",
                "conservation",
                "primer3_quality",
                "robustness",
                "original_compatible_count",
                "evaluation_sequence_count",
                "inclusivity_fraction",
                "compatible_off_target_hit_count",
                "plausible_off_target_count",
                "detectable_off_target_count",
                "pair_penalty",
                "reason_codes",
            ):
                with self.subTest(column=column):
                    self.assertIn(column, header)

            report = json.loads(result.ranking_report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["schema_version"], 1)
            self.assertEqual(report["status"], "COMPLETE")
            self.assertEqual(report["counts"]["assays"], 1)
            self.assertEqual(report["counts"]["in_silico_pass"], 1)
            self.assertEqual(report["counts"]["review"], 0)
            self.assertEqual(report["counts"]["high_risk"], 0)
            self.assertEqual(report["counts"]["complete_score"], 1)
            self.assertEqual(report["counts"]["incomplete_score"], 0)
            self.assertEqual(report["artifacts"]["ranking_tsv"], "ranking/assay_ranking.tsv")
            self.assertEqual(report["artifacts"]["html_report"], "report.html")
            self.assertEqual(report["assays"][0]["assay_id"], "a1")
            self.assertIn("components", report["assays"][0])
            self.assertIn("reasons", report["assays"][0])
            self.assertIn("weights", report["config"])

    def test_zero_assays_publish_header_only_tsv_zero_counts_and_empty_html(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir)
            primer = make_primer_result(assays=(), candidates=())
            inclusivity = make_inclusivity_result(primer)
            specificity = make_specificity_result(primer)
            result = evaluate_ranking(
                primer,
                inclusivity,
                specificity,
                RankingConfig(enabled=True),
                output,
                target_name="target",
            )
            self.assertEqual(len(result.ranking_tsv_path.read_text(encoding="utf-8").splitlines()), 1)
            report = json.loads(result.ranking_report_path.read_text(encoding="utf-8"))
            self.assertEqual(report["counts"]["assays"], 0)
            self.assertIn("No assay candidates", result.html_report_path.read_text(encoding="utf-8"))

    def test_json_and_tsv_are_deterministic_across_output_directories(self):
        primer = make_primer_result()
        inclusivity = make_inclusivity_result(primer)
        specificity = make_specificity_result(primer)
        outputs = []
        with tempfile.TemporaryDirectory() as first, tempfile.TemporaryDirectory() as second:
            for directory in (first, second):
                result = evaluate_ranking(
                    primer,
                    inclusivity,
                    specificity,
                    RankingConfig(enabled=True),
                    Path(directory),
                    target_name="target",
                )
                outputs.append(
                    (
                        result.ranking_tsv_path.read_bytes(),
                        result.ranking_report_path.read_bytes(),
                    )
                )
        self.assertEqual(outputs[0], outputs[1])

    def test_stage_scoring_uses_retention_total_even_when_retained_hits_are_truncated(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            primer = make_primer_result()
            inclusivity = make_inclusivity_result(primer)
            specificity = make_specificity_result(
                primer, hit_totals={("off", "a1", "FORWARD"): 25}
            )
            specificity = replace(
                specificity,
                retention=tuple(
                    replace(row, retained_hit_count=1, truncated=True)
                    if row.role == "FORWARD"
                    else row
                    for row in specificity.retention
                ),
            )
            result = evaluate_ranking(
                primer,
                inclusivity,
                specificity,
                RankingConfig(enabled=True),
                Path(tmpdir),
                target_name="target",
            )
            self.assertAlmostEqual(result.assays[0].components.specificity, 0.80)


if __name__ == "__main__":
    unittest.main()
