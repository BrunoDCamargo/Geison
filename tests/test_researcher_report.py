from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest

from qpcr_pipeline.researcher_report import (
    ResearcherReportError,
    load_researcher_report_data,
    scientific_outcome,
)


class ScientificOutcomeTests(unittest.TestCase):
    def test_completed_with_pass_has_positive_outcome(self) -> None:
        title, body = scientific_outcome(
            "COMPLETED",
            {"assays": [{"classification": "IN SILICO PASS", "reasons": []}]},
        )

        self.assertEqual(title, "In-silico candidate(s) identified")
        self.assertIn("IN SILICO PASS", body)

    def test_completed_without_pass_but_with_review_requires_review(self) -> None:
        title, body = scientific_outcome(
            "COMPLETED",
            {"assays": [{"classification": "REVIEW", "reasons": []}]},
        )

        self.assertEqual(title, "No in-silico pass; candidate(s) require review")
        self.assertIn("REVIEW", body)

    def test_completed_all_high_risk_has_negative_outcome(self) -> None:
        title, body = scientific_outcome(
            "COMPLETED",
            {
                "assays": [
                    {
                        "classification": "HIGH_RISK",
                        "reasons": [{"code": "DETECTABLE_OFF_TARGET"}],
                    }
                ]
            },
        )

        self.assertEqual(title, "No in-silico acceptable assay candidates identified")
        self.assertIn("HIGH_RISK", body)
        self.assertIn("detectable off-target", body.lower())

    def test_partial_is_inconclusive(self) -> None:
        title, _ = scientific_outcome("PARTIAL", None)
        self.assertEqual(title, "Inconclusive - insufficient evidence")

    def test_failed_has_no_conclusive_scientific_outcome(self) -> None:
        title, body = scientific_outcome("FAILED", None)
        self.assertEqual(title, "Execution failed - no conclusive scientific outcome")
        self.assertIn("failed", body.lower())

    def test_completed_with_zero_assays_is_inconclusive(self) -> None:
        title, _ = scientific_outcome("COMPLETED", {"assays": []})
        self.assertEqual(title, "Inconclusive - insufficient evidence")


class ResearcherReportDataTests(unittest.TestCase):
    def _write_json(self, root: Path, relative: str, payload: object) -> None:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload), encoding="utf-8")

    def test_loads_known_artifacts_and_keeps_missing_sections_none(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write_json(root, "run_manifest.json", {"status": "PARTIAL"})
            self._write_json(root, "run_summary.json", {"sequence_count": 4})
            self._write_json(root, "panel/approved_panel.json", {"status": "APPROVED"})
            self._write_json(
                root,
                "ranking/ranking_report.json",
                {"status": "COMPLETE", "assays": []},
            )

            data = load_researcher_report_data(root)

        self.assertEqual(data.run_manifest["status"], "PARTIAL")
        self.assertEqual(data.run_summary["sequence_count"], 4)
        self.assertEqual(data.panel["status"], "APPROVED")
        self.assertEqual(data.ranking["status"], "COMPLETE")
        self.assertIsNone(data.conservation)
        self.assertIsNone(data.contrastive)
        self.assertIsNone(data.primer_design)
        self.assertIsNone(data.inclusivity)
        self.assertIsNone(data.specificity)

    def test_requires_run_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            with self.assertRaisesRegex(ResearcherReportError, "run_manifest.json"):
                load_researcher_report_data(Path(temporary))

    def test_malformed_known_json_is_not_silently_ignored(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            self._write_json(root, "run_manifest.json", {"status": "COMPLETED"})
            path = root / "specificity" / "specificity_report.json"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("{not-json", encoding="utf-8")

            with self.assertRaisesRegex(
                ResearcherReportError, "specificity/specificity_report.json"
            ):
                load_researcher_report_data(root)


if __name__ == "__main__":
    unittest.main()
