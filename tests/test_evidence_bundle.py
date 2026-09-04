from __future__ import annotations

import json
from pathlib import Path
import tempfile
import unittest
import zipfile

from qpcr_pipeline.evidence_bundle import create_evidence_bundle


class EvidenceBundleTests(unittest.TestCase):
    def _write(self, root: Path, relative: str, text: str = "x") -> Path:
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")
        return path

    def test_bundle_contains_published_evidence_and_excludes_unrelated_files(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            self._write(output, "report.html", "<html>report</html>")
            self._write(output, "run_manifest.json", json.dumps({"status": "COMPLETED"}))
            self._write(output, "run_summary.json", json.dumps({"status": "COMPLETED"}))
            self._write(output, "qc_report.json", json.dumps({"records": []}))
            self._write(output, "panel/approved_panel.json", json.dumps({"status": "APPROVED"}))
            self._write(output, "contrastive_conservation/window_metrics.tsv")
            self._write(output, "contrastive_conservation/contrastive_conservation_report.json")
            self._write(output, "primer_design/primer_design_report.json")
            self._write(output, "inclusivity/inclusivity_report.json")
            self._write(output, "specificity/specificity_report.json")
            self._write(output, "ranking/assay_ranking.tsv")
            self._write(output, "ranking/ranking_report.json")
            self._write(output, ".checkpoints/primer_design/manifest.json")
            self._write(output, ".checkpoints/primer_design/result.json", "must-not-be-bundled")
            self._write(output, "unrelated-secret.txt", "do-not-package")
            approved_config = self._write(root, "config-approved.yaml", "target: synthetic")

            destination = output / "evidence_bundle.zip"
            bundle = create_evidence_bundle(
                output,
                destination,
                extra_files=(approved_config,),
            )

            with zipfile.ZipFile(bundle) as archive:
                names = archive.namelist()
                content = {name: archive.read(name) for name in names}

        self.assertEqual(names, sorted(names))
        self.assertIn("report.html", names)
        self.assertIn("run_manifest.json", names)
        self.assertIn("run_summary.json", names)
        self.assertIn("qc_report.json", names)
        self.assertIn("panel/approved_panel.json", names)
        self.assertIn("contrastive_conservation/window_metrics.tsv", names)
        self.assertIn("primer_design/primer_design_report.json", names)
        self.assertIn("inclusivity/inclusivity_report.json", names)
        self.assertIn("specificity/specificity_report.json", names)
        self.assertIn("ranking/ranking_report.json", names)
        self.assertIn(".checkpoints/primer_design/manifest.json", names)
        self.assertIn("inputs/config-approved.yaml", names)
        self.assertNotIn(".checkpoints/primer_design/result.json", names)
        self.assertNotIn("unrelated-secret.txt", names)
        self.assertNotIn("evidence_bundle.zip", names)
        self.assertEqual(content["inputs/config-approved.yaml"], b"target: synthetic")

    def test_missing_optional_artifacts_do_not_prevent_bundle_creation(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            self._write(output, "run_manifest.json", json.dumps({"status": "PARTIAL"}))

            bundle = create_evidence_bundle(
                output,
                root / "partial-evidence.zip",
            )

            with zipfile.ZipFile(bundle) as archive:
                names = archive.namelist()

        self.assertEqual(names, ["run_manifest.json"])

    def test_extra_file_must_be_a_regular_file(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "output"
            self._write(output, "run_manifest.json", "{}")
            missing = root / "missing-config.yaml"

            with self.assertRaisesRegex(ValueError, "extra file"):
                create_evidence_bundle(
                    output,
                    root / "bundle.zip",
                    extra_files=(missing,),
                )


if __name__ == "__main__":
    unittest.main()
