from __future__ import annotations

import csv
import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml

from qpcr_pipeline.evidence_bundle import create_evidence_bundle


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "examples" / "guided_demo" / "generate_demo_data.py"
TOOLS_AVAILABLE = shutil.which("mafft") is not None and shutil.which("primer3_core") is not None


def _run(command: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        cwd=cwd,
        capture_output=True,
        text=True,
        check=False,
    )


@unittest.skipUnless(TOOLS_AVAILABLE, "MAFFT and Primer3 are required for guided integration")
class GuidedContrastiveDemoIntegrationTest(unittest.TestCase):
    def test_installed_cli_completes_synthetic_guided_workflow(self):
        with tempfile.TemporaryDirectory() as temporary:
            demo = Path(temporary) / "demo"
            generated = _run([sys.executable, str(GENERATOR), str(demo)], cwd=REPO_ROOT)
            self.assertEqual(generated.returncode, 0, generated.stderr or generated.stdout)

            output = demo / "output"
            proposal = _run(
                [
                    "qpcr-pipeline",
                    "run",
                    str(demo / "config-proposal.yaml"),
                    "--outdir",
                    str(output),
                ],
                cwd=demo,
            )
            self.assertEqual(proposal.returncode, 3, proposal.stderr or proposal.stdout)
            self.assertIn("PANEL_APPROVAL_REQUIRED", proposal.stdout)
            self.assertFalse((output / ".checkpoints" / "input" / "manifest.json").exists())

            approved_panel = demo / "approved_panel.json"
            approval = _run(
                [
                    "qpcr-pipeline",
                    "panel",
                    "approve",
                    str(output / "panel_proposal.yaml"),
                    "--output",
                    str(approved_panel),
                ],
                cwd=demo,
            )
            self.assertEqual(approval.returncode, 0, approval.stderr or approval.stdout)

            template = yaml.safe_load(
                (demo / "config-approved-template.yaml").read_text(encoding="utf-8")
            )
            template["panel"]["frozen_manifest"] = str(approved_panel)
            approved_config = demo / "config-approved.yaml"
            approved_config.write_text(
                yaml.safe_dump(template, sort_keys=False),
                encoding="utf-8",
            )

            resumed = _run(
                [
                    "qpcr-pipeline",
                    "run",
                    str(approved_config),
                    "--outdir",
                    str(output),
                    "--resume",
                ],
                cwd=demo,
            )
            self.assertEqual(
                resumed.returncode,
                0,
                "stdout:\n" + resumed.stdout + "\nstderr:\n" + resumed.stderr,
            )

            contrast_dir = output / "contrastive_conservation"
            for name in (
                "window_metrics.tsv",
                "dataset_metrics.tsv",
                "candidate_regions.tsv",
                "contrastive_conservation_report.json",
                "report.html",
            ):
                self.assertTrue((contrast_dir / name).is_file(), name)

            with (contrast_dir / "candidate_regions.tsv").open(
                encoding="utf-8", newline=""
            ) as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(candidates)
            self.assertTrue(
                any(len(json.loads(row["contributing_windows"])) >= 2 for row in candidates),
                "expected at least one consolidated region to retain multiple raw windows",
            )

            primer_report = json.loads(
                (output / "primer_design" / "primer_design_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(
                primer_report["candidate_source"], "CONTRASTIVE_CONSERVATION"
            )
            self.assertTrue(
                primer_report["assays"],
                "expected at least one Primer3 assay constrained by contrastive evidence",
            )
            candidates_by_id = {
                candidate["region_id"]: candidate
                for candidate in primer_report["candidates"]
            }
            for assay in primer_report["assays"]:
                candidate = candidates_by_id[assay["region_id"]]
                self.assertLessEqual(
                    assay["forward_primer"]["reference_start"],
                    candidate["peak_start"],
                    assay["assay_id"],
                )
                self.assertGreaterEqual(
                    assay["reverse_primer"]["reference_end"],
                    candidate["peak_end"],
                    assay["assay_id"],
                )

            inclusivity_report = json.loads(
                (output / "inclusivity" / "inclusivity_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(inclusivity_report["status"], "COMPLETE")
            self.assertEqual(
                inclusivity_report["counts"]["evaluation_sequences"],
                4,
            )

            specificity_report = json.loads(
                (output / "specificity" / "specificity_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(specificity_report["status"], "COMPLETE")

            ranking_report = json.loads(
                (output / "ranking" / "ranking_report.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(ranking_report["status"], "COMPLETE")
            self.assertTrue(ranking_report["assays"])
            classifications = {
                assay["classification"] for assay in ranking_report["assays"]
            }
            self.assertTrue(
                classifications & {"IN SILICO PASS", "REVIEW"},
                "the corrected synthetic demo must not end with every assay HIGH_RISK",
            )
            self.assertTrue(
                any(
                    all(
                        reason.get("code") != "DETECTABLE_OFF_TARGET"
                        for reason in assay.get("reasons", [])
                    )
                    for assay in ranking_report["assays"]
                ),
                "expected at least one anchored assay without detectable off-target amplification",
            )

            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "COMPLETED")

            report_path = output / "report.html"
            self.assertTrue(report_path.is_file())
            report_html = report_path.read_text(encoding="utf-8")
            self.assertIn("Geison Researcher Report", report_html)
            if "IN SILICO PASS" in classifications:
                self.assertIn("In-silico candidate(s) identified", report_html)
            else:
                self.assertIn(
                    "No in-silico pass; candidate(s) require review",
                    report_html,
                )

            bundle = create_evidence_bundle(
                output,
                demo / "evidence_bundle.zip",
                extra_files=(approved_config, approved_panel),
            )
            self.assertTrue(bundle.is_file())
            self.assertGreater(bundle.stat().st_size, 0)


if __name__ == "__main__":
    unittest.main()
