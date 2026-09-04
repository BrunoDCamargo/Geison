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
            self.assertTrue((output / "specificity" / "specificity_report.json").is_file())
            self.assertTrue((output / "ranking" / "ranking_report.json").is_file())

            manifest = json.loads(
                (output / "run_manifest.json").read_text(encoding="utf-8")
            )
            self.assertEqual(manifest["status"], "COMPLETED")


if __name__ == "__main__":
    unittest.main()
