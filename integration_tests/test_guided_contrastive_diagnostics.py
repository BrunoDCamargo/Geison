from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import yaml


REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "examples" / "guided_demo" / "generate_demo_data.py"


def _run(command: list[str], *, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, cwd=cwd, capture_output=True, text=True, check=False)


def _generate_demo() -> tuple[tempfile.TemporaryDirectory[str], Path, Path]:
    temporary = tempfile.TemporaryDirectory()
    demo = Path(temporary.name) / "demo"
    generated = _run([sys.executable, str(GENERATOR), str(demo)], cwd=REPO_ROOT)
    if generated.returncode != 0:
        temporary.cleanup()
        raise AssertionError(generated.stderr or generated.stdout)
    return temporary, demo, demo / "output"


def _proposal_and_approve(demo: Path, output: Path) -> Path:
    proposal = _run(
        ["qpcr-pipeline", "run", str(demo / "config-proposal.yaml"), "--outdir", str(output)],
        cwd=demo,
    )
    if proposal.returncode != 3 or "PANEL_APPROVAL_REQUIRED" not in proposal.stdout:
        raise AssertionError(
            f"proposal rc={proposal.returncode}\nstdout:\n{proposal.stdout}\nstderr:\n{proposal.stderr}"
        )
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
    if approval.returncode != 0:
        raise AssertionError(
            f"approval rc={approval.returncode}\nstdout:\n{approval.stdout}\nstderr:\n{approval.stderr}"
        )
    return approved_panel


def _resume(demo: Path, output: Path, approved_panel: Path) -> subprocess.CompletedProcess[str]:
    template = yaml.safe_load((demo / "config-approved-template.yaml").read_text(encoding="utf-8"))
    template["panel"]["frozen_manifest"] = str(approved_panel)
    approved_config = demo / "config-approved.yaml"
    approved_config.write_text(yaml.safe_dump(template, sort_keys=False), encoding="utf-8")
    return _run(
        ["qpcr-pipeline", "run", str(approved_config), "--outdir", str(output), "--resume"],
        cwd=demo,
    )


def _completed_workflow() -> tuple[tempfile.TemporaryDirectory[str], Path, Path]:
    temporary, demo, output = _generate_demo()
    approved = _proposal_and_approve(demo, output)
    resumed = _resume(demo, output, approved)
    if resumed.returncode != 0:
        temporary.cleanup()
        raise AssertionError(
            f"resume rc={resumed.returncode}\nstdout:\n{resumed.stdout}\nstderr:\n{resumed.stderr}"
        )
    return temporary, demo, output


class GuidedContrastiveDiagnostics(unittest.TestCase):
    def test_gate_and_approval(self):
        temporary, demo, output = _generate_demo()
        try:
            approved = _proposal_and_approve(demo, output)
            self.assertTrue(approved.is_file())
            self.assertFalse((output / ".checkpoints" / "input" / "manifest.json").exists())
        finally:
            temporary.cleanup()

    def test_resume_reaches_completed(self):
        temporary, _, output = _completed_workflow()
        try:
            manifest = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(manifest["status"], "COMPLETED")
        finally:
            temporary.cleanup()

    def test_contrast_report_is_published(self):
        temporary, _, output = _completed_workflow()
        try:
            report = output / "contrastive_conservation" / "report.html"
            self.assertTrue(report.is_file())
            self.assertIn("Target vs non-target contrast", report.read_text(encoding="utf-8"))
        finally:
            temporary.cleanup()

    def test_candidate_retains_multiple_contributing_windows(self):
        temporary, _, output = _completed_workflow()
        try:
            candidate_path = output / "contrastive_conservation" / "candidate_regions.tsv"
            with candidate_path.open(encoding="utf-8", newline="") as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            self.assertTrue(candidates)
            self.assertTrue(
                any(len(json.loads(row["contributing_windows"])) >= 2 for row in candidates)
            )
        finally:
            temporary.cleanup()


if __name__ == "__main__":
    unittest.main()
