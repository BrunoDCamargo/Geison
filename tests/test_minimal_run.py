import json
import os
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_FASTA = REPO_ROOT / "tests" / "fixtures" / "target_small.fasta"
LOCAL_PACKAGE_ENV = {
    **os.environ,
    "PYTHONPATH": str(REPO_ROOT) + os.pathsep + os.environ.get("PYTHONPATH", ""),
}


class MinimalPipelineRunTests(unittest.TestCase):
    def test_run_creates_completed_summary_for_fixture(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(executable, "qpcr-pipeline console command is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            config_path = tmp_path / "config.yaml"
            outdir = tmp_path / "run"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                f"  fasta: {FIXTURE_FASTA.as_posix()}\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path), "--outdir", str(outdir)],
                capture_output=True,
                text=True,
                check=False,
                env=LOCAL_PACKAGE_ENV,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("COMPLETED", result.stdout)

            summary_path = outdir / "run_summary.json"
            self.assertTrue(summary_path.exists(), "run_summary.json was not created")
            summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["status"], "COMPLETED")
        self.assertEqual(summary["target_name"], "synthetic-target")
        self.assertEqual(summary["sequence_count"], 3)
        self.assertEqual(summary["sequence_ids"], ["seq-1", "seq-2", "seq-3"])

    def test_run_writes_traceable_qc_report_for_approved_sequences(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(executable, "qpcr-pipeline console command is not installed")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            fasta_path = tmp_path / "target.fasta"
            config_path = tmp_path / "config.yaml"
            outdir = tmp_path / "run"
            fasta_path.write_text(
                ">accepted-1\n"
                "ACGTACGT\n"
                ">invalid\n"
                "ACGTXCGT\n"
                ">accepted-2\n"
                "ACGTACGA\n",
                encoding="utf-8",
            )
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                f"  fasta: {fasta_path.as_posix()}\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path), "--outdir", str(outdir)],
                capture_output=True,
                text=True,
                check=False,
                env=LOCAL_PACKAGE_ENV,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            summary = json.loads((outdir / "run_summary.json").read_text(encoding="utf-8"))
            qc_report_path = outdir / "qc_report.json"
            self.assertTrue(qc_report_path.exists(), "qc_report.json was not created")
            qc_report = json.loads(qc_report_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["sequence_ids"], ["accepted-1", "accepted-2"])
        self.assertEqual(
            qc_report["records"],
            [
                {"sequence_id": "accepted-1", "status": "ACCEPTED", "reason_codes": []},
                {"sequence_id": "invalid", "status": "REJECTED", "reason_codes": ["INVALID_NUCLEOTIDE"]},
                {"sequence_id": "accepted-2", "status": "ACCEPTED", "reason_codes": []},
            ],
        )
        self.assertEqual(qc_report["target_sequence_set"]["sequence_ids"], ["accepted-1", "accepted-2"])
        self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["accepted-1", "accepted-2"])


if __name__ == "__main__":
    unittest.main()
