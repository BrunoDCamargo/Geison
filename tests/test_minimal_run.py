import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURE_FASTA = REPO_ROOT / "tests" / "fixtures" / "target_small.fasta"


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


if __name__ == "__main__":
    unittest.main()
