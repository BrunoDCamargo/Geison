import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class PipelineCliTests(unittest.TestCase):
    def test_run_command_loads_minimal_configuration(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(
            executable,
            "qpcr-pipeline console command is not installed",
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [executable, "run", str(config_path)],
                capture_output=True,
                text=True,
                check=False,
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("synthetic-target", result.stdout)


if __name__ == "__main__":
    unittest.main()
