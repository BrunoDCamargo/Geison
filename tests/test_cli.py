import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path


class PipelineCliTests(unittest.TestCase):
    def _executable(self):
        executable = shutil.which("qpcr-pipeline")
        self.assertIsNotNone(
            executable,
            "qpcr-pipeline console command is not installed",
        )
        return executable

    @staticmethod
    def _config(tmpdir):
        config_path = Path(tmpdir) / "config.yaml"
        config_path.write_text(
            "target:\n"
            "  name: synthetic-target\n"
            "input:\n"
            "  fasta: tests/fixtures/target_small.fasta\n",
            encoding="utf-8",
        )
        return config_path

    def _run(self, tmpdir, *args):
        return subprocess.run(
            [self._executable(), "run", str(self._config(tmpdir)), *args],
            capture_output=True,
            text=True,
            check=False,
        )

    def test_run_command_loads_minimal_configuration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(tmpdir)

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("synthetic-target", result.stdout)

    def test_resume_runs_with_outdir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(tmpdir, "--outdir", str(outdir), "--resume")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("COMPLETED", result.stdout)

    def test_resume_force_step_is_accepted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(
                tmpdir,
                "--outdir",
                str(outdir),
                "--resume",
                "--force-step",
                "specificity",
            )

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("COMPLETED", result.stdout)

    def test_resume_and_from_step_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(
                tmpdir,
                "--outdir",
                str(outdir),
                "--resume",
                "--from-step",
                "specificity",
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("resume", result.stderr.lower())
        self.assertIn("from-step", result.stderr.lower())

    def test_from_step_and_force_step_are_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(
                tmpdir,
                "--outdir",
                str(outdir),
                "--from-step",
                "specificity",
                "--force-step",
                "specificity",
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("from-step", result.stderr.lower())
        self.assertIn("force-step", result.stderr.lower())

    def test_force_step_without_resume_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(
                tmpdir,
                "--outdir",
                str(outdir),
                "--force-step",
                "specificity",
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("force-step", result.stderr.lower())
        self.assertIn("resume", result.stderr.lower())

    def test_unknown_stage_is_rejected_by_argparse(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "run"
            result = self._run(
                tmpdir,
                "--outdir",
                str(outdir),
                "--from-step",
                "unknown",
            )

        self.assertNotEqual(result.returncode, 0)
        self.assertIn("invalid choice", result.stderr.lower())

    def test_resume_controls_require_outdir(self):
        for flag in ("--resume", "--from-step", "--force-step"):
            with self.subTest(flag=flag), tempfile.TemporaryDirectory() as tmpdir:
                args = [flag]
                if flag != "--resume":
                    args.append("specificity")
                result = self._run(tmpdir, *args)

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("outdir", result.stderr.lower())


if __name__ == "__main__":
    unittest.main()
