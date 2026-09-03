import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from qpcr_pipeline.alignment import SubprocessMafftRunner
from qpcr_pipeline.config import AlignmentConfig


class MafftDirectionModeTests(unittest.TestCase):
    def test_runner_uses_fast_direction_adjustment_not_accurate_mode(self):
        completed = subprocess.CompletedProcess(
            args=[], returncode=0, stdout=">aligned\nACGT\n", stderr=""
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            input_path = Path(tmpdir) / "input.fasta"
            output_path = Path(tmpdir) / "output.fasta"

            with patch(
                "qpcr_pipeline.alignment.shutil.which", return_value="/usr/bin/mafft"
            ), patch(
                "qpcr_pipeline.alignment.subprocess.run", return_value=completed
            ) as run:
                SubprocessMafftRunner().run(
                    input_path,
                    output_path,
                    AlignmentConfig(enabled=True, threads=2),
                )

        command = run.call_args.args[0]
        self.assertIn("--adjustdirection", command)
        self.assertNotIn("--adjustdirectionaccurately", command)


if __name__ == "__main__":
    unittest.main()
