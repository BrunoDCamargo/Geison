import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from qpcr_pipeline.alignment import SubprocessMafftRunner, align_discovery
from qpcr_pipeline.config import AlignmentConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


class _FakeMafftRunner:
    def run(self, input_path, output_path, config):
        del input_path, config
        Path(output_path).write_text(
            ">geison-00000000\nACGT\n>geison-00000001\nACGA\n",
            encoding="utf-8",
        )


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

    def test_alignment_report_records_fast_direction_adjustment(self):
        records = (
            LocalSequenceRecord(sequence_id="first", sequence="ACGT"),
            LocalSequenceRecord(sequence_id="second", sequence="ACGA"),
        )
        discovery = DiscoverySet(("first", "second"))

        with tempfile.TemporaryDirectory() as tmpdir:
            result = align_discovery(
                records,
                discovery,
                AlignmentConfig(enabled=True, threads=2),
                Path(tmpdir),
                runner=_FakeMafftRunner(),
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

        parameters = report["tool"]["parameters"]
        self.assertTrue(parameters["adjust_direction"])
        self.assertNotIn("adjust_direction_accurately", parameters)


if __name__ == "__main__":
    unittest.main()
