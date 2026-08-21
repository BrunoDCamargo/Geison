import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import PipelineConfig, load_config


class PipelineConfigTests(unittest.TestCase):
    def test_loads_minimal_yaml_configuration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n",
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIsInstance(config, PipelineConfig)
        self.assertEqual(config.target_name, "synthetic-target")
        self.assertEqual(config.input_fasta, Path("tests/fixtures/target_small.fasta"))

    def test_loads_genbank_input_and_optional_qc_thresholds(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  genbank: tests/fixtures/target.gb\n"
                "qc:\n"
                "  min_length: 100\n"
                "  max_ambiguous_fraction: 0.05\n"
                "  expected_length: 150\n"
                "  length_tolerance_fraction: 0.10\n",
                encoding="utf-8",
            )

            config = load_config(config_path)

        self.assertIsNone(config.input_fasta)
        self.assertEqual(config.input_genbank, Path("tests/fixtures/target.gb"))
        self.assertEqual(config.selected_input, (Path("tests/fixtures/target.gb"), "genbank"))
        self.assertEqual(config.qc.min_length, 100)
        self.assertEqual(config.qc.max_ambiguous_fraction, 0.05)
        self.assertEqual(config.qc.expected_length, 150)
        self.assertEqual(config.qc.length_tolerance_fraction, 0.10)

    def test_rejects_fractional_min_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "qc:\n"
                "  min_length: 100.5\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "qc.min_length.*integer"):
                load_config(config_path)

    def test_rejects_fractional_expected_length(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "qc:\n"
                "  expected_length: 150.5\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "qc.expected_length.*integer"):
                load_config(config_path)

    def test_requires_exactly_one_local_sequence_input(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(
                "target:\n"
                "  name: synthetic-target\n"
                "input:\n"
                "  fasta: tests/fixtures/target_small.fasta\n"
                "  genbank: tests/fixtures/target.gb\n",
                encoding="utf-8",
            )

            with self.assertRaisesRegex(ValueError, "Exactly one local sequence input"):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
