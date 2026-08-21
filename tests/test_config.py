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


if __name__ == "__main__":
    unittest.main()
