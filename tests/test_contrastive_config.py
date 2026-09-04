import tempfile
import unittest
from pathlib import Path

import qpcr_pipeline.config as config_module
from qpcr_pipeline.config import load_config


class ContrastiveConservationConfigTests(unittest.TestCase):
    def _load_yaml(self, text: str):
        with tempfile.TemporaryDirectory() as tmpdir:
            config_path = Path(tmpdir) / "config.yaml"
            config_path.write_text(text, encoding="utf-8")
            return load_config(config_path)

    @staticmethod
    def _base_yaml() -> str:
        return (
            "target:\n"
            "  name: target\n"
            "input:\n"
            "  fasta: tests/fixtures/target_small.fasta\n"
            "alignment:\n"
            "  enabled: true\n"
            "conservation:\n"
            "  enabled: true\n"
        )

    def test_config_type_exists_and_defaults_disabled(self):
        config_type = getattr(config_module, "ContrastiveConservationConfig", None)
        self.assertIsNotNone(config_type)

        config = self._load_yaml(self._base_yaml())
        self.assertEqual(config.contrastive_conservation, config_type())
        self.assertFalse(config.contrastive_conservation.enabled)

    def test_loads_enabled_contrastive_configuration(self):
        config = self._load_yaml(
            self._base_yaml()
            + "panel:\n"
            "  frozen_manifest: approved.json\n"
            "off_targets:\n"
            "  - name: challenge-a\n"
            "    fasta: challenge-a.fasta\n"
            "contrastive_conservation:\n"
            "  enabled: true\n"
        )

        self.assertTrue(config.contrastive_conservation.enabled)

    def test_rejects_non_boolean_enabled(self):
        with self.assertRaisesRegex(ValueError, "Contrastive conservation enabled must be a boolean"):
            self._load_yaml(
                self._base_yaml()
                + "contrastive_conservation:\n"
                "  enabled: 1\n"
            )

    def test_rejects_enabled_without_conservation(self):
        yaml = self._base_yaml().replace("conservation:\n  enabled: true\n", "")
        yaml += (
            "panel:\n"
            "  frozen_manifest: approved.json\n"
            "off_targets:\n"
            "  - name: challenge-a\n"
            "    fasta: challenge-a.fasta\n"
            "contrastive_conservation:\n"
            "  enabled: true\n"
        )
        with self.assertRaisesRegex(
            ValueError,
            "Enabled contrastive conservation requires enabled conservation",
        ):
            self._load_yaml(yaml)

    def test_rejects_enabled_without_frozen_panel(self):
        with self.assertRaisesRegex(
            ValueError,
            "Enabled contrastive conservation requires an approved frozen panel",
        ):
            self._load_yaml(
                self._base_yaml()
                + "off_targets:\n"
                "  - name: challenge-a\n"
                "    fasta: challenge-a.fasta\n"
                "contrastive_conservation:\n"
                "  enabled: true\n"
            )

    def test_rejects_enabled_without_off_targets(self):
        with self.assertRaisesRegex(
            ValueError,
            "Enabled contrastive conservation requires off-target datasets",
        ):
            self._load_yaml(
                self._base_yaml()
                + "panel:\n"
                "  frozen_manifest: approved.json\n"
                "contrastive_conservation:\n"
                "  enabled: true\n"
            )

    def test_rejects_unknown_contrastive_fields(self):
        with self.assertRaisesRegex(
            ValueError,
            "contrastive_conservation.*extra.*unrecognized",
        ):
            self._load_yaml(
                self._base_yaml()
                + "contrastive_conservation:\n"
                "  extra: true\n"
            )


if __name__ == "__main__":
    unittest.main()
