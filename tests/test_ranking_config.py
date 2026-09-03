import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import (
    AlignmentConfig,
    ConservationConfig,
    PanelConfig,
    PipelineConfig,
    PrimerDesignConfig,
    RankingConfig,
    RankingWeights,
    load_config,
    validate_ranking_config,
)


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class RankingConfigTests(unittest.TestCase):
    def _load_yaml(self, text: str) -> PipelineConfig:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "config.yaml"
            path.write_text(text, encoding="utf-8")
            return load_config(path)

    def _enabled_base(self) -> str:
        return (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
            "primer_design:\n  enabled: true\n"
            "panel:\n  frozen_manifest: approved.json\n"
        )

    def test_defaults_are_disabled_with_stable_weights(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.ranking, RankingConfig())
        self.assertEqual(config.ranking.weights, RankingWeights())

    def test_loads_complete_ranking_configuration(self):
        config = self._load_yaml(
            self._enabled_base()
            + "ranking:\n"
            "  enabled: true\n"
            "  min_inclusivity_for_pass: 0.98\n"
            "  min_inclusivity_before_high_risk: 0.85\n"
            "  weights:\n"
            "    inclusivity: 0.40\n"
            "    specificity: 0.30\n"
            "    conservation: 0.15\n"
            "    primer3_quality: 0.10\n"
            "    robustness: 0.05\n"
        )
        self.assertEqual(
            config.ranking,
            RankingConfig(
                enabled=True,
                min_inclusivity_for_pass=0.98,
                min_inclusivity_before_high_risk=0.85,
                weights=RankingWeights(
                    inclusivity=0.40,
                    specificity=0.30,
                    conservation=0.15,
                    primer3_quality=0.10,
                    robustness=0.05,
                ),
            ),
        )

    def test_rejects_invalid_ranking_values_and_unknown_fields(self):
        invalid = (
            ("ranking:\n  enabled: nope\n", "enabled"),
            ("ranking:\n  min_inclusivity_for_pass: 1.1\n", "min_inclusivity_for_pass"),
            (
                "ranking:\n  min_inclusivity_before_high_risk: -0.1\n",
                "min_inclusivity_before_high_risk",
            ),
            (
                "ranking:\n"
                "  min_inclusivity_for_pass: 0.8\n"
                "  min_inclusivity_before_high_risk: 0.9\n",
                "min_inclusivity_before_high_risk",
            ),
            (
                "ranking:\n"
                "  weights:\n"
                "    inclusivity: 0.50\n"
                "    specificity: 0.50\n"
                "    conservation: 0.50\n"
                "    primer3_quality: 0.00\n"
                "    robustness: 0.00\n",
                "sum",
            ),
            ("ranking:\n  surprise: 1\n", "unrecognized"),
            ("ranking:\n  weights:\n    surprise: 1\n", "unrecognized"),
        )
        for suffix, message in invalid:
            with self.subTest(suffix=suffix):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_yaml(self._enabled_base() + suffix)

    def test_enabled_ranking_requires_enabled_primer_design(self):
        with self.assertRaisesRegex(
            ValueError, "ranking.*requires enabled primer design"
        ):
            PipelineConfig(
                target_name="target",
                input_fasta=FIXTURE_FASTA,
                ranking=RankingConfig(enabled=True),
            ).selected_input

    def test_direct_validator_rejects_wrong_type(self):
        with self.assertRaisesRegex(ValueError, "Ranking configuration"):
            validate_ranking_config(object())

    def test_direct_config_accepts_valid_enabled_dependency_chain(self):
        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            alignment=AlignmentConfig(enabled=True),
            conservation=ConservationConfig(enabled=True),
            primer_design=PrimerDesignConfig(enabled=True),
            panel=PanelConfig(frozen_manifest=Path("approved.json")),
            ranking=RankingConfig(enabled=True),
        )
        self.assertEqual(config.selected_input, (FIXTURE_FASTA, "fasta"))


if __name__ == "__main__":
    unittest.main()
