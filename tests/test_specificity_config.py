import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import (
    AlignmentConfig,
    ConservationConfig,
    OffTargetConfig,
    PipelineConfig,
    PrimerDesignConfig,
    SpecificityConfig,
    load_config,
    validate_off_target_config,
    validate_specificity_config,
)


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class SpecificityConfigTests(unittest.TestCase):
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
        )

    def test_defaults_are_disabled_and_have_no_off_targets(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.off_targets, ())
        self.assertEqual(config.specificity, SpecificityConfig())

    def test_loads_fasta_and_frozen_off_targets_in_configured_order(self):
        config = self._load_yaml(
            self._enabled_base()
            + "off_targets:\n"
            "  - name: z-first\n    fasta: data/human.fasta\n"
            "  - name: a-second\n    frozen_dataset: data/neighbors\n"
            "specificity:\n"
            "  enabled: true\n"
            "  max_hits_per_oligo_per_dataset: 7\n"
            "  max_primer_mismatches: 1\n"
            "  max_probe_mismatches: 0\n"
            "  reject_primer_3_prime_mismatch: false\n"
            "  primer_3_prime_bases: 4\n"
            "  max_amplicon_size: 600\n"
        )
        self.assertEqual(
            config.off_targets,
            (
                OffTargetConfig(name="z-first", fasta=Path("data/human.fasta")),
                OffTargetConfig(name="a-second", frozen_dataset=Path("data/neighbors")),
            ),
        )
        self.assertEqual(
            config.specificity,
            SpecificityConfig(
                enabled=True,
                max_hits_per_oligo_per_dataset=7,
                max_primer_mismatches=1,
                max_probe_mismatches=0,
                reject_primer_3_prime_mismatch=False,
                primer_3_prime_bases=4,
                max_amplicon_size=600,
            ),
        )

    def test_rejects_invalid_off_target_sources_and_duplicate_names(self):
        invalid = (
            (
                "off_targets:\n  - name: x\n",
                "exactly one",
            ),
            (
                "off_targets:\n  - name: x\n    fasta: a.fa\n    frozen_dataset: frozen\n",
                "exactly one",
            ),
            (
                "off_targets:\n  - name: x\n    fasta: a.fa\n  - name: x\n    fasta: b.fa\n",
                "unique",
            ),
            (
                "off_targets:\n  - name: '   '\n    fasta: a.fa\n",
                "name",
            ),
        )
        for suffix, message in invalid:
            with self.subTest(suffix=suffix):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_yaml(self._enabled_base() + suffix)

    def test_rejects_invalid_specificity_values_and_unknown_fields(self):
        invalid = (
            ("specificity:\n  enabled: nope\n", "enabled"),
            ("specificity:\n  max_primer_mismatches: -1\n", "max_primer_mismatches"),
            ("specificity:\n  max_probe_mismatches: true\n", "max_probe_mismatches"),
            ("specificity:\n  max_hits_per_oligo_per_dataset: 0\n", "max_hits_per_oligo_per_dataset"),
            ("specificity:\n  primer_3_prime_bases: 0\n", "primer_3_prime_bases"),
            ("specificity:\n  max_amplicon_size: 0\n", "max_amplicon_size"),
            ("specificity:\n  surprise: 1\n", "unrecognized"),
        )
        for suffix, message in invalid:
            with self.subTest(suffix=suffix):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_yaml(self._enabled_base() + suffix)

    def test_enabled_specificity_requires_primer_design_and_off_targets(self):
        with self.assertRaisesRegex(ValueError, "specificity.*requires enabled primer design"):
            PipelineConfig(
                target_name="target",
                input_fasta=FIXTURE_FASTA,
                off_targets=(OffTargetConfig(name="x", fasta=Path("x.fa")),),
                specificity=SpecificityConfig(enabled=True),
            ).selected_input

        with self.assertRaisesRegex(ValueError, "specificity.*at least one off-target"):
            PipelineConfig(
                target_name="target",
                input_fasta=FIXTURE_FASTA,
                alignment=AlignmentConfig(enabled=True),
                conservation=ConservationConfig(enabled=True),
                primer_design=PrimerDesignConfig(enabled=True),
                specificity=SpecificityConfig(enabled=True),
            ).selected_input

    def test_direct_validators_reject_wrong_types(self):
        with self.assertRaisesRegex(ValueError, "Off-target configuration"):
            validate_off_target_config(object())
        with self.assertRaisesRegex(ValueError, "Specificity configuration"):
            validate_specificity_config(object())


if __name__ == "__main__":
    unittest.main()
