import tempfile
from pathlib import Path

import pytest

from qpcr_pipeline.config import ContrastiveConservationConfig, load_config
from qpcr_pipeline.execution import STAGE_DEPENDENCIES, STAGE_ORDER


def _load_yaml(text: str):
    with tempfile.TemporaryDirectory() as tmpdir:
        path = Path(tmpdir) / "config.yaml"
        path.write_text(text, encoding="utf-8")
        return load_config(path)


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


def test_contrastive_config_defaults_to_disabled():
    config = _load_yaml(_base_yaml())
    assert config.contrastive_conservation == ContrastiveConservationConfig()


def test_contrastive_config_parses_enabled_with_frozen_panel_and_off_target():
    config = _load_yaml(
        _base_yaml()
        + "panel:\n"
        "  frozen_manifest: approved.json\n"
        "off_targets:\n"
        "  - name: challenge\n"
        "    fasta: tests/fixtures/target_small.fasta\n"
        "contrastive_conservation:\n"
        "  enabled: true\n"
    )
    assert config.contrastive_conservation.enabled is True


def test_enabled_contrastive_requires_conservation():
    yaml = _base_yaml().replace("conservation:\n  enabled: true\n", "")
    with pytest.raises(ValueError, match="contrastive conservation requires enabled conservation"):
        _load_yaml(
            yaml
            + "panel:\n"
            "  frozen_manifest: approved.json\n"
            "off_targets:\n"
            "  - name: challenge\n"
            "    fasta: tests/fixtures/target_small.fasta\n"
            "contrastive_conservation:\n"
            "  enabled: true\n"
        )


def test_enabled_contrastive_requires_frozen_panel():
    with pytest.raises(ValueError, match="requires an approved frozen panel"):
        _load_yaml(
            _base_yaml()
            + "off_targets:\n"
            "  - name: challenge\n"
            "    fasta: tests/fixtures/target_small.fasta\n"
            "contrastive_conservation:\n"
            "  enabled: true\n"
        )


def test_enabled_contrastive_requires_off_targets():
    with pytest.raises(ValueError, match="requires off-target datasets"):
        _load_yaml(
            _base_yaml()
            + "panel:\n"
            "  frozen_manifest: approved.json\n"
            "contrastive_conservation:\n"
            "  enabled: true\n"
        )


def test_contrastive_enabled_must_be_boolean():
    with pytest.raises(ValueError, match="contrastive conservation enabled must be a boolean"):
        _load_yaml(_base_yaml() + "contrastive_conservation:\n  enabled: 1\n")


def test_stage_graph_places_contrast_between_conservation_and_primer_design():
    assert STAGE_ORDER == (
        "panel",
        "input",
        "qc",
        "clustering",
        "alignment",
        "conservation",
        "contrastive_conservation",
        "primer_design",
        "inclusivity",
        "specificity",
        "ranking",
    )
    assert STAGE_DEPENDENCIES["panel"] == ()
    assert STAGE_DEPENDENCIES["input"] == ()
    assert STAGE_DEPENDENCIES["contrastive_conservation"] == ("panel", "conservation")
    assert STAGE_DEPENDENCIES["primer_design"] == ("contrastive_conservation",)
