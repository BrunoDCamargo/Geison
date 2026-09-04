from pathlib import Path
from unittest.mock import patch

import pytest

from qpcr_pipeline.challenge_panel import resolve_challenge_datasets
from qpcr_pipeline.config import OffTargetConfig
from qpcr_pipeline.off_targets import OffTargetDataset
from qpcr_pipeline.panel import (
    DiagnosticContext,
    PanelDefinition,
    PanelNonTarget,
    PanelTarget,
    TargetGroup,
)
from qpcr_pipeline.panel_manifest import ApprovedPanelManifest


def _target() -> PanelTarget:
    return PanelTarget(
        name="target",
        taxid=None,
        mode="broad_detection",
        subtype=None,
        groups=(
            TargetGroup(
                name="target-group",
                required=True,
                dataset_roles=("DESIGN",),
                reasons=("fixture",),
                proposed_by=("test",),
            ),
        ),
    )


def _non_target(
    name: str,
    criticality: str = "IMPORTANT",
    roles=("CHALLENGE",),
) -> PanelNonTarget:
    return PanelNonTarget(
        name=name,
        taxid=None,
        criticality=criticality,
        dataset_roles=roles,
        reasons=("fixture",),
        proposed_by=("test",),
    )


def _manifest(*items: PanelNonTarget) -> ApprovedPanelManifest:
    return ApprovedPanelManifest(
        schema_version=1,
        status="APPROVED",
        approved_by_user=True,
        proposal_sha256="sha256:" + "a" * 64,
        definition=PanelDefinition(
            target=_target(),
            non_targets=tuple(items),
            diagnostic_context=DiagnosticContext(),
        ),
    )


def _dataset(name: str) -> OffTargetDataset:
    return OffTargetDataset(
        name=name,
        source_type="FASTA",
        source_path=Path(f"{name}.fasta"),
        sha256="b" * 64,
        sequence_ids=(),
        records=(),
        frozen_manifest_path=None,
        frozen_manifest=None,
    )


def test_resolver_uses_only_challenge_entries_in_panel_order_and_criticality():
    manifest = _manifest(
        _non_target("Challenge A", "CRITICAL"),
        _non_target("Design only", "BACKGROUND", ("DESIGN",)),
        _non_target("Challenge B", "IMPORTANT"),
    )
    configs = (
        OffTargetConfig(name="challenge b", fasta=Path("b.fasta")),
        OffTargetConfig(name=" challenge a ", fasta=Path("a.fasta")),
        OffTargetConfig(name="Design only", fasta=Path("design.fasta")),
    )

    def fake_loader(config: OffTargetConfig) -> OffTargetDataset:
        return _dataset(config.name)

    with patch(
        "qpcr_pipeline.challenge_panel.load_off_target_dataset",
        side_effect=fake_loader,
    ) as loader:
        resolved = resolve_challenge_datasets(manifest, configs)

    assert [item.name for item in resolved] == ["Challenge A", "Challenge B"]
    assert [item.criticality for item in resolved] == ["CRITICAL", "IMPORTANT"]
    assert [item.dataset.name for item in resolved] == [" challenge a ", "challenge b"]
    assert loader.call_count == 2


def test_missing_challenge_mapping_is_rejected():
    manifest = _manifest(_non_target("Missing", "CRITICAL"))
    with pytest.raises(
        ValueError,
        match="Challenge dataset missing for approved panel entry 'Missing'",
    ):
        resolve_challenge_datasets(manifest, ())


def test_duplicate_normalized_off_target_names_are_rejected():
    manifest = _manifest(_non_target("Challenge", "CRITICAL"))
    configs = (
        OffTargetConfig(name="Challenge", fasta=Path("a.fasta")),
        OffTargetConfig(name=" challenge ", fasta=Path("b.fasta")),
    )
    with pytest.raises(
        ValueError,
        match="unique after normalization",
    ):
        resolve_challenge_datasets(manifest, configs)


def test_approved_panel_with_zero_challenge_entries_is_rejected():
    manifest = _manifest(
        _non_target("Design only", "BACKGROUND", ("DESIGN",)),
    )
    configs = (OffTargetConfig(name="Design only", fasta=Path("design.fasta")),)
    with pytest.raises(
        ValueError,
        match="at least one approved CHALLENGE",
    ):
        resolve_challenge_datasets(manifest, configs)
