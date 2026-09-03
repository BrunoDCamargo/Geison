from pathlib import Path

from qpcr_pipeline.config import PanelConfig
from qpcr_pipeline.panel import PanelDefinition, PanelTarget, TargetGroup
from qpcr_pipeline.panel_manifest import (
    PanelProposal,
    approve_panel_proposal,
    write_panel_proposal,
)


def proposal_panel_config(target_name: str) -> PanelConfig:
    return PanelConfig(
        proposal=PanelDefinition(
            target=PanelTarget(
                name=target_name,
                taxid=None,
                mode="broad_detection",
                subtype=None,
                groups=(
                    TargetGroup(
                        name="all",
                        required=True,
                        dataset_roles=("DESIGN",),
                        reasons=("test_fixture",),
                        proposed_by=("manual",),
                    ),
                ),
            ),
            non_targets=(),
        )
    )


def approved_panel_config(directory: Path, target_name: str) -> PanelConfig:
    proposal_config = proposal_panel_config(target_name)
    assert proposal_config.proposal is not None
    proposal_path = Path(directory) / "test_panel_proposal.yaml"
    approved_path = Path(directory) / "test_approved_panel.json"
    write_panel_proposal(
        PanelProposal(
            schema_version=1,
            status="PROPOSED",
            definition=proposal_config.proposal,
        ),
        proposal_path,
    )
    approve_panel_proposal(proposal_path, approved_path)
    return PanelConfig(frozen_manifest=approved_path)
