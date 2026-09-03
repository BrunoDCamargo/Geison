from qpcr_pipeline.config import PanelConfig
from qpcr_pipeline.panel import PanelDefinition, PanelTarget, TargetGroup


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
