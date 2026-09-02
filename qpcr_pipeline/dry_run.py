from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from qpcr_pipeline.checkpoint_stages import ToolIdentityProvider
from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.diagnostics import EnvironmentInspector, EnvironmentReport
from qpcr_pipeline.execution import ExecutionPolicy, StageDecision
from qpcr_pipeline.planning import plan_pipeline


@dataclass(frozen=True, slots=True)
class DryRunReport:
    target_name: str
    decisions: tuple[StageDecision, ...]
    environment: EnvironmentReport


class _EnvironmentToolIdentityProvider:
    def __init__(self, environment: EnvironmentReport) -> None:
        self.environment = environment

    def identity(self, tool_name: str) -> Mapping[str, str]:
        component = self.environment.tools.get(tool_name)
        if component is None or not component.installed or not component.version:
            raise RuntimeError(
                f"Required external tool {tool_name!r} is unavailable in the dry-run environment report."
            )
        return {"name": tool_name, "version": component.version}


def dry_run_pipeline(
    config: PipelineConfig,
    outdir: str | Path | None = None,
    *,
    execution: ExecutionPolicy | None = None,
    inspector: EnvironmentInspector | None = None,
    tool_identity_provider: ToolIdentityProvider | None = None,
) -> DryRunReport:
    environment = (inspector or EnvironmentInspector()).inspect(config)
    provider = tool_identity_provider or _EnvironmentToolIdentityProvider(environment)
    plan = plan_pipeline(
        config,
        outdir,
        execution=execution,
        tool_identity_provider=provider,
    )
    return DryRunReport(
        target_name=config.target_name,
        decisions=plan.decisions,
        environment=environment,
    )
