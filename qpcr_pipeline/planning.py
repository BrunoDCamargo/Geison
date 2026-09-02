from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from qpcr_pipeline.checkpoint_stages import (
    STAGE_DEFINITIONS,
    SubprocessToolIdentityProvider,
    ToolIdentityProvider,
    stage_request,
)
from qpcr_pipeline.checkpointing import CheckpointManager, CheckpointManifest
from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.execution import (
    STAGE_ORDER,
    ExecutionPolicy,
    StageDecision,
    plan_from_validity,
    required_reuse_boundary,
    transitive_descendants,
)


@dataclass(frozen=True, slots=True)
class PipelineExecutionPlan:
    decisions: tuple[StageDecision, ...]
    reused_results: Mapping[str, object]
    reused_manifests: Mapping[str, CheckpointManifest]


def plan_pipeline(
    config: PipelineConfig,
    outdir: str | Path | None,
    *,
    execution: ExecutionPolicy | None = None,
    tool_identity_provider: ToolIdentityProvider | None = None,
) -> PipelineExecutionPlan:
    policy = execution or ExecutionPolicy()
    if outdir is None or (not policy.resume and policy.from_step is None):
        return PipelineExecutionPlan(
            decisions=plan_from_validity(policy, {}),
            reused_results={},
            reused_manifests={},
        )

    manager = CheckpointManager(Path(outdir))
    provider = tool_identity_provider or SubprocessToolIdentityProvider()
    results: dict[str, object] = {}
    manifests: dict[str, CheckpointManifest] = {}
    reusable: dict[str, bool] = {}

    if policy.from_step is not None:
        required = set(required_reuse_boundary(policy.from_step))
        blocked: set[str] = set()
        for stage in STAGE_ORDER:
            if stage not in required:
                continue
            if stage in blocked:
                reusable[stage] = False
                continue
            try:
                request = stage_request(stage, config, manifests, results, provider)
            except Exception:
                reusable[stage] = False
                blocked.update(transitive_descendants(stage))
                continue
            validation = manager.validate(request, STAGE_DEFINITIONS[stage].codec)
            reusable[stage] = validation.valid
            if not validation.valid:
                blocked.update(transitive_descendants(stage))
                continue
            assert validation.loaded is not None
            results[stage] = validation.loaded.state
            manifests[stage] = validation.loaded.manifest

        return PipelineExecutionPlan(
            decisions=plan_from_validity(policy, reusable),
            reused_results=dict(results),
            reused_manifests=dict(manifests),
        )

    blocked: set[str] = set()
    for stage in STAGE_ORDER:
        if stage in blocked:
            reusable[stage] = False
            continue
        try:
            request = stage_request(stage, config, manifests, results, provider)
        except Exception:
            reusable[stage] = False
            blocked.update(transitive_descendants(stage))
            continue
        validation = manager.validate(request, STAGE_DEFINITIONS[stage].codec)
        reusable[stage] = validation.valid
        if not validation.valid:
            blocked.update(transitive_descendants(stage))
            continue
        assert validation.loaded is not None
        results[stage] = validation.loaded.state
        manifests[stage] = validation.loaded.manifest

    decisions = plan_from_validity(policy, reusable)
    reused_stages = {item.stage for item in decisions if item.action == "REUSE"}
    return PipelineExecutionPlan(
        decisions=decisions,
        reused_results={stage: results[stage] for stage in reused_stages if stage in results},
        reused_manifests={stage: manifests[stage] for stage in reused_stages if stage in manifests},
    )
