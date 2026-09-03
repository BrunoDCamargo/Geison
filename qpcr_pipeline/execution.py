"""Deterministic execution planning for resumable Geison pipeline stages."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Mapping


StageName = Literal[
    "panel",
    "input",
    "qc",
    "clustering",
    "alignment",
    "conservation",
    "primer_design",
    "inclusivity",
    "specificity",
    "ranking",
]
StageAction = Literal["RUN", "REUSE", "FORCED"]

STAGE_ORDER: tuple[StageName, ...] = (
    "panel",
    "input",
    "qc",
    "clustering",
    "alignment",
    "conservation",
    "primer_design",
    "inclusivity",
    "specificity",
    "ranking",
)

STAGE_DEPENDENCIES: dict[StageName, tuple[StageName, ...]] = {
    "panel": (),
    "input": ("panel",),
    "qc": ("input",),
    "clustering": ("qc",),
    "alignment": ("clustering",),
    "conservation": ("alignment",),
    "primer_design": ("conservation",),
    "inclusivity": ("primer_design", "qc"),
    "specificity": ("primer_design",),
    "ranking": ("primer_design", "inclusivity", "specificity"),
}

_STAGE_SET = frozenset(STAGE_ORDER)


@dataclass(frozen=True, slots=True)
class ExecutionPolicy:
    resume: bool = False
    from_step: str | None = None
    force_step: str | None = None

    def __post_init__(self) -> None:
        if self.resume and self.from_step is not None:
            raise ValueError("--resume cannot be combined with --from-step")
        if self.from_step is not None and self.force_step is not None:
            raise ValueError("--from-step cannot be combined with --force-step")
        if self.force_step is not None and not self.resume:
            raise ValueError("--force-step requires --resume")
        for value in (self.from_step, self.force_step):
            if value is not None and value not in _STAGE_SET:
                raise ValueError(f"unknown pipeline stage: {value}")


@dataclass(frozen=True, slots=True)
class StageDecision:
    stage: StageName
    action: StageAction
    reason: str


def _validate_stage(stage: str) -> StageName:
    if stage not in _STAGE_SET:
        raise ValueError(f"unknown pipeline stage: {stage}")
    return stage  # type: ignore[return-value]


def transitive_descendants(stage: str) -> tuple[StageName, ...]:
    root = _validate_stage(stage)
    found: set[StageName] = set()
    changed = True
    while changed:
        changed = False
        for candidate in STAGE_ORDER:
            if candidate == root or candidate in found:
                continue
            dependencies = STAGE_DEPENDENCIES[candidate]
            if root in dependencies or any(dependency in found for dependency in dependencies):
                found.add(candidate)
                changed = True
    return tuple(candidate for candidate in STAGE_ORDER if candidate in found)


def _transitive_dependencies(stages: set[StageName]) -> set[StageName]:
    required: set[StageName] = set()
    frontier = list(stages)
    while frontier:
        current = frontier.pop()
        for dependency in STAGE_DEPENDENCIES[current]:
            if dependency not in required and dependency not in stages:
                required.add(dependency)
                frontier.append(dependency)
            elif dependency in stages:
                frontier.append(dependency)
    return required


def required_reuse_boundary(stage: str) -> tuple[StageName, ...]:
    root = _validate_stage(stage)
    forced = {root, *transitive_descendants(root)}
    required = _transitive_dependencies(set(forced))
    return tuple(candidate for candidate in STAGE_ORDER if candidate in required)


def plan_from_validity(
    policy: ExecutionPolicy,
    reusable: Mapping[str, bool],
) -> tuple[StageDecision, ...]:
    if policy.from_step is not None:
        root = _validate_stage(policy.from_step)
        forced = {root, *transitive_descendants(root)}
        boundary = set(required_reuse_boundary(root))
        blocking = [
            stage for stage in STAGE_ORDER
            if stage in boundary and not bool(reusable.get(stage, False))
        ]
        if blocking:
            names = ", ".join(blocking)
            raise ValueError(
                "--from-step requires valid reusable checkpoints for: " + names
            )
        decisions = []
        for stage in STAGE_ORDER:
            if stage in forced:
                decisions.append(StageDecision(stage, "FORCED", "selected --from-step subgraph"))
            elif stage in boundary:
                decisions.append(StageDecision(stage, "REUSE", "required valid boundary checkpoint"))
            else:
                decisions.append(StageDecision(stage, "REUSE", "unaffected prerequisite checkpoint"))
        return tuple(decisions)

    if not policy.resume:
        return tuple(
            StageDecision(stage, "RUN", "normal execution recomputes every stage")
            for stage in STAGE_ORDER
        )

    forced: set[StageName] = set()
    if policy.force_step is not None:
        root = _validate_stage(policy.force_step)
        forced.add(root)
        forced.update(transitive_descendants(root))

    decisions: list[StageDecision] = []
    for stage in STAGE_ORDER:
        if stage in forced:
            decisions.append(
                StageDecision(
                    stage,
                    "FORCED",
                    "explicitly forced or dependent on a forced/recomputed stage",
                )
            )
            continue
        if bool(reusable.get(stage, False)):
            decisions.append(StageDecision(stage, "REUSE", "checkpoint is valid"))
            continue
        decisions.append(StageDecision(stage, "RUN", "checkpoint is invalid or missing"))
        forced.update(transitive_descendants(stage))

    return tuple(decisions)
