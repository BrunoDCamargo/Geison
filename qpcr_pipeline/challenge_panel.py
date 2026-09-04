"""Resolve approved CHALLENGE panel entries to configured off-target datasets."""

from __future__ import annotations

from dataclasses import dataclass

from qpcr_pipeline.config import OffTargetConfig
from qpcr_pipeline.off_targets import OffTargetDataset, load_off_target_dataset
from qpcr_pipeline.panel import Criticality
from qpcr_pipeline.panel_manifest import ApprovedPanelManifest


@dataclass(frozen=True, slots=True)
class ChallengeDatasetBinding:
    name: str
    criticality: Criticality
    dataset: OffTargetDataset


def _normalized_name(value: str) -> str:
    return value.strip().casefold()


def resolve_challenge_datasets(
    manifest: ApprovedPanelManifest,
    configs: tuple[OffTargetConfig, ...],
) -> tuple[ChallengeDatasetBinding, ...]:
    """Resolve approved CHALLENGE entries in panel order.

    The approved panel owns semantic membership and criticality. ``off_targets``
    only provides the physical datasets. A contrastive run with no approved
    CHALLENGE entries is invalid rather than silently completing without
    contrast evidence.
    """
    if not isinstance(manifest, ApprovedPanelManifest):
        raise ValueError("Challenge resolution requires an approved panel manifest.")
    if not isinstance(configs, tuple):
        raise ValueError("Off-target configurations must be a tuple.")

    by_name: dict[str, OffTargetConfig] = {}
    for config in configs:
        if not isinstance(config, OffTargetConfig):
            raise ValueError("Off-target configurations must contain OffTargetConfig values.")
        key = _normalized_name(config.name)
        if key in by_name:
            raise ValueError(
                "Off-target dataset names must be unique after normalization."
            )
        by_name[key] = config

    challenge_entries = tuple(
        item
        for item in manifest.definition.non_targets
        if "CHALLENGE" in item.dataset_roles
    )
    if not challenge_entries:
        raise ValueError(
            "Enabled contrastive conservation requires at least one approved CHALLENGE dataset."
        )

    resolved: list[ChallengeDatasetBinding] = []
    for item in challenge_entries:
        config = by_name.get(_normalized_name(item.name))
        if config is None:
            raise ValueError(
                f"Challenge dataset missing for approved panel entry {item.name!r}."
            )
        resolved.append(
            ChallengeDatasetBinding(
                name=item.name,
                criticality=item.criticality,
                dataset=load_off_target_dataset(config),
            )
        )
    return tuple(resolved)
