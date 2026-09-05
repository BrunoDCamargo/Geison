"""Low-intrusion guided project setup backed by Geison's NCBI acquisition."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
from typing import Mapping

import yaml

from qpcr_pipeline.config import NcbiInputConfig
from qpcr_pipeline.ncbi import NcbiClient, acquire_ncbi_dataset
from qpcr_pipeline.panel_manifest import load_approved_panel_manifest


KNOWLEDGE_VERSION = "2026-09-05"


@dataclass(frozen=True, slots=True)
class GuidedChallengePreset:
    name: str
    criticality: str
    reason: str
    query: str
    max_records: int = 20


@dataclass(frozen=True, slots=True)
class GuidedTargetPreset:
    name: str
    query: str
    max_records: int
    challenges: tuple[GuidedChallengePreset, ...]
    syndrome: str | None = None
    vector: str | None = None


_WNV = GuidedTargetPreset(
    name="West Nile virus",
    query='"West Nile virus"[Organism] AND complete genome[Title]',
    max_records=50,
    challenges=(
        GuidedChallengePreset(
            name="Usutu virus",
            criticality="CRITICAL",
            reason="phylogenetic_neighbor",
            query='"Usutu virus"[Organism] AND complete genome[Title]',
        ),
        GuidedChallengePreset(
            name="Japanese encephalitis virus",
            criticality="CRITICAL",
            reason="phylogenetic_neighbor",
            query='"Japanese encephalitis virus"[Organism] AND complete genome[Title]',
        ),
        GuidedChallengePreset(
            name="Dengue virus",
            criticality="IMPORTANT",
            reason="clinical_differential",
            query='"Dengue virus"[Organism] AND complete genome[Title]',
        ),
    ),
    syndrome="arboviral febrile disease",
    vector="mosquito",
)

_PRESETS: Mapping[str, GuidedTargetPreset] = {
    _WNV.name.casefold(): _WNV,
}


def supported_guided_targets() -> tuple[str, ...]:
    """Return canonical targets with curated guided panel knowledge."""
    return tuple(preset.name for preset in _PRESETS.values())


def _preset(target_name: str) -> GuidedTargetPreset:
    if not isinstance(target_name, str) or not target_name.strip():
        raise ValueError("Guided target must be a non-blank string.")
    preset = _PRESETS.get(target_name.strip().casefold())
    if preset is None:
        supported = ", ".join(supported_guided_targets())
        raise ValueError(
            f"Unsupported guided target {target_name!r}. Supported guided targets: {supported}."
        )
    return preset


def _slug(value: str) -> str:
    slug = re.sub(r"[^a-z0-9]+", "-", value.casefold()).strip("-")
    return slug or "dataset"


def _source_marker() -> str:
    return f"geison_guided_knowledge@{KNOWLEDGE_VERSION}"


def _challenge_relative_path(index: int, challenge: GuidedChallengePreset) -> str:
    return f"guided_challenges/{index:03d}-{_slug(challenge.name)}"


def _ncbi_mapping(query: str, max_records: int) -> dict[str, object]:
    return {
        "query": query,
        "batch_size": 20,
        "retries": 5,
        "max_records": max_records,
    }


def build_guided_proposal_config(target_name: str) -> dict[str, object]:
    """Build a normal pipeline proposal configuration without network access."""
    preset = _preset(target_name)
    proposed_by = [_source_marker()]
    non_targets = [
        {
            "name": challenge.name,
            "taxid": None,
            "criticality": challenge.criticality,
            "dataset_roles": ["CHALLENGE"],
            "reasons": [challenge.reason],
            "proposed_by": proposed_by,
            "sequence_selection": [],
        }
        for challenge in preset.challenges
    ]
    off_targets = [
        {
            "name": challenge.name,
            "frozen_dataset": _challenge_relative_path(index, challenge),
        }
        for index, challenge in enumerate(preset.challenges, 1)
    ]
    return {
        "target": {"name": preset.name},
        "input": {"ncbi": _ncbi_mapping(preset.query, preset.max_records)},
        "alignment": {"enabled": True, "threads": 2},
        "conservation": {"enabled": True, "window_size": 100, "step_size": 20},
        "primer_design": {
            "enabled": True,
            "max_candidate_regions": 6,
            "assays_per_region": 3,
        },
        "inclusivity": {"enabled": True},
        "off_targets": off_targets,
        "specificity": {"enabled": True},
        "ranking": {"enabled": True},
        "panel": {
            "proposal": {
                "target": {
                    "name": preset.name,
                    "taxid": None,
                    "mode": "broad_detection",
                    "subtype": None,
                    "groups": [
                        {
                            "name": "broad target diversity",
                            "required": True,
                            "dataset_roles": ["DESIGN"],
                            "reasons": ["broad_detection_target"],
                            "proposed_by": proposed_by,
                            "sequence_selection": [],
                        }
                    ],
                },
                "non_targets": non_targets,
                "diagnostic_context": {
                    "syndrome": preset.syndrome,
                    "geography": None,
                    "sample_type": None,
                    "vector": preset.vector,
                },
            }
        },
        "contrastive_conservation": {"enabled": False},
    }


def _approved_challenge_presets(
    preset: GuidedTargetPreset,
    approved_panel_path: Path,
) -> tuple[tuple[object, GuidedChallengePreset], ...]:
    manifest = load_approved_panel_manifest(approved_panel_path)
    if manifest.definition.target.name.strip().casefold() != preset.name.casefold():
        raise ValueError(
            f"Approved panel target {manifest.definition.target.name!r} does not match "
            f"guided target {preset.name!r}."
        )
    by_name = {challenge.name.casefold(): challenge for challenge in preset.challenges}
    selected: list[tuple[object, GuidedChallengePreset]] = []
    for item in manifest.definition.non_targets:
        if "CHALLENGE" not in item.dataset_roles:
            continue
        challenge = by_name.get(item.name.strip().casefold())
        if challenge is None:
            raise ValueError(
                f"Approved challenge {item.name!r} is not available in guided knowledge "
                f"version {KNOWLEDGE_VERSION}."
            )
        selected.append((item, challenge))
    if not selected:
        raise ValueError("Approved guided panel must contain at least one CHALLENGE non-target.")
    return tuple(selected)


def _dataset_manifest_summary(dataset_dir: Path) -> tuple[int, str, str]:
    manifest_path = dataset_dir / "dataset_manifest.json"
    records_path = dataset_dir / "records.gb"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    consolidated = manifest.get("consolidated")
    if not isinstance(consolidated, dict):
        raise ValueError(f"Guided NCBI dataset {dataset_dir} has no consolidated manifest.")
    count = consolidated.get("record_count")
    if type(count) is not int or count < 1:
        raise ValueError(f"Guided NCBI dataset {dataset_dir} has invalid record_count.")
    return (
        count,
        hashlib.sha256(records_path.read_bytes()).hexdigest(),
        hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
    )


def finalize_guided_project(
    target_name: str,
    approved_panel_path: str | Path,
    workspace: str | Path,
    *,
    ncbi_client: NcbiClient | None = None,
) -> Path:
    """Freeze approved guided challenge data and write a standard run config."""
    preset = _preset(target_name)
    approved_path = Path(approved_panel_path).resolve()
    project_dir = Path(workspace).resolve()
    project_dir.mkdir(parents=True, exist_ok=True)
    selected = _approved_challenge_presets(preset, approved_path)

    challenge_root = project_dir / "guided_challenges"
    challenge_root.mkdir(parents=True, exist_ok=True)
    off_targets: list[dict[str, object]] = []
    acquisition_rows: list[dict[str, object]] = []

    for index, (panel_item, challenge) in enumerate(selected, 1):
        dataset_dir = challenge_root / f"{index:03d}-{_slug(challenge.name)}"
        acquired = acquire_ncbi_dataset(
            NcbiInputConfig(
                query=challenge.query,
                batch_size=20,
                retries=5,
                max_records=challenge.max_records,
            ),
            dataset_dir,
            client=ncbi_client,
        )
        record_count, records_sha256, manifest_sha256 = _dataset_manifest_summary(dataset_dir)
        off_targets.append(
            {
                "name": challenge.name,
                "frozen_dataset": str(dataset_dir),
            }
        )
        acquisition_rows.append(
            {
                "name": challenge.name,
                "criticality": getattr(panel_item, "criticality"),
                "query": challenge.query,
                "frozen_dataset": str(dataset_dir),
                "record_count": record_count,
                "records_sha256": records_sha256,
                "manifest_sha256": manifest_sha256,
            }
        )
        if acquired.records_path != dataset_dir / "records.gb":
            raise ValueError("Guided challenge acquisition returned an unexpected records path.")

    approved_config = deepcopy(build_guided_proposal_config(preset.name))
    approved_config["panel"] = {"frozen_manifest": str(approved_path)}
    approved_config["contrastive_conservation"] = {"enabled": True}
    approved_config["off_targets"] = off_targets

    config_path = project_dir / "config-approved.yaml"
    config_path.write_text(
        yaml.safe_dump(approved_config, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )

    acquisition_manifest = {
        "schema_version": 1,
        "knowledge_version": KNOWLEDGE_VERSION,
        "target": preset.name,
        "approved_panel_sha256": hashlib.sha256(approved_path.read_bytes()).hexdigest(),
        "datasets": acquisition_rows,
        "limitations": [
            "Guided max_records values bound the workbench acquisition and are not a representative-sampling claim."
        ],
    }
    (project_dir / "guided_acquisition_manifest.json").write_text(
        json.dumps(acquisition_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return config_path
