"""Extended pipeline configuration with contrastive-conservation support.

The original configuration implementation is kept in ``legacy`` so existing
parsing and validation behavior stays unchanged while this branch adds the
new stage contract.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import yaml

from . import legacy as _legacy


for _name in dir(_legacy):
    if not _name.startswith("_"):
        globals()[_name] = getattr(_legacy, _name)

_LegacyPipelineConfig = _legacy.PipelineConfig
_legacy_load_config = _legacy.load_config
_legacy_validate_pipeline_config = _legacy.validate_pipeline_config


@dataclass(frozen=True, slots=True)
class ContrastiveConservationConfig:
    enabled: bool = False


@dataclass(frozen=True, slots=True)
class PipelineConfig(_LegacyPipelineConfig):
    contrastive_conservation: ContrastiveConservationConfig = field(
        default_factory=ContrastiveConservationConfig
    )


def validate_contrastive_conservation_config(
    config: ContrastiveConservationConfig,
) -> None:
    if not isinstance(config, ContrastiveConservationConfig):
        raise ValueError(
            "Contrastive conservation configuration must be a "
            "ContrastiveConservationConfig."
        )
    if not isinstance(config.enabled, bool):
        raise ValueError("contrastive conservation enabled must be a boolean.")


def _parse_contrastive_conservation_config(raw: Any) -> ContrastiveConservationConfig:
    if not isinstance(raw, dict):
        raise ValueError(
            "Configuration section 'contrastive_conservation' must be a mapping."
        )
    allowed_fields = {"enabled"}
    unknown_fields = set(raw) - allowed_fields
    if unknown_fields:
        rendered = ", ".join(sorted(str(field) for field in unknown_fields))
        raise ValueError(
            "Configuration section 'contrastive_conservation' fields "
            f"{rendered} are unrecognized."
        )
    config = ContrastiveConservationConfig(enabled=raw.get("enabled", False))
    validate_contrastive_conservation_config(config)
    return config


def validate_pipeline_config(config: PipelineConfig) -> None:
    _legacy_validate_pipeline_config(config)
    validate_contrastive_conservation_config(config.contrastive_conservation)
    if config.contrastive_conservation.enabled and not config.conservation.enabled:
        raise ValueError(
            "Enabled contrastive conservation requires enabled conservation."
        )
    if config.contrastive_conservation.enabled:
        if config.panel is None or config.panel.frozen_manifest is None:
            raise ValueError(
                "Enabled contrastive conservation requires an approved frozen panel."
            )
        if not config.off_targets:
            raise ValueError(
                "Enabled contrastive conservation requires off-target datasets."
            )


def load_config(path: str | Path) -> PipelineConfig:
    base = _legacy_load_config(path)
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("Configuration root must be a mapping.")
    contrastive = _parse_contrastive_conservation_config(
        raw.get("contrastive_conservation", {})
    )
    config = replace(base, contrastive_conservation=contrastive)
    validate_pipeline_config(config)
    return config


# The legacy functions resolve these names dynamically from their module globals.
# Point them at the extended class/validator so inherited ``selected_input`` and
# legacy parsing continue to validate the new configuration shape.
_legacy.PipelineConfig = PipelineConfig
_legacy.validate_pipeline_config = validate_pipeline_config

__all__ = [
    name for name in globals()
    if not name.startswith("_") and name not in {"Any", "Path", "yaml"}
]
