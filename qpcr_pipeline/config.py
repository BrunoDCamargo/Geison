from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path


def load_config(path: str | Path) -> PipelineConfig:
    config_path = Path(path)
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    if not isinstance(raw, dict):
        raise ValueError("Configuration root must be a mapping.")

    target = _mapping(raw, "target")
    input_config = _mapping(raw, "input")

    target_name = _required_string(target, "name", section="target")
    input_fasta = Path(_required_string(input_config, "fasta", section="input"))

    return PipelineConfig(target_name=target_name, input_fasta=input_fasta)


def _mapping(raw: dict[str, Any], key: str) -> dict[str, Any]:
    value = raw.get(key)
    if not isinstance(value, dict):
        raise ValueError(f"Configuration section '{key}' must be a mapping.")
    return value


def _required_string(raw: dict[str, Any], key: str, *, section: str) -> str:
    value = raw.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Configuration value '{section}.{key}' must be a non-empty string.")
    return value
