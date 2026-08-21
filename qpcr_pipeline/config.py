from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

import yaml


@dataclass(frozen=True, slots=True)
class QCConfig:
    min_length: int | None = None
    max_ambiguous_fraction: float | None = None
    expected_length: int | None = None
    length_tolerance_fraction: float | None = None


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path | None = None
    input_genbank: Path | None = None
    qc: QCConfig = field(default_factory=QCConfig)

    @property
    def selected_input(self) -> tuple[Path, Literal["fasta", "genbank"]]:
        if self.input_fasta is not None:
            return self.input_fasta, "fasta"
        if self.input_genbank is not None:
            return self.input_genbank, "genbank"
        raise ValueError("Exactly one local sequence input must be configured.")


def load_config(path: str | Path) -> PipelineConfig:
    config_path = Path(path)
    raw = yaml.safe_load(config_path.read_text(encoding="utf-8"))

    if not isinstance(raw, dict):
        raise ValueError("Configuration root must be a mapping.")

    target = _mapping(raw, "target")
    input_config = _mapping(raw, "input")
    qc_config = raw.get("qc", {})
    if not isinstance(qc_config, dict):
        raise ValueError("Configuration section 'qc' must be a mapping.")

    target_name = _required_string(target, "name", section="target")
    input_fasta = _optional_path(input_config, "fasta")
    input_genbank = _optional_path(input_config, "genbank")
    if (input_fasta is None) == (input_genbank is None):
        raise ValueError("Exactly one local sequence input must be configured.")

    return PipelineConfig(
        target_name=target_name,
        input_fasta=input_fasta,
        input_genbank=input_genbank,
        qc=QCConfig(
            min_length=_optional_integer(qc_config, "min_length"),
            max_ambiguous_fraction=_optional_number(qc_config, "max_ambiguous_fraction"),
            expected_length=_optional_integer(qc_config, "expected_length"),
            length_tolerance_fraction=_optional_number(qc_config, "length_tolerance_fraction"),
        ),
    )


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


def _optional_path(raw: dict[str, Any], key: str) -> Path | None:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Configuration value 'input.{key}' must be a non-empty string.")
    return Path(value)


def _optional_number(raw: dict[str, Any], key: str) -> int | float | None:
    value = raw.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"Configuration value 'qc.{key}' must be a number.")
    return value


def _optional_integer(raw: dict[str, Any], key: str) -> int | None:
    value = raw.get(key)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Configuration value 'qc.{key}' must be an integer.")
    return value
