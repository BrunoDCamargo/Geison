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
class NcbiInputConfig:
    query: str | None = None
    accessions: tuple[str, ...] = ()
    frozen_dataset: Path | None = None
    batch_size: int = 100
    retries: int = 3
    max_records: int | None = None


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path | None = None
    input_genbank: Path | None = None
    input_ncbi: NcbiInputConfig | None = None
    qc: QCConfig = field(default_factory=QCConfig)

    @property
    def selected_input(
        self,
    ) -> tuple[Path, Literal["fasta", "genbank"]] | NcbiInputConfig:
        if self.input_fasta is not None:
            return self.input_fasta, "fasta"
        if self.input_genbank is not None:
            return self.input_genbank, "genbank"
        if self.input_ncbi is not None:
            return self.input_ncbi
        raise ValueError("Exactly one sequence input must be configured.")


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
    input_ncbi = _parse_ncbi_input(input_config["ncbi"]) if "ncbi" in input_config else None
    input_sources = sum(
        source is not None for source in (input_fasta, input_genbank, input_ncbi)
    )
    if input_sources != 1:
        raise ValueError("Exactly one local sequence input must be configured.")

    return PipelineConfig(
        target_name=target_name,
        input_fasta=input_fasta,
        input_genbank=input_genbank,
        input_ncbi=input_ncbi,
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


def _parse_ncbi_input(raw: Any) -> NcbiInputConfig:
    if not isinstance(raw, dict):
        raise ValueError("Configuration section 'input.ncbi' must be a mapping.")

    query = _optional_ncbi_string(raw, "query")
    accessions = _optional_accessions(raw)
    frozen_dataset = _optional_ncbi_path(raw, "frozen_dataset")
    modes = sum(value is not None for value in (query, accessions, frozen_dataset))
    if modes != 1:
        raise ValueError(
            "Configuration value 'input.ncbi' must specify exactly one of "
            "query, accessions, or frozen_dataset."
        )

    batch_size = _ncbi_integer(raw, "batch_size", default=100, minimum=1, maximum=500)
    retries = _ncbi_integer(raw, "retries", default=3, minimum=0, maximum=10)
    max_records = _ncbi_integer(raw, "max_records", default=None, minimum=1, maximum=None)

    if frozen_dataset is not None:
        for key in ("batch_size", "retries", "max_records"):
            if key in raw:
                raise ValueError(
                    f"Configuration value 'input.ncbi.{key}' cannot be used with "
                    "input.ncbi.frozen_dataset."
                )
    if max_records is not None and query is None:
        raise ValueError(
            "Configuration value 'input.ncbi.max_records' is only valid with input.ncbi.query."
        )

    return NcbiInputConfig(
        query=query,
        accessions=accessions or (),
        frozen_dataset=frozen_dataset,
        batch_size=batch_size,
        retries=retries,
        max_records=max_records,
    )


def _optional_ncbi_string(raw: dict[str, Any], key: str) -> str | None:
    value = raw.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"Configuration value 'input.ncbi.{key}' must be a non-empty string.")
    return value


def _optional_ncbi_path(raw: dict[str, Any], key: str) -> Path | None:
    value = _optional_ncbi_string(raw, key)
    return Path(value) if value is not None else None


def _optional_accessions(raw: dict[str, Any]) -> tuple[str, ...] | None:
    if "accessions" not in raw or raw["accessions"] is None:
        return None
    value = raw["accessions"]
    if (
        not isinstance(value, list)
        or not value
        or any(not isinstance(accession, str) or not accession.strip() for accession in value)
        or len(set(value)) != len(value)
    ):
        raise ValueError(
            "Configuration value 'input.ncbi.accessions' must be a non-empty list "
            "of unique non-empty strings."
        )
    return tuple(value)


def _ncbi_integer(
    raw: dict[str, Any],
    key: str,
    *,
    default: int | None,
    minimum: int,
    maximum: int | None,
) -> int | None:
    value = raw.get(key, default)
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"Configuration value 'input.ncbi.{key}' must be an integer.")
    if value < minimum or (maximum is not None and value > maximum):
        if maximum is None:
            bounds = f"at least {minimum}"
        else:
            bounds = f"between {minimum} and {maximum}"
        raise ValueError(f"Configuration value 'input.ncbi.{key}' must be {bounds}.")
    return value


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
