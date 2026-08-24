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
class ClusteringConfig:
    enabled: bool = False
    identity: float = 0.95
    threads: int = 1
    memory_mb: int = 800


@dataclass(frozen=True, slots=True)
class AlignmentConfig:
    enabled: bool = False
    threads: int = 1
    reference_id: str | None = None


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
    clustering: ClusteringConfig = field(default_factory=ClusteringConfig)
    alignment: AlignmentConfig = field(default_factory=AlignmentConfig)

    @property
    def selected_input(
        self,
    ) -> tuple[Path, Literal["fasta", "genbank"]] | NcbiInputConfig:
        validate_pipeline_config(self)
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
    clustering_config = raw.get("clustering", {})
    clustering = _parse_clustering_config(clustering_config)
    alignment_config = raw.get("alignment", {})
    alignment = _parse_alignment_config(alignment_config)

    target_name = _required_string(target, "name", section="target")
    input_fasta = _optional_path(input_config, "fasta")
    input_genbank = _optional_path(input_config, "genbank")
    input_ncbi = _parse_ncbi_input(input_config["ncbi"]) if "ncbi" in input_config else None
    config = PipelineConfig(
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
        clustering=clustering,
        alignment=alignment,
    )
    validate_pipeline_config(config)
    return config


def validate_pipeline_config(config: PipelineConfig) -> None:
    """Validate a parsed or directly constructed pipeline configuration."""
    if not isinstance(config, PipelineConfig):
        raise ValueError("Pipeline configuration must be a PipelineConfig.")
    if not isinstance(config.target_name, str) or not config.target_name.strip():
        raise ValueError("Pipeline target_name must be a non-empty string.")
    for field_name, path in (
        ("input_fasta", config.input_fasta),
        ("input_genbank", config.input_genbank),
    ):
        if path is not None and not isinstance(path, Path):
            raise ValueError(f"Pipeline {field_name} must be a Path when configured.")
    if config.input_ncbi is not None and not isinstance(
        config.input_ncbi, NcbiInputConfig
    ):
        raise ValueError("Pipeline input_ncbi must be an NcbiInputConfig when configured.")
    source_count = sum(
        source is not None
        for source in (config.input_fasta, config.input_genbank, config.input_ncbi)
    )
    if source_count != 1:
        raise ValueError(
            "Exactly one local sequence input or NCBI input must be configured. "
            "Exactly one sequence input is allowed."
        )
    if not isinstance(config.qc, QCConfig):
        raise ValueError("Pipeline qc must be a QCConfig.")
    if not isinstance(config.clustering, ClusteringConfig):
        raise ValueError("Pipeline clustering must be a ClusteringConfig.")
    if not isinstance(config.alignment, AlignmentConfig):
        raise ValueError("Pipeline alignment must be an AlignmentConfig.")
    validate_clustering_config(config.clustering)
    validate_alignment_config(config.alignment)
    if config.input_ncbi is not None:
        validate_ncbi_input_config(config.input_ncbi)


def validate_ncbi_input_config(
    config: NcbiInputConfig,
) -> Literal["query", "accessions", "frozen"]:
    """Validate NCBI input regardless of whether it came from YAML or Python."""
    if not isinstance(config, NcbiInputConfig):
        raise ValueError("NCBI input must be an NcbiInputConfig.")

    query = config.query
    if query is not None and (not isinstance(query, str) or not query.strip()):
        raise ValueError("NCBI input query must be a non-blank string.")
    if not isinstance(config.accessions, tuple):
        raise ValueError("NCBI input accessions must be a tuple.")
    if any(
        not isinstance(accession, str) or not accession.strip()
        for accession in config.accessions
    ):
        raise ValueError("NCBI input accessions must contain only non-blank strings.")
    if len(set(config.accessions)) != len(config.accessions):
        raise ValueError("NCBI input accessions must be unique.")
    if config.frozen_dataset is not None and not isinstance(
        config.frozen_dataset, Path
    ):
        raise ValueError("NCBI input frozen_dataset must be a Path.")

    has_query = query is not None
    has_accessions = bool(config.accessions)
    has_frozen = config.frozen_dataset is not None
    if sum((has_query, has_accessions, has_frozen)) != 1:
        raise ValueError(
            "NCBI input must specify exactly one query, accession list, or frozen_dataset."
        )

    _validate_ncbi_config_integer(
        config.batch_size, "batch_size", minimum=1, maximum=500
    )
    _validate_ncbi_config_integer(config.retries, "retries", minimum=0, maximum=10)
    if config.max_records is not None and (
        isinstance(config.max_records, bool)
        or not isinstance(config.max_records, int)
        or config.max_records < 1
    ):
        raise ValueError("NCBI input max_records must be a positive integer.")
    if has_frozen:
        if config.batch_size != 100:
            raise ValueError(
                "NCBI frozen_dataset input cannot override the default batch_size."
            )
        if config.retries != 3:
            raise ValueError(
                "NCBI frozen_dataset input cannot override the default retries."
            )
        if config.max_records is not None:
            raise ValueError("NCBI frozen_dataset input cannot specify max_records.")
    elif config.max_records is not None and not has_query:
        raise ValueError(
            "NCBI accession input cannot specify max_records; it is only valid with a query."
        )

    if has_query:
        return "query"
    if has_accessions:
        return "accessions"
    return "frozen"


def validate_clustering_config(config: ClusteringConfig) -> None:
    if not isinstance(config, ClusteringConfig):
        raise ValueError("Clustering configuration must be a ClusteringConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Clustering enabled must be a boolean.")
    if (
        isinstance(config.identity, bool)
        or not isinstance(config.identity, (int, float))
        or not 0.80 <= config.identity <= 1.0
    ):
        raise ValueError("Clustering identity must be a number between 0.80 and 1.0.")
    if (
        isinstance(config.threads, bool)
        or not isinstance(config.threads, int)
        or not 1 <= config.threads <= 256
    ):
        raise ValueError("Clustering threads must be an integer between 1 and 256.")
    if (
        isinstance(config.memory_mb, bool)
        or not isinstance(config.memory_mb, int)
        or config.memory_mb < 1
    ):
        raise ValueError("Clustering memory_mb must be a positive integer.")


def validate_alignment_config(config: AlignmentConfig) -> None:
    if not isinstance(config, AlignmentConfig):
        raise ValueError("Alignment configuration must be an AlignmentConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Alignment enabled must be a boolean.")
    if (
        isinstance(config.threads, bool)
        or not isinstance(config.threads, int)
        or not 1 <= config.threads <= 256
    ):
        raise ValueError("Alignment threads must be an integer between 1 and 256.")
    if config.reference_id is not None and (
        not isinstance(config.reference_id, str) or not config.reference_id.strip()
    ):
        raise ValueError("Alignment reference_id must be a non-blank string when configured.")


def _validate_ncbi_config_integer(
    value: object, field_name: str, *, minimum: int, maximum: int
) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < minimum
        or value > maximum
    ):
        raise ValueError(
            f"NCBI input {field_name} must be an integer between {minimum} and {maximum}."
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


def _parse_clustering_config(raw: Any) -> ClusteringConfig:
    if not isinstance(raw, dict):
        raise ValueError("Configuration section 'clustering' must be a mapping.")
    allowed_fields = {"enabled", "identity", "threads", "memory_mb"}
    unknown_fields = set(raw) - allowed_fields
    if unknown_fields:
        rendered = ", ".join(sorted(str(field) for field in unknown_fields))
        raise ValueError(
            f"Configuration section 'clustering' fields {rendered} are unrecognized."
        )
    config = ClusteringConfig(
        enabled=raw.get("enabled", False),
        identity=raw.get("identity", 0.95),
        threads=raw.get("threads", 1),
        memory_mb=raw.get("memory_mb", 800),
    )
    validate_clustering_config(config)
    return config


def _parse_alignment_config(raw: Any) -> AlignmentConfig:
    if not isinstance(raw, dict):
        raise ValueError("Configuration section 'alignment' must be a mapping.")
    allowed_fields = {"enabled", "threads", "reference_id"}
    unknown_fields = set(raw) - allowed_fields
    if unknown_fields:
        rendered = ", ".join(sorted(str(field) for field in unknown_fields))
        raise ValueError(
            f"Configuration section 'alignment' fields {rendered} are unrecognized."
        )
    config = AlignmentConfig(
        enabled=raw.get("enabled", False),
        threads=raw.get("threads", 1),
        reference_id=raw.get("reference_id"),
    )
    validate_alignment_config(config)
    return config


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

    allowed_fields = {
        "query",
        "accessions",
        "frozen_dataset",
        "batch_size",
        "retries",
        "max_records",
    }
    unknown_fields = set(raw) - allowed_fields
    if unknown_fields:
        rendered = ", ".join(sorted((str(field) for field in unknown_fields)))
        raise ValueError(
            f"Configuration section 'input.ncbi' fields {rendered} are unrecognized."
        )

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

    config = NcbiInputConfig(
        query=query,
        accessions=accessions or (),
        frozen_dataset=frozen_dataset,
        batch_size=batch_size,
        retries=retries,
        max_records=max_records,
    )
    validate_ncbi_input_config(config)
    return config


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
    if key not in raw:
        return default
    value = raw[key]
    if value is None:
        if default is None:
            return None
        raise ValueError(f"Configuration value 'input.ncbi.{key}' must be an integer.")
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
