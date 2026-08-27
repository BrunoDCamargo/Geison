"""Offline off-target dataset loading and provenance for specificity evaluation."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Literal

from qpcr_pipeline.config import OffTargetConfig, validate_off_target_config
from qpcr_pipeline.local_input import LocalSequenceRecord, load_fasta, load_genbank
from qpcr_pipeline.ncbi import validate_frozen_dataset


@dataclass(frozen=True, slots=True)
class OffTargetDataset:
    name: str
    source_type: Literal["FASTA", "NCBI_FROZEN"]
    source_path: Path
    sha256: str
    sequence_ids: tuple[str, ...]
    records: tuple[LocalSequenceRecord, ...]
    frozen_manifest_path: Path | None
    frozen_manifest: dict[str, object] | None


def _sha256(path: Path, *, dataset_name: str) -> str:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError as error:
        raise ValueError(
            f"Off-target dataset {dataset_name!r} could not read {path}."
        ) from error


def _validated_record_ids(
    records: tuple[LocalSequenceRecord, ...], *, dataset_name: str
) -> tuple[str, ...]:
    ids: list[str] = []
    for index, record in enumerate(records, 1):
        if not isinstance(record, LocalSequenceRecord):
            raise ValueError(
                f"Off-target dataset {dataset_name!r} record {index} is invalid."
            )
        if not isinstance(record.sequence_id, str) or not record.sequence_id.strip():
            raise ValueError(
                f"Off-target dataset {dataset_name!r} record {index} must have a non-blank ID."
            )
        ids.append(record.sequence_id)
    if len(set(ids)) != len(ids):
        raise ValueError(
            f"Off-target dataset {dataset_name!r} sequence IDs must be unique."
        )
    return tuple(ids)


def _load_fasta_dataset(config: OffTargetConfig) -> OffTargetDataset:
    assert config.fasta is not None
    path = config.fasta
    try:
        raw = path.read_bytes()
        records = load_fasta(path)
    except (OSError, ValueError) as error:
        raise ValueError(
            f"Off-target dataset {config.name!r} FASTA {path} is unreadable or invalid."
        ) from error
    if raw.strip() and not records:
        raise ValueError(
            f"Off-target dataset {config.name!r} FASTA {path} is invalid."
        )
    sequence_ids = _validated_record_ids(records, dataset_name=config.name)
    return OffTargetDataset(
        name=config.name,
        source_type="FASTA",
        source_path=path,
        sha256=hashlib.sha256(raw).hexdigest(),
        sequence_ids=sequence_ids,
        records=records,
        frozen_manifest_path=None,
        frozen_manifest=None,
    )


def _load_frozen_dataset(config: OffTargetConfig) -> OffTargetDataset:
    assert config.frozen_dataset is not None
    directory = config.frozen_dataset
    try:
        acquired = validate_frozen_dataset(directory)
        records = load_genbank(acquired.records_path)
        manifest = json.loads(acquired.manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as error:
        raise ValueError(
            f"Off-target dataset {config.name!r} frozen dataset {directory} is invalid."
        ) from error
    if not isinstance(manifest, dict):
        raise ValueError(
            f"Off-target dataset {config.name!r} frozen manifest must be a mapping."
        )
    sequence_ids = _validated_record_ids(records, dataset_name=config.name)
    return OffTargetDataset(
        name=config.name,
        source_type="NCBI_FROZEN",
        source_path=directory,
        sha256=_sha256(acquired.records_path, dataset_name=config.name),
        sequence_ids=sequence_ids,
        records=records,
        frozen_manifest_path=acquired.manifest_path,
        frozen_manifest=manifest,
    )


def load_off_target_dataset(config: OffTargetConfig) -> OffTargetDataset:
    """Load exactly one configured dataset without any network acquisition."""
    validate_off_target_config(config)
    if config.fasta is not None:
        return _load_fasta_dataset(config)
    return _load_frozen_dataset(config)


def load_off_target_datasets(
    configs: tuple[OffTargetConfig, ...],
) -> tuple[OffTargetDataset, ...]:
    """Load datasets in configuration order."""
    if not isinstance(configs, tuple):
        raise ValueError("Off-target configurations must be a tuple.")
    return tuple(load_off_target_dataset(config) for config in configs)
