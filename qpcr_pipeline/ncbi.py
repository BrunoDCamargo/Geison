"""Frozen NCBI dataset artifacts and validation."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path

from Bio import SeqIO


MANIFEST_NAME = "dataset_manifest.json"
RECORDS_NAME = "records.gb"
SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class AcquiredNcbiDataset:
    records_path: Path
    manifest_path: Path


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read_manifest(path: Path) -> dict[str, object]:
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("NCBI dataset manifest is not valid JSON.") from error
    if not isinstance(raw, dict):
        raise ValueError("NCBI dataset manifest root must be an object.")
    return raw


def _write_json_atomic(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    temporary.replace(path)


def validate_frozen_dataset(dataset_dir: str | Path) -> AcquiredNcbiDataset:
    """Validate and return a complete, immutable local NCBI dataset."""
    directory = Path(dataset_dir)
    manifest_path = directory / MANIFEST_NAME
    records_path = directory / RECORDS_NAME
    if not manifest_path.is_file():
        raise ValueError(f"NCBI dataset manifest '{MANIFEST_NAME}' is missing.")

    manifest = _read_manifest(manifest_path)
    _validate_manifest_header(manifest)
    consolidated = _consolidated_metadata(manifest)
    _validate_consolidated_filename(consolidated)

    if not records_path.is_file():
        raise ValueError(f"NCBI dataset records file '{RECORDS_NAME}' is missing.")
    record_bytes = records_path.read_bytes()
    _validate_byte_size(consolidated, len(record_bytes))
    _validate_sha256(consolidated, record_bytes)

    record_ids = _parse_record_ids(records_path)
    _validate_record_count(consolidated, len(record_ids))
    _validate_resolved_entry_order(manifest, record_ids)

    return AcquiredNcbiDataset(records_path=records_path, manifest_path=manifest_path)


def _validate_manifest_header(manifest: dict[str, object]) -> None:
    schema_version = manifest.get("schema_version")
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != SCHEMA_VERSION
    ):
        raise ValueError(
            "NCBI dataset manifest field 'schema_version' is unsupported: "
            f"expected {SCHEMA_VERSION}."
        )
    if manifest.get("status") != "COMPLETE":
        raise ValueError("NCBI dataset manifest field 'status' must be 'COMPLETE'.")


def _consolidated_metadata(manifest: dict[str, object]) -> dict[str, object]:
    consolidated = manifest.get("consolidated")
    if not isinstance(consolidated, dict):
        raise ValueError("NCBI dataset manifest field 'consolidated' must be an object.")
    return consolidated


def _validate_consolidated_filename(consolidated: dict[str, object]) -> None:
    filename = consolidated.get("filename")
    if filename != RECORDS_NAME:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.filename' must be "
            f"'{RECORDS_NAME}'."
        )


def _validate_byte_size(consolidated: dict[str, object], actual_size: int) -> None:
    expected_size = consolidated.get("byte_size")
    if (
        isinstance(expected_size, bool)
        or not isinstance(expected_size, int)
        or expected_size < 0
    ):
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.byte_size' must be a "
            "non-negative integer."
        )
    if expected_size != actual_size:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.byte_size' does not match "
            f"'{RECORDS_NAME}'."
        )


def _validate_sha256(consolidated: dict[str, object], record_bytes: bytes) -> None:
    expected_checksum = consolidated.get("sha256")
    if (
        not isinstance(expected_checksum, str)
        or len(expected_checksum) != 64
        or any(character not in "0123456789abcdef" for character in expected_checksum)
    ):
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.sha256' must be a "
            "lowercase SHA-256 digest."
        )
    if expected_checksum != _sha256_bytes(record_bytes):
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.sha256' does not match "
            f"'{RECORDS_NAME}'."
        )


def _parse_record_ids(records_path: Path) -> tuple[str, ...]:
    try:
        with records_path.open(encoding="utf-8") as handle:
            return tuple(record.id for record in SeqIO.parse(handle, "genbank"))
    except (OSError, UnicodeError, ValueError) as error:
        raise ValueError(
            f"NCBI dataset records file '{RECORDS_NAME}' is not valid GenBank."
        ) from error


def _validate_record_count(consolidated: dict[str, object], actual_count: int) -> None:
    expected_count = consolidated.get("record_count")
    if (
        isinstance(expected_count, bool)
        or not isinstance(expected_count, int)
        or expected_count < 0
    ):
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.record_count' must be a "
            "non-negative integer."
        )
    if expected_count != actual_count:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.record_count' does not match "
            f"'{RECORDS_NAME}'."
        )


def _validate_resolved_entry_order(
    manifest: dict[str, object], record_ids: tuple[str, ...]
) -> None:
    entries = manifest.get("resolved_entries")
    if not isinstance(entries, list):
        raise ValueError("NCBI dataset manifest field 'resolved_entries' must be a list.")

    accession_versions: list[str] = []
    for index, entry in enumerate(entries):
        if not isinstance(entry, dict):
            raise ValueError(
                "NCBI dataset manifest field "
                f"'resolved_entries[{index}]' must be an object."
            )
        accession_version = entry.get("accession_version")
        if not isinstance(accession_version, str) or not accession_version:
            raise ValueError(
                "NCBI dataset manifest field "
                f"'resolved_entries[{index}].accession_version' must be a "
                "non-empty string."
            )
        accession_versions.append(accession_version)

    if tuple(accession_versions) != record_ids:
        raise ValueError(
            "NCBI dataset manifest field 'resolved_entries[*].accession_version' "
            "does not match the records.gb order."
        )
