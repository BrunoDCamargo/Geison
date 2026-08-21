"""Frozen NCBI dataset artifacts and acquisition."""

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time
from typing import Callable, Literal, Protocol, TypeVar

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig


MANIFEST_NAME = "dataset_manifest.json"
RECORDS_NAME = "records.gb"
SCHEMA_VERSION = 1
_BATCH_DIRECTORY_NAME = "batches"
_BATCH_FILENAME_FORMAT = "batch-{index:05d}.gb"


class NcbiTransientError(RuntimeError):
    """A sanitized NCBI failure that may be retried."""


@dataclass(frozen=True, slots=True)
class NcbiFetchedRecord:
    request_id: str
    record: SeqRecord


@dataclass(frozen=True, slots=True)
class ResolvedNcbiQuery:
    uids: tuple[str, ...]
    reported_count: int
    query_translation: str


class NcbiClient(Protocol):
    def resolve_query(
        self, query: str, max_records: int | None
    ) -> ResolvedNcbiQuery: ...

    def fetch_records(
        self,
        identifiers: tuple[str, ...],
        *,
        identifier_kind: Literal["uid", "accession"],
    ) -> tuple[NcbiFetchedRecord, ...]: ...


@dataclass(frozen=True, slots=True)
class AcquiredNcbiDataset:
    records_path: Path
    manifest_path: Path


_Result = TypeVar("_Result")


def utc_now() -> str:
    """Return the current timestamp in an unambiguous UTC representation."""
    return datetime.now(timezone.utc).isoformat()


def _with_retries(
    operation: Callable[[], _Result],
    retries: int,
    sleep: Callable[[float], object],
) -> _Result:
    for retry_index in range(retries + 1):
        try:
            return operation()
        except NcbiTransientError:
            if retry_index == retries:
                raise
            sleep(min(2**retry_index, 16))
    raise AssertionError("retry loop did not return or raise")


def acquire_ncbi_dataset(
    config: NcbiInputConfig,
    dataset_dir: str | Path,
    *,
    client: NcbiClient | None = None,
    sleep: Callable[[float], object] = time.sleep,
    clock: Callable[[], str] = utc_now,
) -> AcquiredNcbiDataset:
    """Acquire or resume a frozen accession dataset using the supplied client.

    Query acquisition and the production Entrez adapter are deliberately deferred
    to the next task.  Requiring the client here keeps this offline implementation
    from making unconfigured network requests.
    """
    if not config.accessions:
        raise ValueError("Task 3 acquisition requires explicit NCBI accessions.")
    if client is None:
        raise ValueError(
            "An NcbiClient is required for acquisition until the production adapter exists."
        )

    directory = Path(dataset_dir)
    directory.mkdir(parents=True, exist_ok=True)
    batches_directory = directory / _BATCH_DIRECTORY_NAME
    batches_directory.mkdir(exist_ok=True)
    manifest_path = directory / MANIFEST_NAME
    manifest = _load_or_create_accession_manifest(
        manifest_path, config=config, clock=clock
    )
    batches = _planned_batches(config.accessions, config.batch_size)
    valid_batches, cache_changed = _validated_cached_batches(
        manifest, directory, batches
    )

    if cache_changed:
        _set_partial_manifest(manifest, batches, valid_batches, clock)
        _write_json_atomic(manifest_path, manifest)

    if (
        len(valid_batches) == len(batches)
        and manifest["status"] == "COMPLETE"
        and not cache_changed
    ):
        try:
            return validate_frozen_dataset(directory)
        except ValueError:
            _set_partial_manifest(manifest, batches, valid_batches, clock)
            _write_json_atomic(manifest_path, manifest)

    for batch_index, requested_identifiers in enumerate(batches):
        if batch_index in valid_batches:
            continue
        fetched = _with_retries(
            lambda: client.fetch_records(
                requested_identifiers, identifier_kind="accession"
            ),
            config.retries,
            sleep,
        )
        records = _ordered_client_records(fetched, requested_identifiers)
        metadata = _write_batch_atomically(
            batches_directory,
            batch_index=batch_index,
            requested_identifiers=requested_identifiers,
            records=records,
        )
        valid_batches[batch_index] = metadata
        _set_partial_manifest(manifest, batches, valid_batches, clock)
        _write_json_atomic(manifest_path, manifest)

    _write_consolidated_dataset(directory, batches, valid_batches)
    records_path = directory / RECORDS_NAME
    record_bytes = records_path.read_bytes()
    manifest["consolidated"] = {
        "filename": RECORDS_NAME,
        "record_count": len(config.accessions),
        "byte_size": len(record_bytes),
        "sha256": _sha256_bytes(record_bytes),
    }
    manifest["status"] = "COMPLETE"
    manifest["updated_at"] = clock()
    _write_json_atomic(manifest_path, manifest)
    return validate_frozen_dataset(directory)


def _planned_batches(
    accessions: tuple[str, ...], batch_size: int
) -> tuple[tuple[str, ...], ...]:
    return tuple(
        accessions[index : index + batch_size]
        for index in range(0, len(accessions), batch_size)
    )


def _initial_accession_manifest(
    config: NcbiInputConfig, clock: Callable[[], str]
) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PARTIAL",
        "source": {
            "mode": "accessions",
            "database": "nuccore",
            "requested_accessions": list(config.accessions),
        },
        "batch_size": config.batch_size,
        "retries": config.retries,
        "resolved_entries": [],
        "completed_batches": [],
        "consolidated": None,
        "created_at": clock(),
        "updated_at": clock(),
        "tool": "geison-qpcr",
    }


def _load_or_create_accession_manifest(
    manifest_path: Path, *, config: NcbiInputConfig, clock: Callable[[], str]
) -> dict[str, object]:
    if not manifest_path.exists():
        manifest = _initial_accession_manifest(config, clock)
        _write_json_atomic(manifest_path, manifest)
        return manifest

    manifest = _read_manifest(manifest_path)
    expected_source = {
        "mode": "accessions",
        "database": "nuccore",
        "requested_accessions": list(config.accessions),
    }
    if manifest.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("NCBI dataset manifest schema does not support acquisition resume.")
    if manifest.get("source") != expected_source:
        raise ValueError(
            "NCBI dataset manifest source does not match the requested accessions."
        )
    if manifest.get("batch_size") != config.batch_size:
        raise ValueError("NCBI dataset manifest batch_size does not match this configuration.")
    if manifest.get("retries") != config.retries:
        raise ValueError("NCBI dataset manifest retries does not match this configuration.")
    if manifest.get("status") not in {"PARTIAL", "COMPLETE"}:
        raise ValueError("NCBI dataset manifest status is not resumable.")
    if not isinstance(manifest.get("completed_batches"), list):
        raise ValueError("NCBI dataset manifest completed_batches must be a list.")
    if not isinstance(manifest.get("resolved_entries"), list):
        raise ValueError("NCBI dataset manifest resolved_entries must be a list.")
    if not isinstance(manifest.get("created_at"), str):
        raise ValueError("NCBI dataset manifest created_at must be a string.")
    return manifest


def _validated_cached_batches(
    manifest: dict[str, object],
    directory: Path,
    batches: tuple[tuple[str, ...], ...],
) -> tuple[dict[int, dict[str, object]], bool]:
    completed = manifest["completed_batches"]
    assert isinstance(completed, list)
    metadata_by_filename: dict[str, list[dict[str, object]]] = {}
    for metadata in completed:
        if isinstance(metadata, dict) and isinstance(metadata.get("filename"), str):
            metadata_by_filename.setdefault(metadata["filename"], []).append(metadata)

    valid: dict[int, dict[str, object]] = {}
    for batch_index, requested_identifiers in enumerate(batches):
        filename = _batch_filename(batch_index)
        candidates = metadata_by_filename.get(filename, [])
        if len(candidates) != 1:
            continue
        metadata = candidates[0]
        if _cached_batch_is_valid(directory, metadata, requested_identifiers):
            valid[batch_index] = metadata

    expected_metadata = [valid[index] for index in range(len(batches)) if index in valid]
    resolved_entries = _resolved_entries_from_batches(batches, valid)
    changed = (
        completed != expected_metadata
        or manifest.get("resolved_entries") != resolved_entries
        or manifest.get("status") == "COMPLETE" and len(valid) != len(batches)
    )
    return valid, changed


def _cached_batch_is_valid(
    directory: Path,
    metadata: dict[str, object],
    requested_identifiers: tuple[str, ...],
) -> bool:
    filename = metadata.get("filename")
    requested = metadata.get("requested_identifiers")
    record_count = metadata.get("record_count")
    byte_size = metadata.get("byte_size")
    checksum = metadata.get("sha256")
    record_ids = metadata.get("record_ids")
    resolved_entries = metadata.get("resolved_entries")
    if (
        not isinstance(filename, str)
        or not _is_exact_string_list(requested, requested_identifiers)
        or isinstance(record_count, bool)
        or not isinstance(record_count, int)
        or record_count != len(requested_identifiers)
        or isinstance(byte_size, bool)
        or not isinstance(byte_size, int)
        or byte_size < 0
        or not isinstance(checksum, str)
        or not isinstance(resolved_entries, list)
    ):
        return False

    path = directory / _BATCH_DIRECTORY_NAME / filename
    try:
        contents = path.read_bytes()
        records = _parse_genbank_records(path)
    except (OSError, UnicodeError, ValueError):
        return False
    if byte_size != len(contents) or checksum != _sha256_bytes(contents):
        return False
    parsed_ids = tuple(record.id for record in records)
    if not _is_exact_string_list(record_ids, parsed_ids):
        return False
    if len(records) != record_count:
        return False
    return resolved_entries == _resolved_entries(requested_identifiers, records)


def _is_exact_string_list(value: object, expected: tuple[str, ...]) -> bool:
    return (
        isinstance(value, list)
        and all(isinstance(item, str) for item in value)
        and tuple(value) == expected
    )


def _ordered_client_records(
    fetched: object, requested_identifiers: tuple[str, ...]
) -> tuple[SeqRecord, ...]:
    if not isinstance(fetched, tuple) or len(fetched) != len(requested_identifiers):
        raise ValueError(
            "NCBI client response must contain exactly one record for every requested identifier."
        )
    by_request_id: dict[str, SeqRecord] = {}
    for fetched_record in fetched:
        if (
            not isinstance(fetched_record, NcbiFetchedRecord)
            or not isinstance(fetched_record.request_id, str)
            or not isinstance(fetched_record.record, SeqRecord)
            or fetched_record.request_id in by_request_id
        ):
            raise ValueError(
                "NCBI client response must contain exactly one record for every requested identifier."
            )
        by_request_id[fetched_record.request_id] = fetched_record.record
    if set(by_request_id) != set(requested_identifiers):
        raise ValueError(
            "NCBI client response must contain exactly one record for every requested identifier."
        )
    return tuple(by_request_id[identifier] for identifier in requested_identifiers)


def _write_batch_atomically(
    batches_directory: Path,
    *,
    batch_index: int,
    requested_identifiers: tuple[str, ...],
    records: tuple[SeqRecord, ...],
) -> dict[str, object]:
    filename = _batch_filename(batch_index)
    path = batches_directory / filename
    temporary = path.with_suffix(path.suffix + ".tmp")
    try:
        _write_genbank(temporary, records)
        parsed_records = _parse_genbank_records(temporary)
        if tuple(record.id for record in parsed_records) != tuple(record.id for record in records):
            raise ValueError("NCBI batch records changed while writing GenBank output.")
        contents = temporary.read_bytes()
        temporary.replace(path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
    return _batch_metadata(filename, requested_identifiers, parsed_records, contents)


def _batch_filename(batch_index: int) -> str:
    return _BATCH_FILENAME_FORMAT.format(index=batch_index)


def _batch_metadata(
    filename: str,
    requested_identifiers: tuple[str, ...],
    records: tuple[SeqRecord, ...],
    contents: bytes,
) -> dict[str, object]:
    return {
        "filename": filename,
        "requested_identifiers": list(requested_identifiers),
        "record_count": len(records),
        "byte_size": len(contents),
        "sha256": _sha256_bytes(contents),
        "record_ids": [record.id for record in records],
        "resolved_entries": _resolved_entries(requested_identifiers, records),
    }


def _resolved_entries(
    requested_identifiers: tuple[str, ...], records: tuple[SeqRecord, ...]
) -> list[dict[str, object]]:
    return [
        {
            "requested_accession": requested_accession,
            "uid": None,
            "accession": re.sub(r"\.\d+$", "", record.id),
            "accession_version": record.id,
        }
        for requested_accession, record in zip(requested_identifiers, records, strict=True)
    ]


def _resolved_entries_from_batches(
    batches: tuple[tuple[str, ...], ...], valid_batches: dict[int, dict[str, object]]
) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    for batch_index in range(len(batches)):
        metadata = valid_batches.get(batch_index)
        if metadata is not None:
            batch_entries = metadata["resolved_entries"]
            assert isinstance(batch_entries, list)
            entries.extend(batch_entries)
    return entries


def _set_partial_manifest(
    manifest: dict[str, object],
    batches: tuple[tuple[str, ...], ...],
    valid_batches: dict[int, dict[str, object]],
    clock: Callable[[], str],
) -> None:
    manifest["status"] = "PARTIAL"
    manifest["completed_batches"] = [
        valid_batches[index] for index in range(len(batches)) if index in valid_batches
    ]
    manifest["resolved_entries"] = _resolved_entries_from_batches(batches, valid_batches)
    manifest["consolidated"] = None
    manifest["updated_at"] = clock()


def _write_consolidated_dataset(
    directory: Path,
    batches: tuple[tuple[str, ...], ...],
    valid_batches: dict[int, dict[str, object]],
) -> None:
    records: list[SeqRecord] = []
    for batch_index in range(len(batches)):
        metadata = valid_batches[batch_index]
        path = directory / _BATCH_DIRECTORY_NAME / metadata["filename"]
        assert isinstance(path, Path)
        parsed_records = _parse_genbank_records(path)
        records.extend(parsed_records)
    records_tuple = tuple(records)
    expected_entries = _resolved_entries_from_batches(batches, valid_batches)
    if _resolved_entries(
        tuple(entry["requested_accession"] for entry in expected_entries), records_tuple
    ) != expected_entries:
        raise ValueError("Validated NCBI batch metadata no longer matches its records.")

    records_path = directory / RECORDS_NAME
    temporary = records_path.with_suffix(records_path.suffix + ".tmp")
    try:
        _write_genbank(temporary, records_tuple)
        parsed_records = _parse_genbank_records(temporary)
        if tuple(record.id for record in parsed_records) != tuple(record.id for record in records_tuple):
            raise ValueError("NCBI consolidated records changed while writing GenBank output.")
        temporary.replace(records_path)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _write_genbank(path: Path, records: tuple[SeqRecord, ...]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        SeqIO.write(records, handle, "genbank")


def _parse_genbank_records(path: Path) -> tuple[SeqRecord, ...]:
    try:
        with path.open(encoding="utf-8") as handle:
            return tuple(SeqIO.parse(handle, "genbank"))
    except (OSError, UnicodeError, ValueError) as error:
        raise ValueError(f"NCBI GenBank artifact '{path.name}' is not valid GenBank.") from error


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
