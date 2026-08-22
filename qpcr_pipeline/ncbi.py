"""Frozen NCBI dataset artifacts and acquisition."""

from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import threading
import time
from typing import Callable, Literal, Mapping, Protocol, TypeVar
from urllib.error import HTTPError, URLError
from weakref import WeakKeyDictionary

from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig, validate_ncbi_input_config


MANIFEST_NAME = "dataset_manifest.json"
RECORDS_NAME = "records.gb"
SCHEMA_VERSION = 1
_BATCH_DIRECTORY_NAME = "batches"
_BATCH_FILENAME_FORMAT = "batch-{index:05d}.gb"
_ENTREZ_LOCKS_GUARD = threading.Lock()
_ENTREZ_REQUEST_LOCKS = WeakKeyDictionary()
_ENTREZ_FALLBACK_LOCK = threading.Lock()
_MISSING_ATTRIBUTE = object()
_MANIFEST_FIELDS = frozenset(
    {
        "schema_version",
        "status",
        "source",
        "batch_size",
        "retries",
        "resolved_entries",
        "completed_batches",
        "consolidated",
        "created_at",
        "updated_at",
        "tool",
    }
)
_ACCESSION_SOURCE_FIELDS = frozenset(
    {"mode", "database", "requested_accessions"}
)
_QUERY_SOURCE_FIELDS = frozenset(
    {
        "mode",
        "database",
        "query",
        "resolved_uids",
        "reported_count",
        "selected_count",
        "max_records",
        "query_translation",
    }
)
_RESOLVED_ENTRY_FIELDS = frozenset(
    {"requested_accession", "uid", "accession", "accession_version"}
)
_COMPLETED_BATCH_FIELDS = frozenset(
    {
        "filename",
        "requested_identifiers",
        "record_count",
        "byte_size",
        "sha256",
        "record_ids",
        "resolved_entries",
    }
)
_CONSOLIDATED_FIELDS = frozenset(
    {"filename", "record_count", "byte_size", "sha256"}
)


class NcbiTransientError(RuntimeError):
    """A sanitized NCBI failure that may be retried."""


class NcbiRequestError(RuntimeError):
    """A sanitized NCBI failure that must not be retried."""


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
class BioEntrezClient:
    """Credential-safe Bio.Entrez implementation of the NCBI client boundary."""

    email: str = field(repr=False)
    api_key: str | None = field(default=None, repr=False)
    entrez_module: object = field(default=Entrez, repr=False)
    _request_lock: object = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        email = self.email.strip()
        if not email:
            raise ValueError("NCBI_EMAIL must be set for live NCBI acquisition.")
        if self.api_key is not None and not isinstance(self.api_key, str):
            raise ValueError("NCBI_API_KEY must be a string when set.")
        object.__setattr__(self, "email", email)
        api_key = self.api_key.strip() if self.api_key else ""
        object.__setattr__(self, "api_key", api_key or None)
        object.__setattr__(self, "_request_lock", _entrez_request_lock(self.entrez_module))

    @classmethod
    def from_environment(
        cls,
        environ: Mapping[str, str] = os.environ,
        *,
        entrez_module: object = Entrez,
    ) -> "BioEntrezClient":
        email = environ.get("NCBI_EMAIL", "").strip()
        if not email:
            raise ValueError("NCBI_EMAIL must be set for live NCBI acquisition.")
        return cls(
            email=email,
            api_key=environ.get("NCBI_API_KEY") or None,
            entrez_module=entrez_module,
        )

    def resolve_query(
        self, query: str, max_records: int | None
    ) -> ResolvedNcbiQuery:
        requested_count = min(max_records, 10_000) if max_records is not None else 10_000
        uids: list[str] = []
        reported_count: int | None = None
        query_translation: str | None = None

        while reported_count is None or len(uids) < min(
            reported_count, max_records if max_records is not None else reported_count
        ):
            retmax = requested_count if reported_count is None else min(
                10_000,
                min(
                    reported_count,
                    max_records if max_records is not None else reported_count,
                )
                - len(uids),
            )
            response = self._search_page(query, retstart=len(uids), retmax=retmax)
            page_count, page_uids, page_translation = _validated_search_response(response)
            if reported_count is None:
                reported_count = page_count
                query_translation = page_translation
            elif page_count != reported_count or page_translation != query_translation:
                raise ValueError("NCBI query response composition changed between pages.")
            uids.extend(page_uids)
            if not page_uids:
                break

        assert reported_count is not None
        assert query_translation is not None
        selected_count = min(
            reported_count, max_records if max_records is not None else reported_count
        )
        if len(uids) != selected_count or len(set(uids)) != len(uids):
            raise ValueError("NCBI query response count does not match its UID composition.")
        return ResolvedNcbiQuery(
            uids=tuple(uids),
            reported_count=reported_count,
            query_translation=query_translation,
        )

    def fetch_records(
        self,
        identifiers: tuple[str, ...],
        *,
        identifier_kind: Literal["uid", "accession"],
    ) -> tuple[NcbiFetchedRecord, ...]:
        if identifier_kind not in {"uid", "accession"}:
            raise ValueError("NCBI identifier kind must be 'uid' or 'accession'.")
        if (
            not isinstance(identifiers, tuple)
            or not identifiers
            or any(not isinstance(identifier, str) or not identifier for identifier in identifiers)
            or len(set(identifiers)) != len(identifiers)
        ):
            raise ValueError("NCBI fetch identifiers must be unique non-empty strings.")

        records = self._fetch_genbank(identifiers)
        if len(records) != len(identifiers):
            raise ValueError(
                "NCBI fetch response must contain exactly one record for every identifier."
            )
        if identifier_kind == "uid":
            return tuple(
                NcbiFetchedRecord(request_id=identifier, record=record)
                for identifier, record in zip(identifiers, records, strict=True)
            )
        return _map_accession_records(identifiers, records)

    def _configure(self) -> None:
        setattr(self.entrez_module, "email", self.email)
        setattr(self.entrez_module, "tool", "geison-qpcr")
        setattr(self.entrez_module, "api_key", self.api_key)

    def _create_request_handle(self, operation: Callable[[], object]) -> object:
        with self._request_lock:
            self._configure()
            previous_max_tries = getattr(
                self.entrez_module, "max_tries", _MISSING_ATTRIBUTE
            )
            previous_sleep = getattr(
                self.entrez_module, "sleep_between_tries", _MISSING_ATTRIBUTE
            )
            try:
                setattr(self.entrez_module, "max_tries", 1)
                setattr(self.entrez_module, "sleep_between_tries", 0)
                return operation()
            finally:
                if previous_max_tries is _MISSING_ATTRIBUTE:
                    delattr(self.entrez_module, "max_tries")
                else:
                    setattr(self.entrez_module, "max_tries", previous_max_tries)
                if previous_sleep is _MISSING_ATTRIBUTE:
                    delattr(self.entrez_module, "sleep_between_tries")
                else:
                    setattr(self.entrez_module, "sleep_between_tries", previous_sleep)

    def _search_page(self, query: str, *, retstart: int, retmax: int) -> object:
        def acquire_handle() -> object:
            return self._create_request_handle(
                lambda: self.entrez_module.esearch(
                    db="nuccore",
                    term=query,
                    retstart=retstart,
                    retmax=retmax,
                )
            )

        handle = _run_entrez_request(acquire_handle)
        try:
            return self.entrez_module.read(handle)
        finally:
            handle.close()

    def _fetch_genbank(self, identifiers: tuple[str, ...]) -> tuple[SeqRecord, ...]:
        def acquire_handle() -> object:
            return self._create_request_handle(
                lambda: self.entrez_module.efetch(
                    db="nuccore",
                    id=",".join(identifiers),
                    rettype="gb",
                    retmode="text",
                )
            )

        handle = _run_entrez_request(acquire_handle)
        try:
            return tuple(SeqIO.parse(handle, "genbank"))
        finally:
            handle.close()


@dataclass(frozen=True, slots=True)
class AcquiredNcbiDataset:
    records_path: Path
    manifest_path: Path


_Result = TypeVar("_Result")


def _entrez_request_lock(entrez_module: object) -> object:
    with _ENTREZ_LOCKS_GUARD:
        try:
            request_lock = _ENTREZ_REQUEST_LOCKS.get(entrez_module)
            if request_lock is None:
                request_lock = threading.Lock()
                _ENTREZ_REQUEST_LOCKS[entrez_module] = request_lock
            return request_lock
        except TypeError:
            return _ENTREZ_FALLBACK_LOCK


def _run_entrez_request(operation: Callable[[], _Result]) -> _Result:
    translated: NcbiTransientError | NcbiRequestError
    try:
        return operation()
    except HTTPError as error:
        error.close()
        if error.code == 429 or 500 <= error.code <= 599:
            translated = NcbiTransientError("NCBI request failed temporarily.")
        else:
            translated = NcbiRequestError(
                f"NCBI request failed with HTTP status {error.code}."
            )
    except (URLError, TimeoutError, OSError):
        translated = NcbiTransientError("NCBI network request failed temporarily.")
    raise translated from None


def _validated_search_response(response: object) -> tuple[int, tuple[str, ...], str]:
    if not isinstance(response, Mapping):
        raise ValueError("NCBI query response must be a mapping.")
    raw_count = response.get("Count")
    raw_uids = response.get("IdList")
    translation = response.get("QueryTranslation")
    try:
        reported_count = int(raw_count)
    except (TypeError, ValueError):
        raise ValueError("NCBI query response count must be a non-negative integer.") from None
    if reported_count < 0 or str(reported_count) != str(raw_count):
        raise ValueError("NCBI query response count must be a non-negative integer.")
    if (
        not isinstance(raw_uids, list)
        or any(not isinstance(uid, str) or not uid for uid in raw_uids)
        or not isinstance(translation, str)
    ):
        raise ValueError("NCBI query response has an invalid UID composition.")
    return reported_count, tuple(raw_uids), translation


def _map_accession_records(
    identifiers: tuple[str, ...], records: tuple[SeqRecord, ...]
) -> tuple[NcbiFetchedRecord, ...]:
    records_by_request: dict[str, SeqRecord] = {}
    requested = set(identifiers)
    returned_accessions: set[str] = set()
    for record in records:
        if record.id in returned_accessions:
            raise ValueError(
                "NCBI fetch response accession composition does not match the request."
            )
        returned_accessions.add(record.id)
        base_accession = re.sub(r"\.\d+$", "", record.id)
        if record.id in requested and record.id not in records_by_request:
            request_id = record.id
        elif base_accession in requested and base_accession not in records_by_request:
            request_id = base_accession
        else:
            raise ValueError(
                "NCBI fetch response accession composition does not match the request."
            )
        records_by_request[request_id] = record
    if set(records_by_request) != requested:
        raise ValueError(
            "NCBI fetch response accession composition does not match the request."
        )
    return tuple(
        NcbiFetchedRecord(request_id=identifier, record=records_by_request[identifier])
        for identifier in identifiers
    )


def utc_now() -> str:
    """Return the current timestamp in an unambiguous UTC representation."""
    return datetime.now(timezone.utc).isoformat()


def _with_retries(
    operation: Callable[[], _Result],
    retries: int,
    sleep: Callable[[float], object],
    operation_label: str,
) -> _Result:
    for retry_index in range(retries + 1):
        try:
            return operation()
        except NcbiTransientError:
            if retry_index == retries:
                break
            sleep(min(2**retry_index, 16))
    raise NcbiTransientError(
        f"{operation_label} failed after {retries + 1} attempts."
    ) from None


def acquire_ncbi_dataset(
    config: NcbiInputConfig,
    dataset_dir: str | Path,
    *,
    client: NcbiClient | None = None,
    sleep: Callable[[float], object] = time.sleep,
    clock: Callable[[], str] = utc_now,
) -> AcquiredNcbiDataset:
    """Acquire or resume a frozen query or accession dataset."""
    source_mode = _validate_acquisition_config(config)
    if client is None:
        client = BioEntrezClient.from_environment()

    directory = Path(dataset_dir)
    manifest_path = directory / MANIFEST_NAME
    if not manifest_path.exists():
        directory.mkdir(parents=True, exist_ok=True)
    manifest, identifiers = _load_or_create_manifest(
        manifest_path,
        config=config,
        client=client,
        sleep=sleep,
        clock=clock,
    )
    batches_directory = directory / _BATCH_DIRECTORY_NAME
    batches_directory.mkdir(exist_ok=True)
    identifier_kind: Literal["uid", "accession"] = (
        "uid" if source_mode == "query" else "accession"
    )
    batches = _planned_batches(identifiers, config.batch_size)
    valid_batches, cache_changed = _validated_cached_batches(
        manifest, directory, batches, identifier_kind
    )

    if cache_changed:
        _set_partial_manifest(
            manifest, batches, valid_batches, clock
        )
        _write_json_atomic(manifest_path, manifest)

    if (
        len(valid_batches) == len(batches)
        and manifest["status"] == "COMPLETE"
        and not cache_changed
    ):
        try:
            return validate_frozen_dataset(directory)
        except ValueError:
            _set_partial_manifest(
                manifest, batches, valid_batches, clock
            )
            _write_json_atomic(manifest_path, manifest)

    for batch_index, requested_identifiers in enumerate(batches):
        if batch_index in valid_batches:
            continue
        fetched = _with_retries(
            lambda: client.fetch_records(
                requested_identifiers, identifier_kind=identifier_kind
            ),
            config.retries,
            sleep,
            "NCBI record fetch",
        )
        records = _ordered_client_records(fetched, requested_identifiers)
        metadata = _write_batch_atomically(
            batches_directory,
            batch_index=batch_index,
            requested_identifiers=requested_identifiers,
            records=records,
            identifier_kind=identifier_kind,
        )
        valid_batches[batch_index] = metadata
        _set_partial_manifest(
            manifest, batches, valid_batches, clock
        )
        _write_json_atomic(manifest_path, manifest)

    _write_consolidated_dataset(
        directory, batches, valid_batches, identifier_kind
    )
    records_path = directory / RECORDS_NAME
    record_bytes = records_path.read_bytes()
    manifest["consolidated"] = {
        "filename": RECORDS_NAME,
        "record_count": len(identifiers),
        "byte_size": len(record_bytes),
        "sha256": _sha256_bytes(record_bytes),
    }
    manifest["status"] = "COMPLETE"
    manifest["updated_at"] = clock()
    _write_json_atomic(manifest_path, manifest)
    return validate_frozen_dataset(directory)


def _validate_acquisition_config(
    config: NcbiInputConfig,
) -> Literal["query", "accessions"]:
    source_mode = validate_ncbi_input_config(config)
    if source_mode == "frozen":
        raise ValueError("NCBI acquisition cannot acquire a frozen_dataset source.")
    return source_mode


def _planned_batches(
    identifiers: tuple[str, ...], batch_size: int
) -> tuple[tuple[str, ...], ...]:
    return tuple(
        identifiers[index : index + batch_size]
        for index in range(0, len(identifiers), batch_size)
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


def _initial_query_manifest(
    config: NcbiInputConfig,
    resolution: ResolvedNcbiQuery,
    clock: Callable[[], str],
) -> dict[str, object]:
    return {
        "schema_version": SCHEMA_VERSION,
        "status": "PARTIAL",
        "source": {
            "mode": "query",
            "database": "nuccore",
            "query": config.query,
            "resolved_uids": list(resolution.uids),
            "reported_count": resolution.reported_count,
            "selected_count": len(resolution.uids),
            "max_records": config.max_records,
            "query_translation": resolution.query_translation,
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


def _load_or_create_manifest(
    manifest_path: Path,
    *,
    config: NcbiInputConfig,
    client: NcbiClient,
    sleep: Callable[[float], object],
    clock: Callable[[], str],
) -> tuple[dict[str, object], tuple[str, ...]]:
    if manifest_path.exists():
        manifest = _read_manifest(manifest_path)
        _validate_resumable_manifest(manifest, config)
        source = manifest["source"]
        assert isinstance(source, dict)
        if source["mode"] == "query":
            resolved_uids = source["resolved_uids"]
            assert isinstance(resolved_uids, list)
            return manifest, tuple(resolved_uids)
        return manifest, config.accessions

    if config.query is None:
        manifest = _initial_accession_manifest(config, clock)
        identifiers = config.accessions
    else:
        resolution = _with_retries(
            lambda: client.resolve_query(config.query, config.max_records),
            config.retries,
            sleep,
            "NCBI query resolution",
        )
        _validate_query_resolution(resolution, config.max_records)
        manifest = _initial_query_manifest(config, resolution, clock)
        identifiers = resolution.uids
    _write_json_atomic(manifest_path, manifest)
    return manifest, identifiers


def _validate_query_resolution(
    resolution: object, max_records: int | None
) -> None:
    if not isinstance(resolution, ResolvedNcbiQuery):
        raise ValueError("NCBI query resolution has an invalid composition.")
    if (
        isinstance(resolution.reported_count, bool)
        or not isinstance(resolution.reported_count, int)
        or resolution.reported_count < 0
    ):
        raise ValueError("NCBI query resolution count or UID composition is invalid.")
    selected_count = min(
        resolution.reported_count,
        max_records if max_records is not None else resolution.reported_count,
    )
    if (
        not isinstance(resolution.uids, tuple)
        or any(not isinstance(uid, str) or not uid for uid in resolution.uids)
        or len(set(resolution.uids)) != len(resolution.uids)
        or len(resolution.uids) != selected_count
        or not isinstance(resolution.query_translation, str)
    ):
        raise ValueError("NCBI query resolution count or UID composition is invalid.")


def _validate_resumable_manifest(
    manifest: dict[str, object], config: NcbiInputConfig
) -> None:
    unknown_fields = set(manifest) - _MANIFEST_FIELDS
    missing_fields = _MANIFEST_FIELDS - set(manifest)
    if unknown_fields:
        raise ValueError(
            "NCBI dataset manifest has unrecognized fields: "
            f"{', '.join(sorted(unknown_fields))}."
        )
    if missing_fields:
        raise ValueError(
            "NCBI dataset manifest is missing required fields: "
            f"{', '.join(sorted(missing_fields))}."
        )

    schema_version = manifest["schema_version"]
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != SCHEMA_VERSION
    ):
        raise ValueError("NCBI dataset manifest schema does not support acquisition resume.")
    _validate_resumable_source(manifest["source"], config)
    if (
        isinstance(manifest["batch_size"], bool)
        or not isinstance(manifest["batch_size"], int)
        or manifest["batch_size"] != config.batch_size
    ):
        raise ValueError("NCBI dataset manifest batch_size does not match this configuration.")
    if (
        isinstance(manifest["retries"], bool)
        or not isinstance(manifest["retries"], int)
        or manifest["retries"] != config.retries
    ):
        raise ValueError("NCBI dataset manifest retries does not match this configuration.")
    status = manifest["status"]
    if status not in {"PARTIAL", "COMPLETE"}:
        raise ValueError("NCBI dataset manifest status is not resumable.")
    if manifest["tool"] != "geison-qpcr":
        raise ValueError("NCBI dataset manifest tool must be 'geison-qpcr'.")
    _validate_manifest_timestamp(manifest["created_at"], "created_at")
    _validate_manifest_timestamp(manifest["updated_at"], "updated_at")
    if status == "PARTIAL" and manifest["consolidated"] is not None:
        raise ValueError("NCBI PARTIAL dataset manifest consolidated must be null.")
    if status == "COMPLETE" and not isinstance(manifest["consolidated"], dict):
        raise ValueError("NCBI COMPLETE dataset manifest consolidated must be an object.")
    if not isinstance(manifest["completed_batches"], list):
        raise ValueError("NCBI dataset manifest completed_batches must be a list.")
    if not isinstance(manifest["resolved_entries"], list):
        raise ValueError("NCBI dataset manifest resolved_entries must be a list.")


def _validate_resumable_source(source: object, config: NcbiInputConfig) -> None:
    if config.query is None:
        expected_source = {
            "mode": "accessions",
            "database": "nuccore",
            "requested_accessions": list(config.accessions),
        }
        if source != expected_source:
            raise ValueError(
                "NCBI dataset manifest source does not match the requested accessions."
            )
        return

    expected_fields = {
        "mode",
        "database",
        "query",
        "resolved_uids",
        "reported_count",
        "selected_count",
        "max_records",
        "query_translation",
    }
    if not isinstance(source, dict) or set(source) != expected_fields:
        raise ValueError("NCBI dataset manifest source is not a canonical query source.")
    if (
        source["mode"] != "query"
        or source["database"] != "nuccore"
        or source["query"] != config.query
        or (
            config.max_records is None
            and source["max_records"] is not None
        )
        or (
            config.max_records is not None
            and (
                isinstance(source["max_records"], bool)
                or not isinstance(source["max_records"], int)
                or source["max_records"] != config.max_records
            )
        )
    ):
        raise ValueError("NCBI dataset manifest source does not match this query configuration.")
    if not isinstance(source["resolved_uids"], list):
        raise ValueError(
            "NCBI dataset manifest query source resolved_uids must be a list."
        )
    resolution = ResolvedNcbiQuery(
        uids=tuple(source["resolved_uids"]),
        reported_count=source["reported_count"],
        query_translation=source["query_translation"],
    )
    try:
        _validate_query_resolution(resolution, config.max_records)
    except ValueError:
        raise ValueError("NCBI dataset manifest query source composition is invalid.") from None
    if (
        isinstance(source["selected_count"], bool)
        or not isinstance(source["selected_count"], int)
        or source["selected_count"] != len(resolution.uids)
    ):
        raise ValueError("NCBI dataset manifest query source selected_count is invalid.")


def _validate_manifest_timestamp(value: object, field: str) -> None:
    if not isinstance(value, str):
        raise ValueError(f"NCBI dataset manifest {field} must be a UTC timestamp.")
    try:
        timestamp = datetime.fromisoformat(value)
    except ValueError as error:
        raise ValueError(
            f"NCBI dataset manifest {field} must be a UTC timestamp."
        ) from error
    if timestamp.tzinfo is None or timestamp.utcoffset() != timedelta(0):
        raise ValueError(f"NCBI dataset manifest {field} must be a UTC timestamp.")


def _validated_cached_batches(
    manifest: dict[str, object],
    directory: Path,
    batches: tuple[tuple[str, ...], ...],
    identifier_kind: Literal["uid", "accession"],
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
        if _cached_batch_is_valid(
            directory, metadata, requested_identifiers, identifier_kind
        ):
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
    identifier_kind: Literal["uid", "accession"],
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
    return resolved_entries == _resolved_entries(
        requested_identifiers, records, identifier_kind
    )


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
    identifier_kind: Literal["uid", "accession"],
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
    return _batch_metadata(
        filename,
        requested_identifiers,
        parsed_records,
        contents,
        identifier_kind,
    )


def _batch_filename(batch_index: int) -> str:
    return _BATCH_FILENAME_FORMAT.format(index=batch_index)


def _batch_metadata(
    filename: str,
    requested_identifiers: tuple[str, ...],
    records: tuple[SeqRecord, ...],
    contents: bytes,
    identifier_kind: Literal["uid", "accession"],
) -> dict[str, object]:
    return {
        "filename": filename,
        "requested_identifiers": list(requested_identifiers),
        "record_count": len(records),
        "byte_size": len(contents),
        "sha256": _sha256_bytes(contents),
        "record_ids": [record.id for record in records],
        "resolved_entries": _resolved_entries(
            requested_identifiers, records, identifier_kind
        ),
    }


def _resolved_entries(
    requested_identifiers: tuple[str, ...],
    records: tuple[SeqRecord, ...],
    identifier_kind: Literal["uid", "accession"],
) -> list[dict[str, object]]:
    return [
        {
            "requested_accession": identifier if identifier_kind == "accession" else None,
            "uid": identifier if identifier_kind == "uid" else None,
            "accession": re.sub(r"\.\d+$", "", record.id),
            "accession_version": record.id,
        }
        for identifier, record in zip(requested_identifiers, records, strict=True)
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
    identifier_kind: Literal["uid", "accession"],
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
    requested_identifiers = tuple(
        entry["uid"] if identifier_kind == "uid" else entry["requested_accession"]
        for entry in expected_entries
    )
    if _resolved_entries(
        requested_identifiers, records_tuple, identifier_kind
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
    _require_exact_fields(manifest, _MANIFEST_FIELDS, "manifest")
    _validate_complete_manifest_header(manifest)
    identifier_kind, identifiers = _validate_complete_source(manifest["source"])
    resolved_entries, accession_versions = _validate_complete_resolved_entries(
        manifest["resolved_entries"], identifiers, identifier_kind, "resolved_entries"
    )
    _validate_complete_batches(
        directory,
        manifest["completed_batches"],
        identifiers,
        resolved_entries,
        accession_versions,
        identifier_kind,
        manifest["batch_size"],
    )
    consolidated = _validate_complete_consolidated(
        manifest["consolidated"], len(identifiers)
    )

    if not records_path.is_file():
        raise ValueError(f"NCBI dataset records file '{RECORDS_NAME}' is missing.")
    record_bytes = records_path.read_bytes()
    _validate_artifact_bytes(
        consolidated, record_bytes, "consolidated", RECORDS_NAME
    )
    records = _parse_genbank_records(records_path)
    record_ids = tuple(record.id for record in records)
    if len(records) != consolidated["record_count"]:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.record_count' does not match "
            f"'{RECORDS_NAME}'."
        )
    if record_ids != accession_versions:
        raise ValueError(
            "NCBI dataset manifest field 'resolved_entries[*].accession_version' "
            "does not match the records.gb order."
        )

    return AcquiredNcbiDataset(records_path=records_path, manifest_path=manifest_path)


def _require_exact_fields(
    value: object, expected: frozenset[str], field_path: str
) -> dict[str, object]:
    if not isinstance(value, dict):
        raise ValueError(f"NCBI dataset manifest field '{field_path}' must be an object.")
    unknown = set(value) - expected
    if unknown:
        rendered = ", ".join(
            f"{field_path}.{field}" if field_path != "manifest" else str(field)
            for field in sorted(unknown)
        )
        raise ValueError(
            f"NCBI dataset manifest field(s) '{rendered}' are unrecognized."
        )
    missing = expected - set(value)
    if missing:
        rendered = ", ".join(
            f"{field_path}.{field}" if field_path != "manifest" else str(field)
            for field in sorted(missing)
        )
        raise ValueError(
            f"NCBI dataset manifest field(s) '{rendered}' are required."
        )
    return value


def _validate_complete_manifest_header(manifest: dict[str, object]) -> None:
    schema_version = manifest["schema_version"]
    if (
        isinstance(schema_version, bool)
        or not isinstance(schema_version, int)
        or schema_version != SCHEMA_VERSION
    ):
        raise ValueError(
            "NCBI dataset manifest field 'schema_version' is unsupported: "
            f"expected {SCHEMA_VERSION}."
        )
    if manifest["status"] != "COMPLETE":
        raise ValueError("NCBI dataset manifest field 'status' must be 'COMPLETE'.")
    if manifest["tool"] != "geison-qpcr":
        raise ValueError(
            "NCBI dataset manifest field 'tool' must be 'geison-qpcr'."
        )
    _validate_manifest_timestamp(manifest["created_at"], "created_at")
    _validate_manifest_timestamp(manifest["updated_at"], "updated_at")
    _validate_manifest_integer(
        manifest["batch_size"], "batch_size", minimum=1, maximum=500
    )
    _validate_manifest_integer(
        manifest["retries"], "retries", minimum=0, maximum=10
    )
    if not isinstance(manifest["consolidated"], dict):
        raise ValueError(
            "NCBI COMPLETE dataset manifest field 'consolidated' must be an object."
        )


def _validate_manifest_integer(
    value: object,
    field_path: str,
    *,
    minimum: int,
    maximum: int | None = None,
) -> int:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < minimum
        or (maximum is not None and value > maximum)
    ):
        if maximum is None:
            bounds = f"at least {minimum}"
        else:
            bounds = f"between {minimum} and {maximum}"
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' must be an integer {bounds}."
        )
    return value


def _validate_complete_source(
    value: object,
) -> tuple[Literal["uid", "accession"], tuple[str, ...]]:
    if not isinstance(value, dict):
        raise ValueError("NCBI dataset manifest field 'source' must be an object.")
    mode = value.get("mode")
    if mode == "accessions":
        source = _require_exact_fields(value, _ACCESSION_SOURCE_FIELDS, "source")
        if source["database"] != "nuccore":
            raise ValueError(
                "NCBI dataset manifest field 'source.database' must be 'nuccore'."
            )
        requested = _validate_string_list(
            source["requested_accessions"],
            "source.requested_accessions",
            allow_empty=False,
            require_unique=True,
        )
        return "accession", requested
    if mode == "query":
        source = _require_exact_fields(value, _QUERY_SOURCE_FIELDS, "source")
        if source["database"] != "nuccore":
            raise ValueError(
                "NCBI dataset manifest field 'source.database' must be 'nuccore'."
            )
        if not isinstance(source["query"], str) or not source["query"].strip():
            raise ValueError(
                "NCBI dataset manifest field 'source.query' must be a non-blank string."
            )
        uids = _validate_string_list(
            source["resolved_uids"],
            "source.resolved_uids",
            allow_empty=True,
            require_unique=True,
        )
        reported_count = _validate_manifest_integer(
            source["reported_count"], "source.reported_count", minimum=0
        )
        selected_count = _validate_manifest_integer(
            source["selected_count"], "source.selected_count", minimum=0
        )
        maximum = source["max_records"]
        if maximum is not None:
            maximum = _validate_manifest_integer(
                maximum, "source.max_records", minimum=1
            )
        if not isinstance(source["query_translation"], str):
            raise ValueError(
                "NCBI dataset manifest field 'source.query_translation' must be a string."
            )
        expected_selected = min(
            reported_count, maximum if maximum is not None else reported_count
        )
        if selected_count != len(uids) or selected_count != expected_selected:
            raise ValueError(
                "NCBI dataset manifest query source selected_count, reported_count, "
                "max_records, and resolved_uids composition are inconsistent."
            )
        return "uid", uids
    if "mode" not in value:
        raise ValueError(
            "NCBI dataset manifest field 'source.mode' is required."
        )
    raise ValueError(
        "NCBI dataset manifest field 'source.mode' must be 'accessions' or 'query'."
    )


def _validate_string_list(
    value: object,
    field_path: str,
    *,
    allow_empty: bool,
    require_unique: bool,
) -> tuple[str, ...]:
    if not isinstance(value, list) or (
        not allow_empty and not value
    ) or any(not isinstance(item, str) or not item.strip() for item in value):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' must be "
            f"a {'possibly empty' if allow_empty else 'non-empty'} list of non-blank strings."
        )
    if require_unique and len(set(value)) != len(value):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' must contain unique strings."
        )
    return tuple(value)


def _validate_complete_resolved_entries(
    value: object,
    identifiers: tuple[str, ...],
    identifier_kind: Literal["uid", "accession"],
    field_path: str,
) -> tuple[list[dict[str, object]], tuple[str, ...]]:
    if not isinstance(value, list):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' must be a list."
        )
    if len(value) != len(identifiers):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' does not match the source composition."
        )

    entries: list[dict[str, object]] = []
    association_ids: list[str] = []
    accession_versions: list[str] = []
    for index, raw_entry in enumerate(value):
        entry_path = f"{field_path}[{index}]"
        entry = _require_exact_fields(
            raw_entry, _RESOLVED_ENTRY_FIELDS, entry_path
        )
        requested_accession = entry["requested_accession"]
        uid = entry["uid"]
        if identifier_kind == "accession":
            if not isinstance(requested_accession, str) or not requested_accession.strip():
                raise ValueError(
                    f"NCBI dataset manifest field '{entry_path}.requested_accession' "
                    "must be a non-blank string for an accession source."
                )
            if uid is not None:
                raise ValueError(
                    f"NCBI dataset manifest field '{entry_path}.uid' must be null "
                    "for an accession source."
                )
            association_ids.append(requested_accession)
        else:
            if requested_accession is not None:
                raise ValueError(
                    f"NCBI dataset manifest field '{entry_path}.requested_accession' "
                    "must be null for a query source."
                )
            if not isinstance(uid, str) or not uid.strip():
                raise ValueError(
                    f"NCBI dataset manifest field '{entry_path}.uid' must be a "
                    "non-blank string for a query source."
                )
            association_ids.append(uid)

        accession = entry["accession"]
        accession_version = entry["accession_version"]
        if not isinstance(accession, str) or not accession.strip():
            raise ValueError(
                f"NCBI dataset manifest field '{entry_path}.accession' must be a "
                "non-blank string."
            )
        if not isinstance(accession_version, str) or not accession_version.strip():
            raise ValueError(
                f"NCBI dataset manifest field '{entry_path}.accession_version' must be "
                "a non-blank string."
            )
        if accession != re.sub(r"\.\d+$", "", accession_version):
            raise ValueError(
                f"NCBI dataset manifest field '{entry_path}.accession' does not match "
                f"'{entry_path}.accession_version'."
            )
        entries.append(entry)
        accession_versions.append(accession_version)

    if tuple(association_ids) != identifiers:
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' associations and "
            "accession_version order do not match the source composition."
        )
    if len(set(accession_versions)) != len(accession_versions):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}[*].accession_version' "
            "must be unique."
        )
    return entries, tuple(accession_versions)


def _validate_complete_batches(
    directory: Path,
    value: object,
    identifiers: tuple[str, ...],
    resolved_entries: list[dict[str, object]],
    accession_versions: tuple[str, ...],
    identifier_kind: Literal["uid", "accession"],
    batch_size: object,
) -> None:
    if not isinstance(value, list):
        raise ValueError(
            "NCBI dataset manifest field 'completed_batches' must be a list."
        )
    if isinstance(batch_size, bool) or not isinstance(batch_size, int):
        raise ValueError(
            "NCBI dataset manifest field 'batch_size' must be an integer."
        )
    planned = _planned_batches(identifiers, batch_size)
    if len(value) != len(planned):
        raise ValueError(
            "NCBI dataset manifest field 'completed_batches' does not cover the "
            "complete source composition."
        )

    entry_offset = 0
    for batch_index, expected_identifiers in enumerate(planned):
        batch_path = f"completed_batches[{batch_index}]"
        metadata = _require_exact_fields(
            value[batch_index], _COMPLETED_BATCH_FIELDS, batch_path
        )
        expected_filename = _batch_filename(batch_index)
        if metadata["filename"] != expected_filename:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.filename' must be "
                f"'{expected_filename}'."
            )
        requested = _validate_string_list(
            metadata["requested_identifiers"],
            f"{batch_path}.requested_identifiers",
            allow_empty=False,
            require_unique=True,
        )
        if requested != expected_identifiers:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.requested_identifiers' "
                "does not match the source composition."
            )
        record_count = _validate_manifest_integer(
            metadata["record_count"], f"{batch_path}.record_count", minimum=0
        )
        if record_count != len(expected_identifiers):
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.record_count' does not "
                "match its requested identifiers."
            )
        _validate_manifest_integer(
            metadata["byte_size"], f"{batch_path}.byte_size", minimum=0
        )
        _validate_checksum(metadata["sha256"], f"{batch_path}.sha256")
        batch_entries, batch_versions = _validate_complete_resolved_entries(
            metadata["resolved_entries"],
            expected_identifiers,
            identifier_kind,
            f"{batch_path}.resolved_entries",
        )
        expected_entries = resolved_entries[
            entry_offset : entry_offset + len(expected_identifiers)
        ]
        expected_versions = accession_versions[
            entry_offset : entry_offset + len(expected_identifiers)
        ]
        if batch_entries != expected_entries:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.resolved_entries' does not "
                "match the top-level resolved_entries composition."
            )
        record_ids = _validate_string_list(
            metadata["record_ids"],
            f"{batch_path}.record_ids",
            allow_empty=False,
            require_unique=True,
        )
        if record_ids != expected_versions or batch_versions != expected_versions:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.record_ids' does not match "
                "its resolved_entries."
            )

        artifact_path = directory / _BATCH_DIRECTORY_NAME / expected_filename
        if not artifact_path.is_file():
            raise ValueError(
                f"NCBI dataset batch file '{expected_filename}' is missing."
            )
        artifact_bytes = artifact_path.read_bytes()
        _validate_artifact_bytes(
            metadata, artifact_bytes, batch_path, expected_filename
        )
        parsed_records = _parse_genbank_records(artifact_path)
        parsed_ids = tuple(record.id for record in parsed_records)
        if len(parsed_records) != record_count:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.record_count' does not "
                f"match '{expected_filename}'."
            )
        if parsed_ids != record_ids:
            raise ValueError(
                f"NCBI dataset manifest field '{batch_path}.record_ids' does not "
                f"match '{expected_filename}'."
            )
        entry_offset += len(expected_identifiers)


def _validate_complete_consolidated(
    value: object, expected_count: int
) -> dict[str, object]:
    consolidated = _require_exact_fields(
        value, _CONSOLIDATED_FIELDS, "consolidated"
    )
    if consolidated["filename"] != RECORDS_NAME:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.filename' must be "
            f"'{RECORDS_NAME}'."
        )
    record_count = _validate_manifest_integer(
        consolidated["record_count"], "consolidated.record_count", minimum=0
    )
    if record_count != expected_count:
        raise ValueError(
            "NCBI dataset manifest field 'consolidated.record_count' does not match "
            "the source composition."
        )
    _validate_manifest_integer(
        consolidated["byte_size"], "consolidated.byte_size", minimum=0
    )
    _validate_checksum(consolidated["sha256"], "consolidated.sha256")
    return consolidated


def _validate_checksum(value: object, field_path: str) -> str:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}' must be a lowercase "
            "SHA-256 digest."
        )
    return value


def _validate_artifact_bytes(
    metadata: dict[str, object], contents: bytes, field_path: str, filename: str
) -> None:
    if metadata["byte_size"] != len(contents):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}.byte_size' does not match "
            f"'{filename}'."
        )
    if metadata["sha256"] != _sha256_bytes(contents):
        raise ValueError(
            f"NCBI dataset manifest field '{field_path}.sha256' does not match "
            f"'{filename}'."
        )
