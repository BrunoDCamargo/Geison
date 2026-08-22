import hashlib
import io
import json
import os
import tempfile
import threading
import unittest
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError, URLError

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import qpcr_pipeline.ncbi as ncbi
from qpcr_pipeline.config import NcbiInputConfig
from qpcr_pipeline.ncbi import (
    NcbiFetchedRecord,
    NcbiTransientError,
    ResolvedNcbiQuery,
    acquire_ncbi_dataset,
    validate_frozen_dataset,
)


class FakeNcbiClient:
    def __init__(self, records_by_request, failures=(), resolutions=()):
        self.records_by_request = records_by_request
        self.failures = list(failures)
        self.resolutions = list(resolutions)
        self.fetch_calls = []
        self.resolve_calls = []

    def resolve_query(self, query, max_records):
        self.resolve_calls.append((query, max_records))
        if self.resolutions:
            resolution = self.resolutions.pop(0)
            if isinstance(resolution, BaseException):
                raise resolution
            return resolution
        raise AssertionError("query resolution was not expected")

    def fetch_records(self, identifiers, *, identifier_kind):
        requested = tuple(identifiers)
        self.fetch_calls.append((requested, identifier_kind))
        if self.failures:
            failure = self.failures.pop(0)
            if failure is not None:
                raise failure
        return tuple(
            NcbiFetchedRecord(request_id=identifier, record=self.records_by_request[identifier])
            for identifier in reversed(requested)
        )


class AccessionAcquisitionTests(unittest.TestCase):
    def _record(self, accession_version: str) -> SeqRecord:
        record = SeqRecord(
            Seq("ATGCATGC"), id=accession_version, description=f"{accession_version} record"
        )
        record.annotations["molecule_type"] = "DNA"
        return record

    def _config(self, accessions: tuple[str, ...], *, batch_size: int = 2) -> NcbiInputConfig:
        return NcbiInputConfig(accessions=accessions, batch_size=batch_size, retries=3)

    def _read_manifest(self, directory: Path) -> dict[str, object]:
        return json.loads((directory / "dataset_manifest.json").read_text(encoding="utf-8"))

    def test_fetches_accessions_in_stable_batches_and_consolidates_requested_order(self):
        accessions = ("NC_000001.11", "AB123456.2", "XY_3.4", "MN_5.6", "PQ_7.8")
        client = FakeNcbiClient({accession: self._record(accession) for accession in accessions})

        with tempfile.TemporaryDirectory() as tmpdir:
            dataset = acquire_ncbi_dataset(
                self._config(accessions), Path(tmpdir), client=client, clock=lambda: "now"
            )
            manifest = self._read_manifest(Path(tmpdir))

            self.assertEqual(
                client.fetch_calls,
                [
                    (("NC_000001.11", "AB123456.2"), "accession"),
                    (("XY_3.4", "MN_5.6"), "accession"),
                    (("PQ_7.8",), "accession"),
                ],
            )
            self.assertEqual(client.resolve_calls, [])
            self.assertEqual(
                [entry["requested_accession"] for entry in manifest["resolved_entries"]],
                list(accessions),
            )
            self.assertEqual(
                [entry["accession_version"] for entry in manifest["resolved_entries"]],
                list(accessions),
            )
            with dataset.records_path.open(encoding="utf-8") as handle:
                self.assertEqual(
                    tuple(record.id for record in SeqIO.parse(handle, "genbank")), accessions
                )

    def test_records_secret_free_manifest_and_final_batch_and_consolidated_checksums(self):
        accessions = ("NC_000001.11", "AB123456.2", "XY_3.4")
        client = FakeNcbiClient({accession: self._record(accession) for accession in accessions})

        with tempfile.TemporaryDirectory() as tmpdir, patch.dict(
            os.environ,
            {"NCBI_EMAIL": "private@example.test", "NCBI_API_KEY": "private-key"},
        ):
            directory = Path(tmpdir)
            dataset = acquire_ncbi_dataset(
                self._config(accessions), directory, client=client, clock=lambda: "now"
            )
            manifest = self._read_manifest(directory)
            serialized = json.dumps(manifest)

            self.assertEqual(manifest["status"], "COMPLETE")
            self.assertEqual(
                manifest["source"],
                {
                    "mode": "accessions",
                    "database": "nuccore",
                    "requested_accessions": list(accessions),
                },
            )
            self.assertNotIn("NCBI_EMAIL", serialized)
            self.assertNotIn("NCBI_API_KEY", serialized)
            self.assertNotIn("private@example.test", serialized)
            self.assertNotIn("private-key", serialized)
            self.assertEqual(manifest["created_at"], "now")
            self.assertEqual(manifest["updated_at"], "now")
            self.assertEqual(len(manifest["completed_batches"]), 2)
            for batch in manifest["completed_batches"]:
                batch_path = directory / "batches" / batch["filename"]
                batch_bytes = batch_path.read_bytes()
                self.assertEqual(batch["byte_size"], len(batch_bytes))
                self.assertEqual(batch["sha256"], hashlib.sha256(batch_bytes).hexdigest())
            consolidated = manifest["consolidated"]
            record_bytes = dataset.records_path.read_bytes()
            self.assertEqual(consolidated["filename"], "records.gb")
            self.assertEqual(consolidated["record_count"], 3)
            self.assertEqual(consolidated["byte_size"], len(record_bytes))
            self.assertEqual(consolidated["sha256"], hashlib.sha256(record_bytes).hexdigest())

    def test_retries_only_transient_failures_with_bounded_exponential_delays(self):
        accessions = ("NC_000001.11",)
        client = FakeNcbiClient(
            {accessions[0]: self._record(accessions[0])},
            failures=(NcbiTransientError("temporary"), NcbiTransientError("temporary")),
        )
        sleeps: list[int] = []

        with tempfile.TemporaryDirectory() as tmpdir:
            acquire_ncbi_dataset(
                self._config(accessions),
                Path(tmpdir),
                client=client,
                sleep=sleeps.append,
                clock=lambda: "now",
            )

        self.assertEqual(sleeps, [1, 2])
        self.assertEqual(client.fetch_calls, [((accessions[0],), "accession")] * 3)

    def test_reports_sanitized_fetch_operation_and_attempts_after_retry_exhaustion(self):
        accessions = ("NC_000001.11",)
        secret = "private-api-key"
        client = FakeNcbiClient(
            {accessions[0]: self._record(accessions[0])},
            failures=(NcbiTransientError(secret),) * 4,
        )

        with tempfile.TemporaryDirectory() as tmpdir, self.assertRaisesRegex(
            NcbiTransientError, "NCBI record fetch failed after 4 attempts"
        ) as raised:
            acquire_ncbi_dataset(
                self._config(accessions),
                Path(tmpdir),
                client=client,
                sleep=lambda _: None,
                clock=lambda: "now",
            )

        self.assertEqual(client.fetch_calls, [((accessions[0],), "accession")] * 4)
        self.assertNotIn(secret, str(raised.exception))
        self.assertIsNone(raised.exception.__cause__)
        self.assertTrue(raised.exception.__suppress_context__)

    def test_propagates_non_transient_failures_without_retrying(self):
        accessions = ("NC_000001.11",)
        client = FakeNcbiClient(
            {accessions[0]: self._record(accessions[0])}, failures=(ValueError("invalid"),)
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "invalid"):
                acquire_ncbi_dataset(
                    self._config(accessions), Path(tmpdir), client=client, clock=lambda: "now"
                )
            manifest = self._read_manifest(Path(tmpdir))

        self.assertEqual(client.fetch_calls, [((accessions[0],), "accession")])
        self.assertEqual(manifest["status"], "PARTIAL")
        self.assertEqual(manifest["completed_batches"], [])

    def test_initial_partial_manifest_has_independent_creation_and_update_timestamps(self):
        accessions = ("NC_000001.11",)
        client = FakeNcbiClient(
            {accessions[0]: self._record(accessions[0])}, failures=(ValueError("stopped"),)
        )
        timestamps = iter(("created", "initialized"))

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "stopped"):
                acquire_ncbi_dataset(
                    self._config(accessions),
                    Path(tmpdir),
                    client=client,
                    clock=lambda: next(timestamps),
                )
            manifest = self._read_manifest(Path(tmpdir))

        self.assertEqual(manifest["status"], "PARTIAL")
        self.assertEqual(manifest["created_at"], "created")
        self.assertEqual(manifest["updated_at"], "initialized")

    def test_resumes_after_interruption_with_the_first_completed_batch_reusable(self):
        accessions = ("NC_000001.11", "AB123456.2", "XY_3.4")
        records = {accession: self._record(accession) for accession in accessions}
        interrupted_client = FakeNcbiClient(records, failures=(None, ValueError("stopped")))

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            with self.assertRaisesRegex(ValueError, "stopped"):
                acquire_ncbi_dataset(
                    self._config(accessions),
                    directory,
                    client=interrupted_client,
                    clock=lambda: "2026-08-21T00:00:00+00:00",
                )
            partial = self._read_manifest(directory)
            resumed_client = FakeNcbiClient(records)
            acquire_ncbi_dataset(
                self._config(accessions),
                directory,
                client=resumed_client,
                clock=lambda: "2026-08-21T01:00:00+00:00",
            )
            complete = self._read_manifest(directory)

            self.assertEqual(partial["status"], "PARTIAL")
            self.assertEqual(partial["consolidated"], None)
            self.assertEqual(partial["source"]["requested_accessions"], list(accessions))
            self.assertEqual(len(partial["completed_batches"]), 1)
            self.assertEqual(
                resumed_client.fetch_calls, [(("XY_3.4",), "accession")]
            )
            self.assertEqual(complete["status"], "COMPLETE")
            self.assertEqual(complete["created_at"], "2026-08-21T00:00:00+00:00")
            self.assertEqual(complete["updated_at"], "2026-08-21T01:00:00+00:00")

    def test_selectively_refetches_a_corrupted_completed_batch(self):
        accessions = ("NC_000001.11", "AB123456.2", "XY_3.4", "MN_5.6")
        records = {accession: self._record(accession) for accession in accessions}

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            acquire_ncbi_dataset(
                self._config(accessions),
                directory,
                client=FakeNcbiClient(records),
                clock=lambda: "2026-08-21T00:00:00+00:00",
            )
            (directory / "batches" / "batch-00001.gb").write_text("corrupted", encoding="utf-8")
            refetch_client = FakeNcbiClient(records)
            acquire_ncbi_dataset(
                self._config(accessions),
                directory,
                client=refetch_client,
                clock=lambda: "2026-08-21T01:00:00+00:00",
            )

            self.assertEqual(
                refetch_client.fetch_calls, [(("XY_3.4", "MN_5.6"), "accession")]
            )
            self.assertEqual(self._read_manifest(directory)["status"], "COMPLETE")

    def test_rejects_ambiguous_client_response_composition(self):
        accessions = ("NC_000001.11", "AB123456.2")
        records = {accession: self._record(accession) for accession in accessions}

        class DuplicateResponseClient(FakeNcbiClient):
            def fetch_records(self, identifiers, *, identifier_kind):
                requested = tuple(identifiers)
                self.fetch_calls.append((requested, identifier_kind))
                return (
                    NcbiFetchedRecord(request_id=requested[0], record=records[requested[0]]),
                    NcbiFetchedRecord(request_id=requested[0], record=records[requested[1]]),
                )

        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "exactly one.*requested identifier"):
                acquire_ncbi_dataset(
                    self._config(accessions),
                    Path(tmpdir),
                    client=DuplicateResponseClient(records),
                    clock=lambda: "now",
                )

    def test_default_client_requires_environment_email_before_creating_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset_dir = Path(tmpdir) / "dataset"
            with patch.dict(os.environ, {}, clear=True), self.assertRaisesRegex(
                ValueError, "NCBI_EMAIL"
            ):
                acquire_ncbi_dataset(
                    self._config(("NC_000001.11",)), dataset_dir, clock=lambda: "now"
                )
            self.assertFalse(dataset_dir.exists())

    def test_rejects_noncanonical_resume_manifests_before_reuse(self):
        accessions = ("NC_000001.11",)
        records = {accessions[0]: self._record(accessions[0])}
        cases = (
            ("credential field", "NCBI_API_KEY", "secret", "unrecognized"),
            ("floating schema version", "schema_version", 1.0, "schema"),
            ("floating batch size", "batch_size", 2.0, "batch_size"),
            ("floating retries", "retries", 3.0, "retries"),
            ("wrong tool", "tool", "other-tool", "tool"),
            ("malformed update timestamp", "updated_at", "not-a-timestamp", "updated_at"),
            ("partial with consolidated records", "status", "PARTIAL", "consolidated"),
            ("complete without consolidated records", "consolidated", None, "consolidated"),
        )

        for name, field, value, error in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                directory = Path(tmpdir)
                acquire_ncbi_dataset(
                    self._config(accessions),
                    directory,
                    client=FakeNcbiClient(records),
                    clock=lambda: "2026-08-21T00:00:00+00:00",
                )
                manifest = self._read_manifest(directory)
                manifest[field] = value
                (directory / "dataset_manifest.json").write_text(
                    json.dumps(manifest), encoding="utf-8"
                )
                resumed_client = FakeNcbiClient(records)

                with self.assertRaisesRegex(ValueError, error):
                    acquire_ncbi_dataset(
                        self._config(accessions),
                        directory,
                        client=resumed_client,
                        clock=lambda: "2026-08-21T01:00:00+00:00",
                    )

                self.assertEqual(resumed_client.fetch_calls, [])
                self.assertEqual(self._read_manifest(directory), manifest)

    def test_rejects_invalid_direct_accession_configs_before_creating_output(self):
        records = {
            "A": self._record("A.1"),
            "B": self._record("B.1"),
            " ": self._record("BLANK.1"),
        }
        cases = (
            ("query", NcbiInputConfig(accessions=("A",), query="virus")),
            ("frozen dataset", NcbiInputConfig(accessions=("A",), frozen_dataset=Path("frozen"))),
            ("maximum records", NcbiInputConfig(accessions=("A",), max_records=1)),
            ("duplicate accessions", NcbiInputConfig(accessions=("A", "A"))),
            ("blank accession", NcbiInputConfig(accessions=(" ",))),
            ("empty accessions", NcbiInputConfig(accessions=())),
            ("small batch", NcbiInputConfig(accessions=("A",), batch_size=0)),
            ("large batch", NcbiInputConfig(accessions=("A",), batch_size=501)),
            ("negative retries", NcbiInputConfig(accessions=("A",), retries=-1)),
            ("large retries", NcbiInputConfig(accessions=("A",), retries=11)),
        )

        for name, config in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                dataset_dir = Path(tmpdir) / "dataset"
                client = FakeNcbiClient(records)

                with self.assertRaisesRegex(ValueError, "accession|batch_size|retries"):
                    acquire_ncbi_dataset(
                        config,
                        dataset_dir,
                        client=client,
                        clock=lambda: "2026-08-21T00:00:00+00:00",
                    )

                self.assertFalse(dataset_dir.exists())
                self.assertEqual(client.fetch_calls, [])


class QueryAcquisitionTests(unittest.TestCase):
    def _record(self, accession_version: str) -> SeqRecord:
        record = SeqRecord(
            Seq("ATGCATGC"), id=accession_version, description=f"{accession_version} record"
        )
        record.annotations["molecule_type"] = "DNA"
        return record

    def _config(
        self,
        query: str = "example[Organism]",
        *,
        batch_size: int = 2,
        retries: int = 3,
        max_records: int | None = None,
    ) -> NcbiInputConfig:
        return NcbiInputConfig(
            query=query,
            batch_size=batch_size,
            retries=retries,
            max_records=max_records,
        )

    def _resolution(self, uids=("101", "102"), *, reported_count=2):
        return ResolvedNcbiQuery(
            uids=tuple(uids),
            reported_count=reported_count,
            query_translation="example organism[Organism]",
        )

    def _records(self):
        return {"101": self._record("NC_1.2"), "102": self._record("NC_2.3")}

    def _read_manifest(self, directory: Path) -> dict[str, object]:
        return json.loads((directory / "dataset_manifest.json").read_text(encoding="utf-8"))

    def test_persists_complete_query_composition_before_fetching(self):
        client = FakeNcbiClient(
            self._records(),
            failures=(ValueError("stopped"),),
            resolutions=(self._resolution(),),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            with self.assertRaisesRegex(ValueError, "stopped"):
                acquire_ncbi_dataset(
                    self._config(), directory, client=client, clock=lambda: "now"
                )
            manifest = self._read_manifest(directory)

        self.assertEqual(client.resolve_calls, [("example[Organism]", None)])
        self.assertEqual(manifest["status"], "PARTIAL")
        self.assertEqual(manifest["source"]["resolved_uids"], ["101", "102"])
        self.assertEqual(manifest["source"]["reported_count"], 2)
        self.assertEqual(manifest["source"]["selected_count"], 2)
        self.assertEqual(manifest["completed_batches"], [])

    def test_applies_maximum_prefix_and_preserves_uid_to_accession_mapping(self):
        resolution = self._resolution(("101", "102"), reported_count=9)
        client = FakeNcbiClient(self._records(), resolutions=(resolution,))

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            acquire_ncbi_dataset(
                self._config(max_records=2),
                directory,
                client=client,
                clock=lambda: "now",
            )
            manifest = self._read_manifest(directory)

        self.assertEqual(client.resolve_calls, [("example[Organism]", 2)])
        self.assertEqual(client.fetch_calls, [(("101", "102"), "uid")])
        self.assertEqual(
            manifest["source"],
            {
                "mode": "query",
                "database": "nuccore",
                "query": "example[Organism]",
                "resolved_uids": ["101", "102"],
                "reported_count": 9,
                "selected_count": 2,
                "max_records": 2,
                "query_translation": "example organism[Organism]",
            },
        )
        self.assertEqual(
            [(entry["uid"], entry["accession_version"]) for entry in manifest["resolved_entries"]],
            [("101", "NC_1.2"), ("102", "NC_2.3")],
        )
        self.assertEqual(
            [entry["requested_accession"] for entry in manifest["resolved_entries"]],
            [None, None],
        )

    def test_resumes_partial_query_without_resolving_again(self):
        records = self._records()
        interrupted = FakeNcbiClient(
            records,
            failures=(None, ValueError("stopped")),
            resolutions=(self._resolution(),),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            with self.assertRaisesRegex(ValueError, "stopped"):
                acquire_ncbi_dataset(
                    self._config(batch_size=1),
                    directory,
                    client=interrupted,
                    clock=lambda: "2026-08-21T00:00:00+00:00",
                )
            resumed = FakeNcbiClient(records)
            acquire_ncbi_dataset(
                self._config(batch_size=1),
                directory,
                client=resumed,
                clock=lambda: "2026-08-21T01:00:00+00:00",
            )

        self.assertEqual(resumed.resolve_calls, [])
        self.assertEqual(resumed.fetch_calls, [(("102",), "uid")])

    def test_rejects_incompatible_partial_query_before_network_or_output_mutation(self):
        cases = (
            ("query", self._config(query="other[Organism]"), "source"),
            ("maximum", self._config(max_records=1), "source"),
            ("batch size", self._config(batch_size=1), "batch_size"),
            ("retry limit", self._config(retries=2), "retries"),
        )

        for name, incompatible, error in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                directory = Path(tmpdir)
                creator = FakeNcbiClient(
                    self._records(),
                    failures=(ValueError("stopped"),),
                    resolutions=(self._resolution(),),
                )
                with self.assertRaisesRegex(ValueError, "stopped"):
                    acquire_ncbi_dataset(
                        self._config(), directory, client=creator, clock=lambda: "now"
                    )
                before = (directory / "dataset_manifest.json").read_bytes()
                resumed = FakeNcbiClient(
                    self._records(), resolutions=(self._resolution(),)
                )

                with self.assertRaisesRegex(ValueError, error):
                    acquire_ncbi_dataset(
                        incompatible, directory, client=resumed, clock=lambda: "later"
                    )

                self.assertEqual(resumed.resolve_calls, [])
                self.assertEqual(resumed.fetch_calls, [])
                self.assertEqual((directory / "dataset_manifest.json").read_bytes(), before)

    def test_rejects_invalid_direct_query_configs_before_creating_output(self):
        cases = (
            ("missing query", NcbiInputConfig(), "query"),
            ("blank query", NcbiInputConfig(query=" "), "query"),
            (
                "accessions",
                NcbiInputConfig(query="virus", accessions=("NC_1",)),
                "query",
            ),
            (
                "noncanonical accession list",
                NcbiInputConfig(query="virus", accessions=["NC_1"]),
                "accession",
            ),
            (
                "frozen dataset",
                NcbiInputConfig(query="virus", frozen_dataset=Path("frozen")),
                "query",
            ),
            ("zero maximum", NcbiInputConfig(query="virus", max_records=0), "max_records"),
            ("boolean maximum", NcbiInputConfig(query="virus", max_records=True), "max_records"),
            ("small batch", NcbiInputConfig(query="virus", batch_size=0), "batch_size"),
            ("large batch", NcbiInputConfig(query="virus", batch_size=501), "batch_size"),
            ("negative retries", NcbiInputConfig(query="virus", retries=-1), "retries"),
            ("large retries", NcbiInputConfig(query="virus", retries=11), "retries"),
        )

        for name, config, error in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                dataset_dir = Path(tmpdir) / "dataset"
                client = FakeNcbiClient({}, resolutions=(self._resolution(),))

                with self.assertRaisesRegex(ValueError, error):
                    acquire_ncbi_dataset(config, dataset_dir, client=client, clock=lambda: "now")

                self.assertFalse(dataset_dir.exists())
                self.assertEqual(client.resolve_calls, [])
                self.assertEqual(client.fetch_calls, [])

    def test_rejects_noncanonical_query_source_fields_before_reuse(self):
        cases = (
            ("floating reported count", "reported_count", 2.0, "composition"),
            ("floating selected count", "selected_count", 2.0, "selected_count"),
            ("floating maximum", "max_records", 2.0, "source"),
            ("credential field", "NCBI_API_KEY", "secret", "canonical"),
        )

        for name, field, value, error in cases:
            with self.subTest(name=name), tempfile.TemporaryDirectory() as tmpdir:
                directory = Path(tmpdir)
                creator = FakeNcbiClient(
                    self._records(),
                    failures=(ValueError("stopped"),),
                    resolutions=(self._resolution(),),
                )
                config = self._config(max_records=2)
                with self.assertRaisesRegex(ValueError, "stopped"):
                    acquire_ncbi_dataset(
                        config,
                        directory,
                        client=creator,
                        clock=lambda: "2026-08-21T00:00:00+00:00",
                    )
                manifest = self._read_manifest(directory)
                manifest["source"][field] = value
                (directory / "dataset_manifest.json").write_text(
                    json.dumps(manifest), encoding="utf-8"
                )
                resumed = FakeNcbiClient(self._records())

                with self.assertRaisesRegex(ValueError, error):
                    acquire_ncbi_dataset(
                        config, directory, client=resumed, clock=lambda: "later"
                    )

                self.assertEqual(resumed.resolve_calls, [])
                self.assertEqual(resumed.fetch_calls, [])

    def test_rejects_nonlist_resolved_uids_even_for_a_zero_count_query(self):
        creator = FakeNcbiClient(
            self._records(),
            failures=(ValueError("stopped"),),
            resolutions=(self._resolution(),),
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            with self.assertRaisesRegex(ValueError, "stopped"):
                acquire_ncbi_dataset(
                    self._config(),
                    directory,
                    client=creator,
                    clock=lambda: "2026-08-21T00:00:00+00:00",
                )
            manifest = self._read_manifest(directory)
            manifest["source"]["resolved_uids"] = None
            manifest["source"]["reported_count"] = 0
            manifest["source"]["selected_count"] = 0
            (directory / "dataset_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )
            resumed = FakeNcbiClient(self._records())

            with self.assertRaisesRegex(ValueError, "resolved_uids|composition"):
                acquire_ncbi_dataset(
                    self._config(), directory, client=resumed, clock=lambda: "later"
                )

            self.assertEqual(resumed.resolve_calls, [])
            self.assertEqual(resumed.fetch_calls, [])

    def test_reports_sanitized_query_operation_and_attempts_after_retry_exhaustion(self):
        secret = "private-api-key"
        client = FakeNcbiClient(
            self._records(),
            resolutions=(NcbiTransientError(secret),) * 3,
        )

        with tempfile.TemporaryDirectory() as tmpdir, self.assertRaisesRegex(
            NcbiTransientError, "NCBI query resolution failed after 3 attempts"
        ) as raised:
            acquire_ncbi_dataset(
                self._config(retries=2),
                Path(tmpdir),
                client=client,
                sleep=lambda _: None,
                clock=lambda: "now",
            )

        self.assertEqual(client.resolve_calls, [("example[Organism]", None)] * 3)
        self.assertNotIn(secret, str(raised.exception))
        self.assertIsNone(raised.exception.__cause__)
        self.assertTrue(raised.exception.__suppress_context__)


class _SearchHandle:
    def __init__(self, payload):
        self.payload = payload
        self.closed = False

    def close(self):
        self.closed = True


class _FakeEntrezModule:
    def __init__(self, uids=("101", "102"), *, reported_count=None):
        self.email = None
        self.tool = None
        self.api_key = None
        self.uids = tuple(uids)
        self.reported_count = len(self.uids) if reported_count is None else reported_count
        self.search_handles = []
        self.fetch_handles = []
        self.search_calls = []
        self.fetch_calls = []
        self.search_error = None
        self.read_error = None
        self.fetch_error = None
        self.fetch_text = ""

    def esearch(self, **kwargs):
        self.search_calls.append((kwargs, self.email, self.tool, self.api_key))
        if self.search_error is not None:
            raise self.search_error
        start = kwargs["retstart"]
        stop = start + kwargs["retmax"]
        handle = _SearchHandle(
            {
                "Count": str(self.reported_count),
                "IdList": list(self.uids[start:stop]),
                "QueryTranslation": "translated query",
            }
        )
        self.search_handles.append(handle)
        return handle

    def read(self, handle):
        if self.read_error is not None:
            raise self.read_error
        return handle.payload

    def efetch(self, **kwargs):
        self.fetch_calls.append((kwargs, self.email, self.tool, self.api_key))
        if self.fetch_error is not None:
            raise self.fetch_error
        handle = io.StringIO(self.fetch_text)
        self.fetch_handles.append(handle)
        return handle


class _CoordinatedEntrezModule(_FakeEntrezModule):
    def __init__(self):
        super().__init__()
        self.first_request_entered = threading.Event()
        self.release_first_request = threading.Event()
        self.second_email_configured = threading.Event()
        self._configuration_count = 0
        self._request_count = 0
        self._coordination_lock = threading.Lock()
        self._track_configuration = True

    def __setattr__(self, name, value):
        object.__setattr__(self, name, value)
        if name == "email" and getattr(self, "_track_configuration", False):
            with self._coordination_lock:
                self._configuration_count += 1
                if self._configuration_count == 2:
                    self.second_email_configured.set()

    def esearch(self, **kwargs):
        with self._coordination_lock:
            self._request_count += 1
            first_request = self._request_count == 1
        if first_request:
            self.first_request_entered.set()
            if not self.release_first_request.wait(timeout=5):
                raise AssertionError("test did not release the first request")
        return super().esearch(**kwargs)


class BioEntrezClientTests(unittest.TestCase):
    def _genbank(self, *accession_versions: str) -> str:
        records = []
        for accession_version in accession_versions:
            record = SeqRecord(Seq("ATGCATGC"), id=accession_version, description="record")
            record.annotations["molecule_type"] = "DNA"
            records.append(record)
        handle = io.StringIO()
        SeqIO.write(records, handle, "genbank")
        return handle.getvalue()

    def _client(self, module, *, api_key="private-api-key"):
        return ncbi.BioEntrezClient(
            email="private@example.test",
            api_key=api_key,
            entrez_module=module,
        )

    def test_resolves_ordered_uids_in_stable_pages_of_at_most_ten_thousand(self):
        uids = tuple(str(index) for index in range(1, 10002))
        module = _FakeEntrezModule(uids)

        result = self._client(module).resolve_query("example[Organism]", None)

        self.assertEqual(result.uids, uids)
        self.assertEqual(result.reported_count, 10001)
        self.assertEqual(result.query_translation, "translated query")
        self.assertEqual(
            [call[0] for call in module.search_calls],
            [
                {
                    "db": "nuccore",
                    "term": "example[Organism]",
                    "retstart": 0,
                    "retmax": 10000,
                },
                {
                    "db": "nuccore",
                    "term": "example[Organism]",
                    "retstart": 10000,
                    "retmax": 1,
                },
            ],
        )
        self.assertTrue(all(handle.closed for handle in module.search_handles))
        self.assertNotIn("private-api-key", repr(result))

    def test_fetches_nuccore_records_and_maps_accessions_by_version_or_base(self):
        module = _FakeEntrezModule()
        module.fetch_text = self._genbank("AB2.3", "NC_1.2")

        result = self._client(module).fetch_records(
            ("NC_1", "AB2.3"), identifier_kind="accession"
        )

        self.assertEqual(
            [(item.request_id, item.record.id) for item in result],
            [("NC_1", "NC_1.2"), ("AB2.3", "AB2.3")],
        )
        self.assertEqual(
            module.fetch_calls[0][0],
            {
                "db": "nuccore",
                "id": "NC_1,AB2.3",
                "rettype": "gb",
                "retmode": "text",
            },
        )
        self.assertTrue(module.fetch_handles[0].closed)

    def test_maps_uid_requests_by_returned_record_order(self):
        module = _FakeEntrezModule()
        module.fetch_text = self._genbank("NC_1.2", "NC_2.3")

        result = self._client(module).fetch_records(("101", "102"), identifier_kind="uid")

        self.assertEqual(
            [(item.request_id, item.record.id) for item in result],
            [("101", "NC_1.2"), ("102", "NC_2.3")],
        )

    def test_configures_credentials_immediately_before_every_request(self):
        module = _FakeEntrezModule()
        module.fetch_text = self._genbank("NC_1.2")
        client = self._client(module)

        client.resolve_query("example", 1)
        module.email = "changed"
        module.tool = "changed"
        module.api_key = "changed"
        client.fetch_records(("101",), identifier_kind="uid")

        for _, email, tool, api_key in module.search_calls + module.fetch_calls:
            self.assertEqual(email, "private@example.test")
            self.assertEqual(tool, "geison-qpcr")
            self.assertEqual(api_key, "private-api-key")

    def test_configuration_and_request_creation_are_atomic_across_clients(self):
        module = _CoordinatedEntrezModule()
        first = ncbi.BioEntrezClient(
            email="first@example.test",
            api_key="first-key",
            entrez_module=module,
        )
        second = ncbi.BioEntrezClient(
            email="second@example.test",
            api_key="second-key",
            entrez_module=module,
        )
        errors = []

        def resolve(client):
            try:
                client.resolve_query("example", 1)
            except BaseException as error:
                errors.append(error)

        first_thread = threading.Thread(target=resolve, args=(first,))
        second_thread = threading.Thread(target=resolve, args=(second,))
        first_thread.start()
        self.assertTrue(module.first_request_entered.wait(timeout=2))
        second_thread.start()
        module.second_email_configured.wait(timeout=1)
        module.release_first_request.set()
        first_thread.join(timeout=5)
        second_thread.join(timeout=5)

        self.assertFalse(first_thread.is_alive(), "first request deadlocked")
        self.assertFalse(second_thread.is_alive(), "second request deadlocked")
        self.assertEqual(errors, [])
        self.assertEqual(
            [call[1:] for call in module.search_calls],
            [
                ("first@example.test", "geison-qpcr", "first-key"),
                ("second@example.test", "geison-qpcr", "second-key"),
            ],
        )

    def test_environment_factory_requires_nonblank_email_before_any_request(self):
        module = _FakeEntrezModule()
        for environ in ({}, {"NCBI_EMAIL": "   "}):
            with self.subTest(environ=environ):
                with self.assertRaisesRegex(ValueError, "NCBI_EMAIL") as raised:
                    ncbi.BioEntrezClient.from_environment(
                        environ=environ, entrez_module=module
                    )
                self.assertNotIn("private-api-key", str(raised.exception))
        self.assertEqual(module.search_calls, [])
        self.assertEqual(module.fetch_calls, [])

    def test_environment_factory_treats_blank_api_key_as_absent(self):
        for api_key in ("", "   "):
            with self.subTest(api_key=api_key):
                module = _FakeEntrezModule()
                client = ncbi.BioEntrezClient.from_environment(
                    environ={
                        "NCBI_EMAIL": " person@example.test ",
                        "NCBI_API_KEY": api_key,
                    },
                    entrez_module=module,
                )

                result = client.resolve_query("example", 1)

                self.assertEqual(result.uids, ("101",))
                self.assertEqual(
                    module.search_calls[0][1:],
                    ("person@example.test", "geison-qpcr", None),
                )

    def test_client_representation_does_not_expose_credentials(self):
        representation = repr(self._client(_FakeEntrezModule()))

        self.assertNotIn("private@example.test", representation)
        self.assertNotIn("private-api-key", representation)

    def test_rejects_incomplete_search_composition(self):
        module = _FakeEntrezModule(("101",), reported_count=2)

        with self.assertRaisesRegex(ValueError, "count|composition"):
            self._client(module).resolve_query("example", None)

        self.assertTrue(all(handle.closed for handle in module.search_handles))

    def test_distinguishes_request_oserror_from_local_parser_oserror(self):
        request_module = _FakeEntrezModule()
        request_module.search_error = OSError("socket unavailable")
        with self.assertRaises(NcbiTransientError):
            self._client(request_module).resolve_query("example", 1)

        parser_error = OSError("local parser failure")
        parser_module = _FakeEntrezModule()
        parser_module.read_error = parser_error
        with self.assertRaises(OSError) as raised:
            self._client(parser_module).resolve_query("example", 1)

        self.assertIs(raised.exception, parser_error)
        self.assertTrue(parser_module.search_handles[0].closed)

    def test_rejects_incomplete_fetch_composition(self):
        module = _FakeEntrezModule()
        module.fetch_text = self._genbank("NC_1.2")

        with self.assertRaisesRegex(ValueError, "exactly one|count|composition"):
            self._client(module).fetch_records(("101", "102"), identifier_kind="uid")

        self.assertTrue(module.fetch_handles[0].closed)

    def test_sanitizes_transient_network_failures_without_raw_exception_chaining(self):
        failures = (
            HTTPError("https://private-api-key@example.test", 429, "secret", {}, None),
            HTTPError("https://private-api-key@example.test", 503, "secret", {}, None),
            URLError("https://private-api-key@example.test"),
            TimeoutError("private-api-key"),
            ConnectionError("private-api-key"),
        )
        for failure in failures:
            with self.subTest(failure=type(failure).__name__):
                module = _FakeEntrezModule()
                module.search_error = failure
                with self.assertRaises(NcbiTransientError) as raised:
                    self._client(module).resolve_query("example", 1)
                self.assertNotIn("private-api-key", str(raised.exception))
                self.assertIsNone(raised.exception.__cause__)
                self.assertTrue(raised.exception.__suppress_context__)

    def test_distinguishes_and_sanitizes_nonretryable_http_failures(self):
        module = _FakeEntrezModule()
        module.fetch_error = HTTPError(
            "https://private-api-key@example.test", 404, "secret", {}, None
        )

        with self.assertRaisesRegex(
            ncbi.NcbiRequestError, "NCBI request failed with HTTP status 404"
        ) as raised:
            self._client(module).fetch_records(("101",), identifier_kind="uid")

        self.assertNotIn("private-api-key", str(raised.exception))
        self.assertIsNone(raised.exception.__cause__)
        self.assertTrue(raised.exception.__suppress_context__)


class FrozenDatasetTests(unittest.TestCase):
    def _create_dataset(
        self,
        directory: Path,
        *,
        manifest_changes: dict[str, object] | None = None,
        write_records: bool = True,
    ) -> tuple[Path, dict[str, object]]:
        records_path = directory / "records.gb"
        if write_records:
            records = (
                SeqRecord(Seq("ATGCATGC"), id="NC_1.2", description="first record"),
                SeqRecord(Seq("GCTAGCTA"), id="NC_2.3", description="second record"),
            )
            for record in records:
                record.annotations["molecule_type"] = "DNA"
            with records_path.open("w", encoding="utf-8") as handle:
                SeqIO.write(records, handle, "genbank")

        record_bytes = records_path.read_bytes() if write_records else b""
        manifest: dict[str, object] = {
            "schema_version": 1,
            "status": "COMPLETE",
            "source": {
                "mode": "accessions",
                "database": "nuccore",
                "requested_accessions": ["NC_1.2", "NC_2.3"],
            },
            "batch_size": 100,
            "retries": 3,
            "resolved_entries": [
                {
                    "requested_accession": "NC_1.2",
                    "uid": None,
                    "accession": "NC_1",
                    "accession_version": "NC_1.2",
                },
                {
                    "requested_accession": "NC_2.3",
                    "uid": None,
                    "accession": "NC_2",
                    "accession_version": "NC_2.3",
                },
            ],
            "completed_batches": [],
            "consolidated": {
                "filename": "records.gb",
                "record_count": 2,
                "byte_size": len(record_bytes),
                "sha256": hashlib.sha256(record_bytes).hexdigest(),
            },
            "created_at": "2026-08-21T00:00:00+00:00",
            "updated_at": "2026-08-21T00:00:00+00:00",
            "tool": "geison-qpcr",
        }
        manifest.update(manifest_changes or {})
        (directory / "dataset_manifest.json").write_text(
            json.dumps(manifest), encoding="utf-8"
        )
        return records_path, manifest

    def test_validates_complete_dataset_and_returns_its_artifacts(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dataset_dir = Path(tmpdir)
            self._create_dataset(dataset_dir)

            dataset = validate_frozen_dataset(dataset_dir)

            self.assertEqual(dataset.records_path, dataset_dir / "records.gb")
            self.assertEqual(dataset.manifest_path, dataset_dir / "dataset_manifest.json")
            with dataset.records_path.open(encoding="utf-8") as handle:
                self.assertEqual(
                    tuple(record.id for record in SeqIO.parse(handle, "genbank")),
                    ("NC_1.2", "NC_2.3"),
                )

    def test_rejects_missing_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "dataset_manifest.json.*missing"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_partial_manifest(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._create_dataset(Path(tmpdir), manifest_changes={"status": "PARTIAL"})

            with self.assertRaisesRegex(ValueError, "status.*COMPLETE"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_unsupported_schema_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._create_dataset(Path(tmpdir), manifest_changes={"schema_version": 2})

            with self.assertRaisesRegex(ValueError, "schema_version.*unsupported"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_floating_point_schema_version(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._create_dataset(Path(tmpdir), manifest_changes={"schema_version": 1.0})

            with self.assertRaisesRegex(ValueError, "schema_version.*unsupported"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_missing_consolidated_records_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._create_dataset(Path(tmpdir), write_records=False)

            with self.assertRaisesRegex(ValueError, "records.gb.*missing"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_wrong_consolidated_byte_size(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self._create_dataset(
                Path(tmpdir),
                manifest_changes={
                    "consolidated": {
                        "filename": "records.gb",
                        "record_count": 2,
                        "byte_size": 1,
                        "sha256": "0" * 64,
                    }
                },
            )

            with self.assertRaisesRegex(ValueError, "consolidated.byte_size"):
                validate_frozen_dataset(tmpdir)

    def test_rejects_wrong_consolidated_sha256(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            records_path, manifest = self._create_dataset(Path(tmpdir))
            consolidated = dict(manifest["consolidated"])
            consolidated["sha256"] = "0" * 64
            manifest["consolidated"] = consolidated
            (Path(tmpdir) / "dataset_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "consolidated.sha256"):
                validate_frozen_dataset(tmpdir)

            self.assertTrue(records_path.exists())

    def test_rejects_wrong_consolidated_record_count(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            records_path, manifest = self._create_dataset(Path(tmpdir))
            consolidated = dict(manifest["consolidated"])
            consolidated["record_count"] = 1
            manifest["consolidated"] = consolidated
            (Path(tmpdir) / "dataset_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "consolidated.record_count"):
                validate_frozen_dataset(tmpdir)

            self.assertTrue(records_path.exists())

    def test_rejects_accession_version_order_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            _, manifest = self._create_dataset(Path(tmpdir))
            entries = list(manifest["resolved_entries"])
            manifest["resolved_entries"] = list(reversed(entries))
            (Path(tmpdir) / "dataset_manifest.json").write_text(
                json.dumps(manifest), encoding="utf-8"
            )

            with self.assertRaisesRegex(ValueError, "resolved_entries.*accession_version"):
                validate_frozen_dataset(tmpdir)


if __name__ == "__main__":
    unittest.main()
