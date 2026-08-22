import hashlib
import json
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import NcbiInputConfig
from qpcr_pipeline.ncbi import (
    NcbiFetchedRecord,
    NcbiTransientError,
    acquire_ncbi_dataset,
    validate_frozen_dataset,
)


class FakeNcbiClient:
    def __init__(self, records_by_request, failures=()):
        self.records_by_request = records_by_request
        self.failures = list(failures)
        self.fetch_calls = []
        self.resolve_calls = []

    def resolve_query(self, query, max_records):
        self.resolve_calls.append((query, max_records))
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

    def test_requires_an_explicit_client_until_the_production_adapter_exists(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertRaisesRegex(ValueError, "NcbiClient"):
                acquire_ncbi_dataset(
                    self._config(("NC_000001.11",)), Path(tmpdir), clock=lambda: "now"
                )

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
