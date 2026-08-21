import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.ncbi import validate_frozen_dataset


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
