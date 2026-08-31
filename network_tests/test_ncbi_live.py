"""Opt-in live NCBI smoke test, excluded from the standard test suite."""

import json
import os
import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO

from qpcr_pipeline.config import NcbiInputConfig
from qpcr_pipeline.ncbi import acquire_ncbi_dataset, validate_frozen_dataset


RUN_NETWORK_TESTS = (
    os.environ.get("GEISON_RUN_NETWORK_TESTS") == "1"
    and bool(os.environ.get("NCBI_EMAIL", "").strip())
)


@unittest.skipUnless(
    RUN_NETWORK_TESTS,
    "set GEISON_RUN_NETWORK_TESTS=1 and a nonblank NCBI_EMAIL to run live NCBI tests",
)
class NcbiLiveTests(unittest.TestCase):
    def test_acquires_and_validates_stable_accession(self):
        accession = "NC_001416.1"

        with tempfile.TemporaryDirectory() as tmpdir:
            dataset_dir = Path(tmpdir) / "dataset"
            acquired = acquire_ncbi_dataset(
                NcbiInputConfig(accessions=(accession,), batch_size=1), dataset_dir
            )
            manifest = json.loads(acquired.manifest_path.read_text(encoding="utf-8"))
            with acquired.records_path.open(encoding="utf-8") as handle:
                records = tuple(SeqIO.parse(handle, "genbank"))
            validated = validate_frozen_dataset(dataset_dir)

        self.assertEqual(manifest["status"], "COMPLETE")
        self.assertEqual([record.id for record in records], [accession])
        self.assertTrue(records[0].seq)
        self.assertEqual(validated.records_path, dataset_dir / "records.gb")


if __name__ == "__main__":
    unittest.main()
