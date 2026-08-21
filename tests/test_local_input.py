import tempfile
import unittest
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import qpcr_pipeline.local_input as local_input
from qpcr_pipeline.qc import QCStatus, evaluate_sequences


class LocalInputTests(unittest.TestCase):
    def test_load_genbank_normalizes_sequence_and_preserves_metadata(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            genbank_path = Path(tmpdir) / "target.gb"
            source = SeqRecord(
                Seq("acgtn"),
                id="GB001",
                name="GENE1",
                description="Synthetic target",
            )
            source.annotations.update(
                {"molecule_type": "DNA", "organism": "Synthetic construct"}
            )
            source.dbxrefs = ["BioProject:PRJ1"]
            source.features = [
                SeqFeature(FeatureLocation(0, 5), type="gene", qualifiers={"gene": ["abc"]})
            ]
            SeqIO.write((source,), genbank_path, "genbank")

            loaded = local_input.load_genbank(genbank_path)

        self.assertEqual(loaded[0].sequence_id, "GB001")
        self.assertEqual(loaded[0].sequence, "ACGTN")
        self.assertEqual(loaded[0].metadata["name"], "GENE1")
        self.assertEqual(loaded[0].metadata["description"], "Synthetic target")
        self.assertEqual(
            loaded[0].metadata["annotations"]["organism"], "Synthetic construct"
        )
        self.assertEqual(loaded[0].metadata["dbxrefs"], ("BioProject:PRJ1",))
        self.assertEqual(loaded[0].metadata["features"][0].type, "gene")

    def test_load_fasta_preserves_available_name_and_description(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(">seq-1 descriptive target\nacgt\n", encoding="utf-8")

            loaded = local_input.load_fasta(fasta_path)

        self.assertEqual(loaded[0].sequence_id, "seq-1")
        self.assertEqual(loaded[0].sequence, "ACGT")
        self.assertEqual(loaded[0].metadata["name"], "seq-1")
        self.assertEqual(loaded[0].metadata["description"], "seq-1 descriptive target")

    def test_missing_optional_genbank_metadata_does_not_reject_sequence(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            genbank_path = Path(tmpdir) / "metadata-light.gb"
            source = SeqRecord(Seq("ACGT"), id="metadata-light")
            source.annotations["molecule_type"] = "DNA"
            SeqIO.write((source,), genbank_path, "genbank")

            result = evaluate_sequences(local_input.load_genbank(genbank_path))

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[0].reason_codes, ())
        self.assertEqual(result.evaluation_set.sequence_ids, ("metadata-light",))

    def test_load_local_sequences_selects_each_supported_format(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            directory = Path(tmpdir)
            fasta_path = directory / "target.fasta"
            genbank_path = directory / "target.gb"
            fasta_path.write_text(">FA001\nacgt\n", encoding="utf-8")
            source = SeqRecord(Seq("tgca"), id="GB002")
            source.annotations["molecule_type"] = "DNA"
            SeqIO.write((source,), genbank_path, "genbank")

            fasta_records = local_input.load_local_sequences(fasta_path, "fasta")
            genbank_records = local_input.load_local_sequences(genbank_path, "genbank")

        self.assertEqual(fasta_records[0].sequence_id, "FA001")
        self.assertEqual(genbank_records[0].sequence_id, "GB002")

    def test_load_local_sequences_rejects_unsupported_format(self):
        with self.assertRaisesRegex(
            ValueError, "Unsupported local sequence format: embl"
        ):
            local_input.load_local_sequences("target.embl", "embl")


if __name__ == "__main__":
    unittest.main()
