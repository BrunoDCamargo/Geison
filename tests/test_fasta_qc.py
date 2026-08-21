import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.local_input import load_fasta
from qpcr_pipeline.qc import QCStatus, evaluate_sequences


class FastaQcTests(unittest.TestCase):
    def test_invalid_nucleotide_is_rejected_from_validated_sets(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">invalid\n"
                "ACGTXCGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records)

        self.assertEqual(result.records[0].sequence_id, "accepted")
        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[0].reason_codes, ())

        self.assertEqual(result.records[1].sequence_id, "invalid")
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("INVALID_NUCLEOTIDE",))

        self.assertEqual(result.target_sequence_set.sequence_ids, ("accepted",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("accepted",))

    def test_sequence_shorter_than_configured_minimum_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = Path(tmpdir) / "target.fasta"
            fasta_path.write_text(
                ">accepted\n"
                "ACGTACGT\n"
                ">too-short\n"
                "ACGT\n",
                encoding="utf-8",
            )

            records = load_fasta(fasta_path)
            result = evaluate_sequences(records, min_length=5)

        self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
        self.assertEqual(result.records[1].status, QCStatus.REJECTED)
        self.assertEqual(result.records[1].reason_codes, ("TOO_SHORT",))
        self.assertEqual(result.target_sequence_set.sequence_ids, ("accepted",))
        self.assertEqual(result.evaluation_set.sequence_ids, ("accepted",))


if __name__ == "__main__":
    unittest.main()
