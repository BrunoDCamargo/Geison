import json
import shutil
import tempfile
import unittest
from pathlib import Path

from Bio.Seq import Seq

from qpcr_pipeline.alignment import align_discovery
from qpcr_pipeline.config import AlignmentConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


@unittest.skipUnless(shutil.which("mafft"), "mafft is not installed")
class MafftIntegrationTests(unittest.TestCase):
    def test_mafft_aligns_reverse_complement_and_publishes_traceable_artifacts(self):
        reference = "ACGTTGCAAGTCCTGATCGATGCTAGCTTACGATCGGATCCGTTACGATGCCATGATCGTACGATGCTAGCATCGTACGGA"
        reverse_complement = str(Seq(reference).reverse_complement())
        forward_variant = reference[:40] + "A" + reference[41:]
        records = (
            LocalSequenceRecord("seq-1", reference),
            LocalSequenceRecord("seq-2", reverse_complement),
            LocalSequenceRecord("seq-3", forward_variant),
        )
        discovery = DiscoverySet(("seq-1", "seq-2", "seq-3"))

        with tempfile.TemporaryDirectory() as tmpdir:
            result = align_discovery(
                records,
                discovery,
                AlignmentConfig(enabled=True, reference_id="seq-1", threads=2),
                Path(tmpdir),
            )
            report = json.loads(result.report_path.read_text(encoding="utf-8"))

            self.assertTrue(result.alignment_fasta_path.is_file())
            self.assertTrue(result.coordinate_map_path.is_file())

        self.assertEqual(result.reference_id, "seq-1")
        self.assertEqual(
            [sequence.orientation for sequence in result.sequences],
            ["forward", "reverse_complemented", "forward"],
        )
        self.assertEqual(result.coordinates[0].alignment_position, 1)
        self.assertEqual(result.coordinates[-1].reference_position, len(reference))
        self.assertEqual([sequence.sequence_id for sequence in result.sequences], ["seq-1", "seq-2", "seq-3"])
        self.assertEqual(report["counts"], {
            "discovery": 3,
            "alignment_length": len(reference),
            "reverse_complemented": 1,
        })


if __name__ == "__main__":
    unittest.main()
