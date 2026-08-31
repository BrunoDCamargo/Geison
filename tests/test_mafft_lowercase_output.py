import tempfile
from pathlib import Path

from qpcr_pipeline.alignment import align_discovery
from qpcr_pipeline.config import AlignmentConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


class LowercaseMafftRunner:
    def run(self, input_path, output_path, config):
        Path(output_path).write_text(
            ">geison-00000000\nacgt\n>geison-00000001\nacga\n",
            encoding="utf-8",
        )


def test_mafft_lowercase_output_is_normalized_before_validation_and_publication():
    records = (
        LocalSequenceRecord("ref", "ACGT"),
        LocalSequenceRecord("other", "ACGA"),
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        result = align_discovery(
            records,
            DiscoverySet(("ref", "other")),
            AlignmentConfig(enabled=True, reference_id="ref"),
            Path(tmpdir),
            runner=LowercaseMafftRunner(),
        )

        assert result.status == "COMPLETE"
        assert [item.aligned_sequence for item in result.sequences] == ["ACGT", "ACGA"]
        published = result.alignment_fasta_path.read_text(encoding="utf-8")
        assert "acgt" not in published
        assert "acga" not in published
        assert "ACGT" in published
        assert "ACGA" in published
