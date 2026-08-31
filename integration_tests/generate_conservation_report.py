"""Generate a deterministic report for manual browser verification."""

import sys
from pathlib import Path

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from Bio.SeqFeature import SeqFeature, SimpleLocation

from qpcr_pipeline.alignment import (
    AlignedSequence,
    AlignmentCoordinate,
    AlignmentResult,
)
from qpcr_pipeline.config import ConservationConfig
from qpcr_pipeline.conservation import analyze_conservation
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet


if len(sys.argv) != 2:
    raise SystemExit("usage: generate_conservation_report.py OUTPUT_DIRECTORY")
output_dir = Path(sys.argv[1]).resolve()

reference = "ACGT" * 100
variant_one = list(reference)
variant_two = list(reference)
for position in range(140, 161):
    variant_one[position] = "A" if reference[position] != "A" else "C"
for position in range(270, 291):
    variant_two[position] = "T" if reference[position] != "T" else "G"

records = (
    LocalSequenceRecord(
        "ref",
        reference,
        metadata={
            "features": (
                SeqFeature(
                    SimpleLocation(20, 90, strand=1),
                    type="gene",
                    qualifiers={"gene": ["alpha"]},
                ),
                SeqFeature(
                    SimpleLocation(120, 190, strand=1),
                    type="CDS",
                    qualifiers={"product": ["beta protein"]},
                ),
                SeqFeature(
                    SimpleLocation(240, 320, strand=-1),
                    type="gene",
                    qualifiers={"gene": ["gamma"]},
                ),
                SeqFeature(
                    SimpleLocation(340, 390, strand=1),
                    type="regulatory",
                    qualifiers={"note": ["terminal region"]},
                ),
            )
        },
    ),
    LocalSequenceRecord("variant-1", "".join(variant_one)),
    LocalSequenceRecord("variant-2", "".join(variant_two)),
)
discovery = DiscoverySet(tuple(record.sequence_id for record in records))
alignment = AlignmentResult(
    status="COMPLETE",
    discovery_set=discovery,
    reference_id="ref",
    reference_mode="explicit",
    sequences=tuple(
        AlignedSequence(record.sequence_id, record.sequence, "forward")
        for record in records
    ),
    coordinates=tuple(
        AlignmentCoordinate(position, position, base)
        for position, base in enumerate(reference, start=1)
    ),
    alignment_fasta_path=None,
    coordinate_map_path=None,
    report_path=output_dir / "alignment" / "alignment_report.json",
)
result = analyze_conservation(
    records,
    alignment,
    ConservationConfig(enabled=True, window_size=50, step_size=10),
    output_dir,
    target_name="Geison browser fixture",
)
print(result.html_report_path.resolve())
