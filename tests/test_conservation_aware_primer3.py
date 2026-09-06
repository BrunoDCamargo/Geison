from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace
import tempfile

from qpcr_pipeline.config import PrimerDesignConfig
from qpcr_pipeline.conservation import PositionConservation
from qpcr_pipeline.primer_design import CandidateRegion, design_primers


class LiteralRunner:
    def __init__(self, response: str) -> None:
        self.response = response
        self.inputs: list[str] = []

    def run(self, input_text: str) -> str:
        self.inputs.append(input_text)
        return self.response


def _position(index: int, *, conservation: float = 1.0) -> PositionConservation:
    return PositionConservation(
        alignment_position=index,
        reference_position=index,
        reference_base="A",
        depth=50,
        coverage=1.0,
        frequency_a=conservation,
        frequency_c=1.0 - conservation,
        frequency_g=0.0,
        frequency_t=0.0,
        gap_frequency=0.0,
        major_allele_frequency=conservation,
        entropy_bits=0.0,
        major_consensus="A",
        iupac_consensus="A",
    )


def _candidate() -> CandidateRegion:
    return CandidateRegion(
        region_id="contrast-region-001",
        rank=1,
        reference_start=501,
        reference_end=800,
        peak_start=601,
        peak_end=700,
        position_count=300,
        usable_length=300,
        usable_fraction=1.0,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def _anchored_output() -> str:
    return (
        "SEQUENCE_ID=contrast-region-001\n"
        "PRIMER_LEFT_NUM_RETURNED=1\n"
        "PRIMER_INTERNAL_NUM_RETURNED=1\n"
        "PRIMER_RIGHT_NUM_RETURNED=1\n"
        "PRIMER_PAIR_NUM_RETURNED=1\n"
        "PRIMER_LEFT_0=579,20\n"
        "PRIMER_LEFT_0_SEQUENCE=AAAAAAAAAAAAAAAAAAAA\n"
        "PRIMER_LEFT_0_TM=60.0\n"
        "PRIMER_LEFT_0_GC_PERCENT=50.0\n"
        "PRIMER_INTERNAL_0=619,25\n"
        "PRIMER_INTERNAL_0_SEQUENCE=CCCCCCCCCCCCCCCCCCCCCCCCC\n"
        "PRIMER_INTERNAL_0_TM=70.0\n"
        "PRIMER_INTERNAL_0_GC_PERCENT=50.0\n"
        "PRIMER_RIGHT_0=719,20\n"
        "PRIMER_RIGHT_0_SEQUENCE=GGGGGGGGGGGGGGGGGGGG\n"
        "PRIMER_RIGHT_0_TM=60.0\n"
        "PRIMER_RIGHT_0_GC_PERCENT=50.0\n"
        "PRIMER_PAIR_0_PRODUCT_SIZE=141\n"
        "=\n"
    )


def test_contrastive_design_excludes_low_conservation_binding_positions() -> None:
    positions = tuple(
        _position(index, conservation=0.50 if 560 <= index <= 562 else 1.0)
        for index in range(1, 1201)
    )
    conservation = SimpleNamespace(
        reference_id="ref",
        major_consensus="A" * 1200,
        positions=positions,
    )
    contrastive = SimpleNamespace(
        status="COMPLETE",
        candidates=(SimpleNamespace(region=_candidate()),),
    )
    runner = LiteralRunner(_anchored_output())

    with tempfile.TemporaryDirectory() as temporary:
        design_primers(
            conservation,
            PrimerDesignConfig(enabled=True, assays_per_region=1),
            Path(temporary),
            contrastive=contrastive,
            runner=runner,
        )

    assert len(runner.inputs) == 1
    text = runner.inputs[0]
    assert "SEQUENCE_EXCLUDED_REGION=559,3\n" in text
    assert "SEQUENCE_INTERNAL_EXCLUDED_REGION=559,3\n" in text
