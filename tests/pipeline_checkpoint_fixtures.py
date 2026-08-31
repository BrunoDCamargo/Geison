from dataclasses import replace
from pathlib import Path

from qpcr_pipeline.alignment import AlignmentResult
from qpcr_pipeline.clustering import ClusteringResult
from qpcr_pipeline.conservation import ConservationResult
from qpcr_pipeline.inclusivity import InclusivityResult
from qpcr_pipeline.models import DiscoverySet
from qpcr_pipeline.primer_design import PrimerDesignResult
from qpcr_pipeline.specificity import SpecificityResult


def _write(path: Path, text: str = "{}\n") -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")
    return path


def checkpoint_clustering(output: Path) -> ClusteringResult:
    discovery = _write(output / "discovery_set.fasta", ">seq-1\nACGT\n")
    report = _write(output / "clustering_report.json")
    return ClusteringResult(
        discovery_set=DiscoverySet(("seq-1",)),
        clusters=(),
        discovery_fasta_path=discovery,
        report_path=report,
        raw_cluster_path=None,
    )


def checkpoint_alignment(output: Path) -> AlignmentResult:
    report = _write(output / "alignment" / "alignment_report.json")
    return AlignmentResult(
        status="COMPLETE",
        discovery_set=DiscoverySet(("seq-1",)),
        reference_id="seq-1",
        reference_mode="automatic",
        sequences=(),
        coordinates=(),
        alignment_fasta_path=None,
        coordinate_map_path=None,
        report_path=report,
    )


def checkpoint_conservation(output: Path) -> ConservationResult:
    report = _write(output / "conservation" / "conservation_report.json")
    return ConservationResult(
        status="COMPLETE",
        reference_id="seq-1",
        positions=(),
        windows=(),
        annotations=(),
        major_consensus="",
        iupac_consensus="",
        position_metrics_path=None,
        window_metrics_path=None,
        major_consensus_path=None,
        iupac_consensus_path=None,
        html_report_path=None,
        report_path=report,
    )


def checkpoint_primer(output: Path, result: PrimerDesignResult) -> PrimerDesignResult:
    report = _write(output / "primer_design" / "primer_design_report.json")
    return replace(
        result,
        candidate_regions_path=None,
        assays_path=None,
        primer3_input_path=None,
        primer3_output_path=None,
        report_path=report,
    )


def checkpoint_inclusivity(output: Path, result: InclusivityResult) -> InclusivityResult:
    report = _write(output / "inclusivity" / "inclusivity_report.json")
    return replace(
        result,
        oligo_matches_path=None,
        assay_inclusivity_path=None,
        oligo_variations_path=None,
        degeneracy_proposals_path=None,
        report_path=report,
    )


def skipped_inclusivity(output: Path) -> InclusivityResult:
    return checkpoint_inclusivity(
        output,
        InclusivityResult(
            status="SKIPPED",
            evaluation_sequence_ids=(),
            oligo_matches=(),
            assay_results=(),
            variations=(),
            proposals=(),
            oligo_matches_path=None,
            assay_inclusivity_path=None,
            oligo_variations_path=None,
            degeneracy_proposals_path=None,
            report_path=output / "inclusivity" / "inclusivity_report.json",
        ),
    )


def checkpoint_specificity(output: Path, result: SpecificityResult) -> SpecificityResult:
    report = _write(output / "specificity" / "specificity_report.json")
    return replace(
        result,
        off_target_hits_path=None,
        plausible_amplicons_path=None,
        report_path=report,
    )
