from pathlib import Path

import pytest
from Bio.SeqFeature import SeqFeature, SimpleLocation

from qpcr_pipeline.alignment import AlignedSequence, AlignmentCoordinate, AlignmentResult
from qpcr_pipeline.checkpoint_codecs import (
    ALIGNMENT_CODEC,
    CLUSTERING_CODEC,
    CONSERVATION_CODEC,
    INCLUSIVITY_CODEC,
    INPUT_CODEC,
    PRIMER_DESIGN_CODEC,
    QC_CODEC,
    RANKING_CODEC,
    SPECIFICITY_CODEC,
)
from qpcr_pipeline.clustering import ClusterMember, ClusteringResult, SequenceCluster
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    ReferenceAnnotation,
    WindowConservation,
)
from qpcr_pipeline.inclusivity import (
    AssayInclusivity,
    DegeneracyProposal,
    InclusivityResult,
    OligoMatch,
    OligoVariation,
    ProposedOligoCompatibility,
)
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import DiscoverySet, EvaluationSet, TargetSequenceSet
from qpcr_pipeline.primer_design import (
    AssayCandidate,
    CandidateRegion,
    DesignedOligo,
    PrimerDesignResult,
)
from qpcr_pipeline.qc import QCRecord, QCResult, QCStatus
from qpcr_pipeline.ranking import RankingReason, RankingResult, RankedAssay, ScoreComponents
from qpcr_pipeline.specificity import (
    HitRetentionSummary,
    OffTargetHit,
    PlausibleAmplicon,
    SpecificityResult,
)


def _round_trip(codec, value, outdir):
    payload = codec.encode(value, outdir)
    return payload, codec.decode(payload, outdir)


def test_input_records_round_trip_and_preserve_required_feature_metadata(tmp_path):
    feature = SeqFeature(
        SimpleLocation(1, 4, strand=1),
        type="CDS",
        qualifiers={"gene": ["abc"], "product": ["protein"]},
    )
    value = (
        LocalSequenceRecord(
            "s1",
            "ACGT",
            {
                "name": "s1",
                "description": "example",
                "annotations": {"molecule_type": "DNA"},
                "dbxrefs": ("DB:1",),
                "features": (feature,),
            },
        ),
    )
    payload, decoded = _round_trip(INPUT_CODEC, value, tmp_path)
    assert decoded[0].sequence_id == "s1"
    assert decoded[0].sequence == "ACGT"
    decoded_feature = decoded[0].metadata["features"][0]
    assert decoded_feature.type == "CDS"
    assert int(decoded_feature.location.start) == 1
    assert int(decoded_feature.location.end) == 4
    assert decoded_feature.location.strand == 1
    assert decoded_feature.qualifiers["gene"] == ["abc"]
    assert isinstance(payload, list)


def test_qc_round_trip_is_exact(tmp_path):
    value = QCResult(
        records=(QCRecord("s1", QCStatus.ACCEPTED, ()),),
        target_sequence_set=TargetSequenceSet(("s1",)),
        evaluation_set=EvaluationSet(("s1",)),
    )
    _, decoded = _round_trip(QC_CODEC, value, tmp_path)
    assert decoded == value


def test_clustering_round_trip_encodes_generated_paths_relative_to_outdir(tmp_path):
    outdir = tmp_path / "run"
    value = ClusteringResult(
        discovery_set=DiscoverySet(("s1",)),
        clusters=(
            SequenceCluster(
                "cluster-1",
                "s1",
                (ClusterMember("s1", True, None, None),),
            ),
        ),
        discovery_fasta_path=outdir / "discovery_set.fasta",
        report_path=outdir / "clustering_report.json",
        raw_cluster_path=outdir / "clustering" / "cd-hit-est.clstr",
    )
    payload, decoded = _round_trip(CLUSTERING_CODEC, value, outdir)
    assert decoded == value
    assert str(tmp_path) not in repr(payload)


def test_alignment_round_trip_preserves_nested_values(tmp_path):
    outdir = tmp_path / "run"
    value = AlignmentResult(
        status="COMPLETE",
        discovery_set=DiscoverySet(("s1",)),
        reference_id="s1",
        reference_mode="explicit",
        sequences=(AlignedSequence("s1", "ACGT", "forward"),),
        coordinates=(AlignmentCoordinate(1, 1, "A"),),
        alignment_fasta_path=outdir / "alignment" / "discovery_alignment.fasta",
        coordinate_map_path=outdir / "alignment" / "coordinate_map.tsv",
        report_path=outdir / "alignment" / "alignment_report.json",
    )
    _, decoded = _round_trip(ALIGNMENT_CODEC, value, outdir)
    assert decoded == value


def test_conservation_round_trip_preserves_metrics_and_annotations(tmp_path):
    outdir = tmp_path / "run"
    position = PositionConservation(1, 1, "A", 1, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, "A", "A")
    window = WindowConservation(1, 1, 1, 1.0, 1.0, 1.0, 0.0, 0.0)
    value = ConservationResult(
        "COMPLETE",
        "s1",
        (position,),
        (window,),
        (ReferenceAnnotation("CDS", 1, 1, 1, "gene"),),
        "A",
        "A",
        outdir / "conservation" / "position_metrics.tsv",
        outdir / "conservation" / "window_metrics.tsv",
        outdir / "conservation" / "major_consensus.fasta",
        outdir / "conservation" / "iupac_consensus.fasta",
        outdir / "conservation" / "report.html",
        outdir / "conservation" / "conservation_report.json",
    )
    _, decoded = _round_trip(CONSERVATION_CODEC, value, outdir)
    assert decoded == value


def _oligo(sequence, start):
    return DesignedOligo(sequence, start, start + len(sequence) - 1, len(sequence), 60.0, 50.0, 0.1, (("x", "y"),))


def _primer_result(outdir):
    region = CandidateRegion("region-1", 1, 1, 100, 20, 80, 100, 100, 1.0, 0.99, 0.98, 1.0, 0.0, 0.01)
    assay = AssayCandidate(
        "assay-1",
        "region-1",
        0,
        _oligo("ACGT", 1),
        _oligo("CGTA", 10),
        _oligo("GTAC", 20),
        23,
        0.2,
        (("pair", "ok"),),
    )
    return PrimerDesignResult(
        "COMPLETE",
        "s1",
        (region,),
        (assay,),
        outdir / "primer_design" / "candidate_regions.tsv",
        outdir / "primer_design" / "assays.tsv",
        outdir / "primer_design" / "primer3_input.txt",
        outdir / "primer_design" / "primer3_output.txt",
        outdir / "primer_design" / "primer_design_report.json",
    )


def test_primer_design_round_trip_preserves_assays(tmp_path):
    outdir = tmp_path / "run"
    value = _primer_result(outdir)
    _, decoded = _round_trip(PRIMER_DESIGN_CODEC, value, outdir)
    assert decoded == value


def test_inclusivity_round_trip_preserves_nested_evidence(tmp_path):
    outdir = tmp_path / "run"
    match = OligoMatch("assay-1", "s1", "FORWARD", "FORWARD", 1, 1, 4, 1, 4, 0, (), 0, True, False, False, True, True)
    proposed = ProposedOligoCompatibility("FORWARD", "ACGT", True, (), 0, False, False, True)
    assay = AssayInclusivity("assay-1", "s1", "FORWARD", True, 1, 23, 23, match, None, None, True, proposed, None, None, True)
    variation = OligoVariation("assay-1", "FORWARD", 1, "A", ("A",), "R", ("A", "G"), ("s1",), 1, 1.0, False)
    proposal = DegeneracyProposal("assay-1", "FORWARD", "A", "R", "ACCEPTED", "variation", 1, 2, (1,), 1, 0, 0.0, 1, 1.0)
    value = InclusivityResult(
        "COMPLETE",
        ("s1",),
        (match,),
        (assay,),
        (variation,),
        (proposal,),
        outdir / "inclusivity" / "oligo_matches.tsv",
        outdir / "inclusivity" / "assay_inclusivity.tsv",
        outdir / "inclusivity" / "oligo_variations.tsv",
        outdir / "inclusivity" / "degeneracy_proposals.tsv",
        outdir / "inclusivity" / "inclusivity_report.json",
    )
    _, decoded = _round_trip(INCLUSIVITY_CODEC, value, outdir)
    assert decoded == value


def test_specificity_round_trip_preserves_hits_amplicons_and_retention(tmp_path):
    outdir = tmp_path / "run"
    hit = OffTargetHit("db", "assay-1", "off1", "FORWARD", "FORWARD", 1, 1, 4, (), 0, True, False, True)
    retention = HitRetentionSummary("db", "assay-1", "FORWARD", 1, 1, False)
    amplicon = PlausibleAmplicon("db", "assay-1", "off1", "FORWARD", 1, 23, 23, 1, 4, 20, 23, ((10, 13),), 1, 1, (1,), True, True)
    value = SpecificityResult(
        "COMPLETE",
        ("db",),
        1,
        1,
        (hit,),
        (amplicon,),
        (retention,),
        outdir / "specificity" / "off_target_hits.tsv",
        outdir / "specificity" / "plausible_amplicons.tsv",
        outdir / "specificity" / "specificity_report.json",
    )
    _, decoded = _round_trip(SPECIFICITY_CODEC, value, outdir)
    assert decoded == value


def test_ranking_round_trip_preserves_reason_evidence_order(tmp_path):
    outdir = tmp_path / "run"
    reason = RankingReason("EVIDENCE_INCOMPLETE", "REVIEW", "ranking", "missing", (("b", 2), ("a", 1)))
    ranked = RankedAssay(1, "assay-1", "region-1", 0, "REVIEW", 90.0, "COMPLETE", ScoreComponents(1.0, 0.8, 0.9, 0.9, 0.9), (reason,), 1, 1, 2, 0, 0, 0.2)
    value = RankingResult(
        "COMPLETE",
        (ranked,),
        outdir / "ranking" / "assay_ranking.tsv",
        outdir / "ranking" / "ranking_report.json",
        outdir / "report.html",
    )
    _, decoded = _round_trip(RANKING_CODEC, value, outdir)
    assert decoded == value
    assert decoded.assays[0].reasons[0].evidence == (("b", 2), ("a", 1))


def test_codec_rejects_unknown_top_level_field(tmp_path):
    value = QCResult((), TargetSequenceSet(()), EvaluationSet(()))
    payload = QC_CODEC.encode(value, tmp_path)
    payload["unexpected"] = True
    with pytest.raises(ValueError, match="fields"):
        QC_CODEC.decode(payload, tmp_path)


def test_codec_rejects_generated_path_that_escapes_outdir(tmp_path):
    outdir = tmp_path / "run"
    value = AlignmentResult("SKIPPED", DiscoverySet(()), None, None, (), (), None, None, outdir / "alignment" / "alignment_report.json")
    payload = ALIGNMENT_CODEC.encode(value, outdir)
    payload["report_path"] = "../escape.json"
    with pytest.raises(ValueError, match="escapes"):
        ALIGNMENT_CODEC.decode(payload, outdir)
