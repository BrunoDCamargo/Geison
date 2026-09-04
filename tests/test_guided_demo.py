from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import yaml

from qpcr_pipeline.config import (
    ContrastiveConservationConfig,
    OffTargetConfig,
    PrimerDesignConfig,
)
from qpcr_pipeline.conservation import (
    ConservationResult,
    PositionConservation,
    WindowConservation,
)
from qpcr_pipeline.contrastive_conservation import analyze_contrastive_conservation
from qpcr_pipeline.contrastive_similarity import BiopythonLocalSimilarityEngine
from qpcr_pipeline.local_input import load_fasta
from qpcr_pipeline.panel_manifest import (
    approve_panel_proposal,
    load_approved_panel_manifest,
)


GENERATOR = Path("examples/guided_demo/generate_demo_data.py")
EXPECTED_FILES = (
    "target.fasta",
    "challenge-related-a.fasta",
    "challenge-related-b.fasta",
    "challenge-context-a.fasta",
    "panel-proposal.yaml",
    "config-proposal.yaml",
    "config-approved-template.yaml",
)
SHARED_START = 201
SHARED_END = 300
DISCRIMINANT_START = 601
DISCRIMINANT_END = 700


def _generate(directory: Path) -> None:
    completed = subprocess.run(
        [sys.executable, str(GENERATOR), str(directory)],
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr or completed.stdout


def _read_yaml(path: Path):
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def _position(index: int, base: str) -> PositionConservation:
    return PositionConservation(
        alignment_position=index,
        reference_position=index,
        reference_base=base,
        depth=4,
        coverage=1.0,
        frequency_a=1.0 if base == "A" else 0.0,
        frequency_c=1.0 if base == "C" else 0.0,
        frequency_g=1.0 if base == "G" else 0.0,
        frequency_t=1.0 if base == "T" else 0.0,
        gap_frequency=0.0,
        major_allele_frequency=1.0,
        entropy_bits=0.0,
        major_consensus=base,
        iupac_consensus=base,
    )


def _window(start: int, end: int) -> WindowConservation:
    return WindowConservation(
        reference_start=start,
        reference_end=end,
        position_count=end - start + 1,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def _conservation(reference: str, report_path: Path) -> ConservationResult:
    return ConservationResult(
        status="COMPLETE",
        reference_id="synthetic-target-reference",
        positions=tuple(_position(index, base) for index, base in enumerate(reference, 1)),
        windows=(
            _window(SHARED_START, SHARED_END),
            _window(DISCRIMINANT_START, DISCRIMINANT_END),
        ),
        annotations=(),
        major_consensus=reference,
        iupac_consensus=reference,
        position_metrics_path=None,
        window_metrics_path=None,
        major_consensus_path=None,
        iupac_consensus_path=None,
        html_report_path=None,
        report_path=report_path,
    )


def test_generator_is_deterministic_and_emits_guided_inputs(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    _generate(first)
    _generate(second)

    assert tuple(sorted(path.name for path in first.iterdir())) == tuple(sorted(EXPECTED_FILES))
    for name in EXPECTED_FILES:
        assert (first / name).read_bytes() == (second / name).read_bytes()

    targets = load_fasta(first / "target.fasta")
    assert len(targets) >= 4
    assert all(len(record.sequence) == 1200 for record in targets)
    reference = targets[0].sequence
    windows = [reference[index : index + 100] for index in range(len(reference) - 99)]
    assert len(windows) == len(set(windows))

    panel = _read_yaml(first / "panel-proposal.yaml")
    non_targets = panel["definition"]["non_targets"]
    criticalities = {item["criticality"] for item in non_targets}
    assert {"CRITICAL", "IMPORTANT"} <= criticalities
    challenge_names = [
        item["name"] for item in non_targets if "CHALLENGE" in item["dataset_roles"]
    ]

    proposal = _read_yaml(first / "config-proposal.yaml")
    approved = _read_yaml(first / "config-approved-template.yaml")
    assert proposal["contrastive_conservation"]["enabled"] is False
    assert "proposal" in proposal["panel"]
    assert approved["contrastive_conservation"]["enabled"] is True
    assert approved["panel"]["frozen_manifest"] == "__APPROVED_PANEL__"

    for stage in ("alignment", "conservation", "primer_design", "inclusivity", "specificity", "ranking"):
        assert approved[stage]["enabled"] is True

    off_target_names = [item["name"] for item in approved["off_targets"]]
    assert off_target_names == challenge_names

    rendered = (first / "config-approved-template.yaml").read_text(encoding="utf-8").lower()
    assert "assay_sequence" not in rendered
    assert "real assay" not in rendered


def test_real_similarity_distinguishes_shared_and_discriminant_and_ranks_discriminant_first(tmp_path):
    demo = tmp_path / "demo"
    _generate(demo)

    target_records = load_fasta(demo / "target.fasta")
    reference = target_records[0].sequence
    shared_query = reference[SHARED_START - 1 : SHARED_END]
    discriminant_query = reference[DISCRIMINANT_START - 1 : DISCRIMINANT_END]

    challenge_paths = (
        demo / "challenge-related-a.fasta",
        demo / "challenge-related-b.fasta",
        demo / "challenge-context-a.fasta",
    )
    engine = BiopythonLocalSimilarityEngine()
    for path in challenge_paths:
        records = load_fasta(path)
        shared = engine.best_match(shared_query, records)
        discriminant = engine.best_match(discriminant_query, records)
        assert shared is not None and discriminant is not None
        assert shared.similarity == 1.0
        assert discriminant.similarity < shared.similarity

    approved_path = demo / "approved-panel.json"
    approve_panel_proposal(demo / "panel-proposal.yaml", approved_path)
    manifest = load_approved_panel_manifest(approved_path)
    off_targets = (
        OffTargetConfig(manifest.definition.non_targets[0].name, fasta=challenge_paths[0]),
        OffTargetConfig(manifest.definition.non_targets[1].name, fasta=challenge_paths[1]),
        OffTargetConfig(manifest.definition.non_targets[2].name, fasta=challenge_paths[2]),
    )
    result = analyze_contrastive_conservation(
        _conservation(reference, demo / "conservation-report.json"),
        manifest,
        off_targets,
        ContrastiveConservationConfig(enabled=True),
        PrimerDesignConfig(
            enabled=True,
            candidate_region_length=200,
            max_candidate_regions=2,
            product_size_max=180,
            min_mean_conservation=0.0,
            min_minimum_conservation=0.0,
            min_mean_coverage=0.0,
            max_mean_gap_frequency=1.0,
            max_mean_entropy_bits=2.0,
            min_usable_fraction=0.0,
        ),
        demo / "analysis",
    )

    by_start = {window.reference_start: window for window in result.windows}
    assert by_start[DISCRIMINANT_START].worst_similarity < by_start[SHARED_START].worst_similarity
    assert result.candidates[0].region.peak_start == DISCRIMINANT_START
    assert result.candidates[0].region.peak_end == DISCRIMINANT_END
