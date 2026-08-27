"""Minimal pipeline orchestration used by the first end-to-end tracer bullet."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path

from qpcr_pipeline.alignment import MafftRunner, align_discovery
from qpcr_pipeline.clustering import CdHitRunner, cluster_sequences
from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.conservation import analyze_conservation
from qpcr_pipeline.inclusivity import evaluate_inclusivity
from qpcr_pipeline.local_input import load_genbank, load_local_sequences
from qpcr_pipeline.ncbi import NcbiClient, acquire_ncbi_dataset, validate_frozen_dataset
from qpcr_pipeline.primer3 import Primer3Runner
from qpcr_pipeline.primer_design import design_primers
from qpcr_pipeline.qc import evaluate_sequences
from qpcr_pipeline.specificity import evaluate_specificity


@dataclass(frozen=True, slots=True)
class RunSummary:
    status: str
    target_name: str
    sequence_count: int
    sequence_ids: list[str]


def run_pipeline(
    config: PipelineConfig,
    outdir: str | Path,
    *,
    ncbi_client: NcbiClient | None = None,
    cdhit_runner: CdHitRunner | None = None,
    mafft_runner: MafftRunner | None = None,
    primer3_runner: Primer3Runner | None = None,
) -> RunSummary:
    selected_input = config.selected_input
    output_dir = Path(outdir)
    if (
        isinstance(selected_input, NcbiInputConfig)
        and selected_input.frozen_dataset is not None
    ):
        _reject_output_inside_frozen_dataset(
            output_dir, selected_input.frozen_dataset
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    if isinstance(selected_input, NcbiInputConfig):
        if selected_input.frozen_dataset is not None:
            acquired = validate_frozen_dataset(selected_input.frozen_dataset)
        else:
            acquired = acquire_ncbi_dataset(
                selected_input,
                output_dir / "ncbi_dataset",
                client=ncbi_client,
            )
        records = load_genbank(acquired.records_path)
        _copy_effective_manifest(
            acquired.manifest_path, output_dir / "ncbi_dataset_manifest.json"
        )
    else:
        input_path, input_format = selected_input
        records = load_local_sequences(input_path, input_format)
    result = evaluate_sequences(
        records,
        min_length=config.qc.min_length,
        max_ambiguous_fraction=config.qc.max_ambiguous_fraction,
        expected_length=config.qc.expected_length,
        length_tolerance_fraction=config.qc.length_tolerance_fraction,
    )
    approved_ids = result.evaluation_set.sequence_ids
    approved_id_set = set(approved_ids)
    approved_records = tuple(
        record for record in records if record.sequence_id in approved_id_set
    )
    clustering = cluster_sequences(
        approved_records,
        result.evaluation_set,
        config.clustering,
        output_dir,
        runner=cdhit_runner,
    )
    discovery_ids = clustering.discovery_set.sequence_ids
    discovery_id_set = set(discovery_ids)
    discovery_records = tuple(
        record for record in approved_records if record.sequence_id in discovery_id_set
    )
    alignment = align_discovery(
        discovery_records,
        clustering.discovery_set,
        config.alignment,
        output_dir,
        runner=mafft_runner,
    )
    conservation = analyze_conservation(
        discovery_records,
        alignment,
        config.conservation,
        output_dir,
        target_name=config.target_name,
    )
    primer_design = design_primers(
        conservation,
        config.primer_design,
        output_dir,
        runner=primer3_runner,
    )
    inclusivity = evaluate_inclusivity(
        approved_records,
        result.evaluation_set,
        primer_design,
        config.inclusivity,
        output_dir,
    )
    specificity = evaluate_specificity(
        primer_design,
        config.off_targets,
        config.specificity,
        output_dir,
    )

    summary = RunSummary(
        status="COMPLETED",
        target_name=config.target_name,
        sequence_count=len(approved_ids),
        sequence_ids=list(approved_ids),
    )

    qc_report = {
        "records": [
            {
                "sequence_id": record.sequence_id,
                "status": record.status.value,
                "reason_codes": list(record.reason_codes),
            }
            for record in result.records
        ],
        "target_sequence_set": {"sequence_ids": list(result.target_sequence_set.sequence_ids)},
        "evaluation_set": {"sequence_ids": list(approved_ids)},
        "discovery_set": {"sequence_ids": list(clustering.discovery_set.sequence_ids)},
        "alignment": {
            "status": alignment.status,
            "reference_id": alignment.reference_id,
            "reference_mode": alignment.reference_mode,
        },
        "conservation": {
            "status": conservation.status,
            "reference_id": conservation.reference_id,
            "position_count": len(conservation.positions),
            "window_count": len(conservation.windows),
        },
        "primer_design": {
            "status": primer_design.status,
            "reference_id": primer_design.reference_id,
            "candidate_region_count": len(primer_design.candidates),
            "assay_count": len(primer_design.assays),
        },
        "inclusivity": {
            "status": inclusivity.status,
            "evaluation_sequence_count": len(inclusivity.evaluation_sequence_ids),
            "assay_count": (
                len(primer_design.assays)
                if inclusivity.status == "COMPLETE"
                else 0
            ),
            "assay_evaluation_count": len(inclusivity.assay_results),
            "original_compatible_count": sum(
                assay.original_compatible for assay in inclusivity.assay_results
            ),
            "proposed_compatible_count": sum(
                assay.proposed_compatible for assay in inclusivity.assay_results
            ),
        },
        "specificity": {
            "status": specificity.status,
            "dataset_count": len(specificity.dataset_names),
            "sequence_count": specificity.sequence_count,
            "assay_count": specificity.assay_count,
            "retained_hit_count": len(specificity.hits),
            "plausible_amplicon_count": len(specificity.amplicons),
            "detectable_off_target_count": sum(
                amplicon.detectable_off_target for amplicon in specificity.amplicons
            ),
        },
    }

    _write_json_atomic(
        output_dir / "run_summary.json", asdict(summary)
    )
    _write_json_atomic(
        output_dir / "qc_report.json", qc_report
    )
    return summary


def _reject_output_inside_frozen_dataset(output_dir: Path, frozen_dir: Path) -> None:
    resolved_output = output_dir.resolve()
    resolved_frozen = frozen_dir.resolve()
    try:
        resolved_output.relative_to(resolved_frozen)
    except ValueError:
        return
    raise ValueError(
        "Pipeline output directory must not equal or be inside the frozen dataset directory."
    )


def _write_json_atomic(destination: Path, value: object) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        os.close(descriptor)
        temporary.write_text(
            json.dumps(value, indent=2) + "\n",
            encoding="utf-8",
        )
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise


def _copy_effective_manifest(source: Path, destination: Path) -> None:
    descriptor, temporary_name = tempfile.mkstemp(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
    )
    temporary = Path(temporary_name)
    try:
        os.close(descriptor)
        shutil.copyfile(source, temporary)
        temporary.replace(destination)
    except BaseException:
        temporary.unlink(missing_ok=True)
        raise
