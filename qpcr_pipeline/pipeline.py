"""Minimal pipeline orchestration used by the first end-to-end tracer bullet."""

from __future__ import annotations

import json
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.local_input import load_genbank, load_local_sequences
from qpcr_pipeline.ncbi import NcbiClient, acquire_ncbi_dataset, validate_frozen_dataset
from qpcr_pipeline.qc import evaluate_sequences


@dataclass(frozen=True, slots=True)
class RunSummary:
    status: str
    target_name: str
    sequence_count: int
    sequence_ids: list[str]


def run_pipeline(
    config: PipelineConfig, outdir: str | Path, *, ncbi_client: NcbiClient | None = None
) -> RunSummary:
    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_input = config.selected_input
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
        shutil.copyfile(acquired.manifest_path, output_dir / "ncbi_dataset_manifest.json")
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
    }

    (output_dir / "run_summary.json").write_text(
        json.dumps(asdict(summary), indent=2) + "\n",
        encoding="utf-8",
    )
    (output_dir / "qc_report.json").write_text(
        json.dumps(qc_report, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary
