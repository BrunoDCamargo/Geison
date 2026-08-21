"""Minimal pipeline orchestration used by the first end-to-end tracer bullet."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.models import TargetSequenceSet


@dataclass(frozen=True, slots=True)
class RunSummary:
    status: str
    target_name: str
    sequence_count: int
    sequence_ids: list[str]


def run_pipeline(config: PipelineConfig, outdir: str | Path) -> RunSummary:
    sequence_ids = _read_fasta_ids(config.input_fasta)
    target = TargetSequenceSet(sequence_ids=tuple(sequence_ids))

    summary = RunSummary(
        status="COMPLETED",
        target_name=config.target_name,
        sequence_count=len(target.sequence_ids),
        sequence_ids=list(target.sequence_ids),
    )

    output_dir = Path(outdir)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "run_summary.json").write_text(
        json.dumps(asdict(summary), indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


def _read_fasta_ids(path: str | Path) -> list[str]:
    fasta_path = Path(path)
    sequence_ids: list[str] = []

    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        if not line.startswith(">"):
            continue
        sequence_id = line[1:].strip().split(maxsplit=1)[0]
        if sequence_id:
            sequence_ids.append(sequence_id)

    return sequence_ids
