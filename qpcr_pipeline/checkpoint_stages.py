"""Stage-specific checkpoint metadata, scientific input identities, and outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import importlib.metadata
import re
import shutil
import subprocess
from pathlib import Path
from typing import Mapping, Protocol

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
from qpcr_pipeline.checkpointing import (
    CheckpointManifest,
    CheckpointRequest,
    file_sha256,
)
from qpcr_pipeline.config import PipelineConfig
from qpcr_pipeline.execution import STAGE_DEPENDENCIES, STAGE_ORDER
from qpcr_pipeline.ncbi import validate_frozen_dataset
from qpcr_pipeline.primer_design import primer3_required


class ToolIdentityProvider(Protocol):
    def identity(self, tool_name: str) -> Mapping[str, str]: ...


@dataclass(frozen=True, slots=True)
class StageCheckpointDefinition:
    name: str
    dependencies: tuple[str, ...]
    codec: object


_CODECS = {
    "input": INPUT_CODEC,
    "qc": QC_CODEC,
    "clustering": CLUSTERING_CODEC,
    "alignment": ALIGNMENT_CODEC,
    "conservation": CONSERVATION_CODEC,
    "primer_design": PRIMER_DESIGN_CODEC,
    "inclusivity": INCLUSIVITY_CODEC,
    "specificity": SPECIFICITY_CODEC,
    "ranking": RANKING_CODEC,
}

STAGE_DEFINITIONS: dict[str, StageCheckpointDefinition] = {
    stage: StageCheckpointDefinition(
        name=stage,
        dependencies=tuple(STAGE_DEPENDENCIES[stage]),
        codec=_CODECS[stage],
    )
    for stage in STAGE_ORDER
}


def geison_version() -> str:
    return importlib.metadata.version("geison-qpcr")


def _off_target_parameters(config: PipelineConfig) -> list[dict[str, str]]:
    values: list[dict[str, str]] = []
    for item in config.off_targets:
        source = "FASTA" if item.fasta is not None else "NCBI_FROZEN"
        values.append({"name": item.name, "source": source})
    return values


def _input_parameters(config: PipelineConfig) -> dict[str, object]:
    selected = config.selected_input
    if isinstance(selected, tuple):
        _, file_format = selected
        return {"mode": "LOCAL", "format": file_format.upper()}
    if selected.frozen_dataset is not None:
        return {"mode": "NCBI_FROZEN"}
    base: dict[str, object] = {
        "batch_size": selected.batch_size,
        "retries": selected.retries,
        "max_records": selected.max_records,
    }
    if selected.query is not None:
        return {"mode": "NCBI_QUERY", "query": selected.query, **base}
    return {
        "mode": "NCBI_ACCESSIONS",
        "accessions": list(selected.accessions),
        **base,
    }


def stage_parameters(stage: str, config: PipelineConfig) -> Mapping[str, object]:
    if stage not in STAGE_DEFINITIONS:
        raise ValueError(f"unknown checkpoint stage: {stage}")
    if stage == "input":
        return _input_parameters(config)
    if stage == "qc":
        return asdict(config.qc)
    if stage == "clustering":
        return asdict(config.clustering)
    if stage == "alignment":
        return asdict(config.alignment)
    if stage == "conservation":
        return {"target_name": config.target_name, "config": asdict(config.conservation)}
    if stage == "primer_design":
        return asdict(config.primer_design)
    if stage == "inclusivity":
        return asdict(config.inclusivity)
    if stage == "specificity":
        parameters: dict[str, object] = {"config": asdict(config.specificity)}
        if config.specificity.enabled:
            parameters["off_targets"] = _off_target_parameters(config)
        return parameters
    if stage == "ranking":
        return {"target_name": config.target_name, "config": asdict(config.ranking)}
    raise AssertionError(stage)


def stage_input_identities(stage: str, config: PipelineConfig) -> Mapping[str, object]:
    if stage not in STAGE_DEFINITIONS:
        raise ValueError(f"unknown checkpoint stage: {stage}")
    if stage == "input":
        selected = config.selected_input
        if isinstance(selected, tuple):
            path, file_format = selected
            return {
                "target_source": {
                    "mode": "LOCAL",
                    "format": file_format.upper(),
                    "sha256": file_sha256(path),
                }
            }
        if selected.frozen_dataset is not None:
            acquired = validate_frozen_dataset(selected.frozen_dataset)
            return {
                "target_source": {
                    "mode": "NCBI_FROZEN",
                    "records_sha256": file_sha256(acquired.records_path),
                    "manifest_sha256": file_sha256(acquired.manifest_path),
                }
            }
        return {}

    if stage == "specificity":
        if not config.specificity.enabled:
            return {}
        identities: list[dict[str, object]] = []
        for item in config.off_targets:
            if item.fasta is not None:
                identities.append(
                    {
                        "name": item.name,
                        "source": "FASTA",
                        "records_sha256": file_sha256(item.fasta),
                    }
                )
            else:
                assert item.frozen_dataset is not None
                acquired = validate_frozen_dataset(item.frozen_dataset)
                identities.append(
                    {
                        "name": item.name,
                        "source": "NCBI_FROZEN",
                        "records_sha256": file_sha256(acquired.records_path),
                        "manifest_sha256": file_sha256(acquired.manifest_path),
                    }
                )
        return {"off_targets": identities}
    return {}


def stage_tool_identities(
    stage: str,
    config: PipelineConfig,
    stage_context: Mapping[str, object],
    provider: ToolIdentityProvider,
) -> Mapping[str, object]:
    if stage == "clustering":
        qc = stage_context.get("qc")
        if config.clustering.enabled and qc is not None:
            evaluation_set = getattr(qc, "evaluation_set", None)
            if evaluation_set is not None and getattr(evaluation_set, "sequence_ids", ()):
                return {"cd-hit-est": dict(provider.identity("cd-hit-est"))}
        return {}
    if stage == "alignment":
        clustering = stage_context.get("clustering")
        if config.alignment.enabled and clustering is not None:
            discovery_set = getattr(clustering, "discovery_set", None)
            if discovery_set is not None and len(getattr(discovery_set, "sequence_ids", ())) > 1:
                return {"mafft": dict(provider.identity("mafft"))}
        return {}
    if stage == "primer_design":
        conservation = stage_context.get("conservation")
        if conservation is not None and primer3_required(conservation, config.primer_design):
            return {"primer3_core": dict(provider.identity("primer3_core"))}
        return {}
    if stage not in STAGE_DEFINITIONS:
        raise ValueError(f"unknown checkpoint stage: {stage}")
    return {}


class SubprocessToolIdentityProvider:
    _COMMANDS = {
        "cd-hit-est": ("-h",),
        "mafft": ("--version",),
        "primer3_core": ("--version",),
    }

    def identity(self, tool_name: str) -> Mapping[str, str]:
        if tool_name not in self._COMMANDS:
            raise ValueError(f"unsupported external tool identity: {tool_name}")
        executable = shutil.which(tool_name)
        if executable is None:
            raise RuntimeError(f"Required external tool {tool_name!r} was not found on PATH.")
        completed = subprocess.run(
            [executable, *self._COMMANDS[tool_name]],
            capture_output=True,
            text=True,
            check=False,
        )
        text = "\n".join(
            part for part in (completed.stdout, completed.stderr) if part
        )[:10_000]
        version = self._normalized_version(tool_name, text)
        if not version:
            raise RuntimeError(f"Could not determine version for external tool {tool_name!r}.")
        return {"name": tool_name, "version": version}

    @staticmethod
    def _normalized_version(tool_name: str, text: str) -> str:
        lines = [" ".join(line.split()) for line in text.splitlines() if line.strip()]
        if tool_name == "cd-hit-est":
            for line in lines:
                match = re.search(r"CD-HIT\s+version\s+([^\s,;]+)", line, re.IGNORECASE)
                if match:
                    return match.group(1)
        if lines:
            return lines[0][:500]
        return ""


def stage_request(
    stage: str,
    config: PipelineConfig,
    dependency_manifests: Mapping[str, CheckpointManifest],
    stage_context: Mapping[str, object],
    tool_provider: ToolIdentityProvider,
) -> CheckpointRequest:
    definition = STAGE_DEFINITIONS.get(stage)
    if definition is None:
        raise ValueError(f"unknown checkpoint stage: {stage}")
    dependencies: dict[str, str] = {}
    for dependency in definition.dependencies:
        manifest = dependency_manifests.get(dependency)
        if manifest is None or manifest.result_fingerprint is None:
            raise ValueError(
                f"checkpoint stage {stage!r} requires completed dependency {dependency!r}"
            )
        dependencies[dependency] = manifest.result_fingerprint
    return CheckpointRequest(
        stage=stage,
        dependencies=dependencies,
        inputs=stage_input_identities(stage, config),
        parameters=stage_parameters(stage, config),
        software={"geison": geison_version()},
        tools=stage_tool_identities(stage, config, stage_context, tool_provider),
    )


def _inside_outdir(path: Path, outdir: Path) -> Path:
    resolved_outdir = Path(outdir).resolve()
    resolved = Path(path).resolve()
    try:
        resolved.relative_to(resolved_outdir)
    except ValueError as error:
        raise ValueError("checkpoint stage output must be inside the output directory") from error
    return Path(path)


def _paths(*values: object) -> tuple[Path, ...]:
    return tuple(value for value in values if isinstance(value, Path))


def stage_outputs(stage: str, result: object, outdir: Path) -> tuple[Path, ...]:
    outdir = Path(outdir)
    if stage == "input":
        candidates = (
            outdir / "ncbi_dataset" / "records.gb",
            outdir / "ncbi_dataset" / "dataset_manifest.json",
            outdir / "ncbi_dataset_manifest.json",
        )
        values = tuple(path for path in candidates if path.is_file())
    elif stage == "qc":
        values = ()
    elif stage == "clustering":
        values = _paths(
            getattr(result, "discovery_fasta_path", None),
            getattr(result, "report_path", None),
            getattr(result, "raw_cluster_path", None),
        )
    elif stage == "alignment":
        values = _paths(
            getattr(result, "alignment_fasta_path", None),
            getattr(result, "coordinate_map_path", None),
            getattr(result, "report_path", None),
        )
    elif stage == "conservation":
        values = _paths(
            getattr(result, "position_metrics_path", None),
            getattr(result, "window_metrics_path", None),
            getattr(result, "major_consensus_path", None),
            getattr(result, "iupac_consensus_path", None),
            getattr(result, "report_path", None),
        )
    elif stage == "primer_design":
        values = _paths(
            getattr(result, "candidate_regions_path", None),
            getattr(result, "assays_path", None),
            getattr(result, "primer3_input_path", None),
            getattr(result, "primer3_output_path", None),
            getattr(result, "report_path", None),
        )
    elif stage == "inclusivity":
        values = _paths(
            getattr(result, "oligo_matches_path", None),
            getattr(result, "assay_inclusivity_path", None),
            getattr(result, "oligo_variations_path", None),
            getattr(result, "degeneracy_proposals_path", None),
            getattr(result, "report_path", None),
        )
    elif stage == "specificity":
        values = _paths(
            getattr(result, "off_target_hits_path", None),
            getattr(result, "plausible_amplicons_path", None),
            getattr(result, "report_path", None),
        )
    elif stage == "ranking":
        values = _paths(
            getattr(result, "ranking_tsv_path", None),
            getattr(result, "ranking_report_path", None),
            getattr(result, "html_report_path", None),
        )
    else:
        raise ValueError(f"unknown checkpoint stage: {stage}")

    seen: set[Path] = set()
    checked: list[Path] = []
    for path in values:
        _inside_outdir(path, outdir)
        resolved = path.resolve()
        if resolved not in seen:
            checked.append(path)
            seen.add(resolved)
    return tuple(checked)
