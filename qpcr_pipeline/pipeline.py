"""Checkpoint-aware orchestration for the Geison qPCR pipeline."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Mapping

from qpcr_pipeline.alignment import MafftRunner, align_discovery
from qpcr_pipeline.checkpoint_stages import (
    STAGE_DEFINITIONS,
    SubprocessToolIdentityProvider,
    ToolIdentityProvider,
    stage_outputs,
    stage_request,
)
from qpcr_pipeline.checkpointing import CheckpointManager, CheckpointManifest
from qpcr_pipeline.clustering import CdHitRunner, cluster_sequences
from qpcr_pipeline.config import NcbiInputConfig, PipelineConfig
from qpcr_pipeline.conservation import analyze_conservation
from qpcr_pipeline.diagnostics import EnvironmentInspector
from qpcr_pipeline.execution import (
    STAGE_ORDER,
    ExecutionPolicy,
    required_reuse_boundary,
    transitive_descendants,
)
from qpcr_pipeline.inclusivity import evaluate_inclusivity
from qpcr_pipeline.local_input import LocalSequenceRecord, load_genbank, load_local_sequences
from qpcr_pipeline.ncbi import NcbiClient, acquire_ncbi_dataset, validate_frozen_dataset
from qpcr_pipeline.primer3 import Primer3Runner
from qpcr_pipeline.primer_design import design_primers
from qpcr_pipeline.qc import evaluate_sequences
from qpcr_pipeline.ranking import evaluate_ranking
from qpcr_pipeline.run_recording import (
    RunRecorder,
    assess_final_completeness,
    assess_pre_ranking_completeness,
)
from qpcr_pipeline.specificity import evaluate_specificity


@dataclass(frozen=True, slots=True)
class StageActionSummary:
    stage: str
    action: str


@dataclass(frozen=True, slots=True)
class RunSummary:
    status: str
    target_name: str
    sequence_count: int
    sequence_ids: list[str]
    stage_actions: tuple[StageActionSummary, ...] = ()


class _EffectiveToolIdentityProvider:
    """Use stable injected-runner identities before consulting real binaries."""

    def __init__(
        self,
        base: ToolIdentityProvider,
        *,
        cdhit_runner: CdHitRunner | None,
        mafft_runner: MafftRunner | None,
        primer3_runner: Primer3Runner | None,
    ) -> None:
        self.base = base
        self.injected = {
            "cd-hit-est": cdhit_runner,
            "mafft": mafft_runner,
            "primer3_core": primer3_runner,
        }

    def identity(self, tool_name: str) -> Mapping[str, str]:
        runner = self.injected.get(tool_name)
        if runner is None:
            return self.base.identity(tool_name)
        runner_type = type(runner)
        identity = f"injected:{runner_type.__module__}.{runner_type.__qualname__}"
        return {"name": tool_name, "version": identity}


def run_pipeline(
    config: PipelineConfig,
    outdir: str | Path,
    *,
    ncbi_client: NcbiClient | None = None,
    cdhit_runner: CdHitRunner | None = None,
    mafft_runner: MafftRunner | None = None,
    primer3_runner: Primer3Runner | None = None,
    execution: ExecutionPolicy | None = None,
    tool_identity_provider: ToolIdentityProvider | None = None,
) -> RunSummary:
    policy = execution or ExecutionPolicy()
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

    recorder = RunRecorder(output_dir)
    recorder.begin_attempt(
        config.target_name,
        config,
        policy,
        EnvironmentInspector().inspect(config),
        [],
    )

    base_tool_provider = tool_identity_provider or SubprocessToolIdentityProvider()
    effective_tool_provider = _EffectiveToolIdentityProvider(
        base_tool_provider,
        cdhit_runner=cdhit_runner,
        mafft_runner=mafft_runner,
        primer3_runner=primer3_runner,
    )
    manager = CheckpointManager(output_dir)
    results: dict[str, object] = {}
    manifests: dict[str, CheckpointManifest] = {}
    actions: list[StageActionSummary] = []

    if policy.from_step is not None:
        _preflight_from_step(
            policy.from_step,
            config,
            manager,
            results,
            manifests,
            effective_tool_provider,
        )
        forced = {policy.from_step, *transitive_descendants(policy.from_step)}
        boundary = set(required_reuse_boundary(policy.from_step))
        for stage in STAGE_ORDER:
            if stage in boundary:
                actions.append(StageActionSummary(stage, "REUSE"))
                continue
            if stage not in forced:
                raise RuntimeError(
                    f"Execution graph did not classify stage {stage!r} for --from-step."
                )
            request = stage_request(
                stage, config, manifests, results, effective_tool_provider
            )
            result, manifest = _run_and_checkpoint_stage(
                stage,
                request,
                manager,
                config,
                output_dir,
                results,
                ncbi_client=ncbi_client,
                cdhit_runner=cdhit_runner,
                mafft_runner=mafft_runner,
                primer3_runner=primer3_runner,
                refresh_online_input=(stage == "input"),
            )
            results[stage] = result
            manifests[stage] = manifest
            actions.append(StageActionSummary(stage, "FORCED"))
    elif policy.resume:
        forced: set[str] = set()
        if policy.force_step is not None:
            forced.add(policy.force_step)
            forced.update(transitive_descendants(policy.force_step))

        for stage in STAGE_ORDER:
            request = stage_request(
                stage, config, manifests, results, effective_tool_provider
            )
            if stage not in forced:
                validation = manager.validate(
                    request, STAGE_DEFINITIONS[stage].codec
                )
                if validation.valid:
                    assert validation.loaded is not None
                    results[stage] = validation.loaded.state
                    manifests[stage] = validation.loaded.manifest
                    actions.append(StageActionSummary(stage, "REUSE"))
                    continue

                action = "RUN"
                forced.update(transitive_descendants(stage))
            else:
                action = "FORCED"

            explicit_input_refresh = (
                stage == "input" and policy.force_step == "input"
            )
            result, manifest = _run_and_checkpoint_stage(
                stage,
                request,
                manager,
                config,
                output_dir,
                results,
                ncbi_client=ncbi_client,
                cdhit_runner=cdhit_runner,
                mafft_runner=mafft_runner,
                primer3_runner=primer3_runner,
                refresh_online_input=explicit_input_refresh,
            )
            results[stage] = result
            manifests[stage] = manifest
            actions.append(StageActionSummary(stage, action))
    else:
        for stage in STAGE_ORDER:
            request = stage_request(
                stage, config, manifests, results, effective_tool_provider
            )
            result, manifest = _run_and_checkpoint_stage(
                stage,
                request,
                manager,
                config,
                output_dir,
                results,
                ncbi_client=ncbi_client,
                cdhit_runner=cdhit_runner,
                mafft_runner=mafft_runner,
                primer3_runner=primer3_runner,
                refresh_online_input=(stage == "input"),
            )
            results[stage] = result
            manifests[stage] = manifest
            actions.append(StageActionSummary(stage, "RUN"))

    qc_result = results["qc"]
    clustering = results["clustering"]
    alignment = results["alignment"]
    conservation = results["conservation"]
    primer_design = results["primer_design"]
    inclusivity = results["inclusivity"]
    specificity = results["specificity"]
    ranking = results["ranking"]

    pre_ranking_completeness = assess_pre_ranking_completeness(
        evaluation_sequence_count=len(qc_result.evaluation_set.sequence_ids),
        assay_count=len(primer_design.assays),
        inclusivity_status=inclusivity.status,
        specificity_status=specificity.status,
    )
    scientific_completeness = assess_final_completeness(
        pre_ranking_completeness,
        ranking_status=ranking.status,
    )
    final_status = "COMPLETED" if scientific_completeness.complete else "PARTIAL"
    if not scientific_completeness.complete and any(
        assay.classification == "IN SILICO PASS" for assay in ranking.assays
    ):
        raise RuntimeError("Incomplete run evidence cannot produce IN SILICO PASS.")

    approved_ids = qc_result.evaluation_set.sequence_ids
    summary = RunSummary(
        status=final_status,
        target_name=config.target_name,
        sequence_count=len(approved_ids),
        sequence_ids=list(approved_ids),
        stage_actions=tuple(actions),
    )

    top_recommended_assay_id = next(
        (
            assay.assay_id
            for assay in ranking.assays
            if assay.classification == "IN SILICO PASS"
        ),
        None,
    )
    qc_report = {
        "records": [
            {
                "sequence_id": record.sequence_id,
                "status": record.status.value,
                "reason_codes": list(record.reason_codes),
            }
            for record in qc_result.records
        ],
        "target_sequence_set": {
            "sequence_ids": list(qc_result.target_sequence_set.sequence_ids)
        },
        "evaluation_set": {"sequence_ids": list(approved_ids)},
        "discovery_set": {
            "sequence_ids": list(clustering.discovery_set.sequence_ids)
        },
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
                amplicon.detectable_off_target
                for amplicon in specificity.amplicons
            ),
        },
        "ranking": {
            "status": ranking.status,
            "assay_count": len(ranking.assays),
            "in_silico_pass_count": sum(
                assay.classification == "IN SILICO PASS"
                for assay in ranking.assays
            ),
            "review_count": sum(
                assay.classification == "REVIEW" for assay in ranking.assays
            ),
            "high_risk_count": sum(
                assay.classification == "HIGH_RISK"
                for assay in ranking.assays
            ),
            "complete_score_count": sum(
                assay.score_status == "COMPLETE" for assay in ranking.assays
            ),
            "incomplete_score_count": sum(
                assay.score_status == "INCOMPLETE" for assay in ranking.assays
            ),
            "top_recommended_assay_id": top_recommended_assay_id,
        },
    }

    _write_json_atomic(output_dir / "run_summary.json", asdict(summary))
    _write_json_atomic(output_dir / "qc_report.json", qc_report)
    for action in actions:
        checkpoint_path = output_dir / ".checkpoints" / action.stage / "manifest.json"
        if action.action == "REUSE":
            recorder.stage_reused(action.stage, action.action, checkpoint_path)
        else:
            recorder.stage_started(action.stage, action.action)
            recorder.stage_completed(action.stage, action.action, checkpoint_path)
    recorder.complete(
        final_status,
        scientific_completeness,
        {},
        {"id": alignment.reference_id, "mode": alignment.reference_mode},
    )
    return summary


def _preflight_from_step(
    from_step: str,
    config: PipelineConfig,
    manager: CheckpointManager,
    results: dict[str, object],
    manifests: dict[str, CheckpointManifest],
    tool_provider: ToolIdentityProvider,
) -> None:
    boundary = set(required_reuse_boundary(from_step))
    blocking: list[str] = []
    for stage in STAGE_ORDER:
        if stage not in boundary:
            continue
        missing_dependencies = [
            dependency
            for dependency in STAGE_DEFINITIONS[stage].dependencies
            if dependency in boundary and dependency not in manifests
        ]
        if missing_dependencies:
            blocking.append(
                f"{stage}: DEPENDENCY_UNAVAILABLE ({', '.join(missing_dependencies)})"
            )
            continue
        try:
            request = stage_request(
                stage, config, manifests, results, tool_provider
            )
        except Exception as error:
            blocking.append(
                f"{stage}: REQUEST_INVALID ({type(error).__name__}: {error})"
            )
            continue
        validation = manager.validate(request, STAGE_DEFINITIONS[stage].codec)
        if not validation.valid:
            category = (
                validation.invalidity.value
                if validation.invalidity is not None
                else "INVALID_CHECKPOINT"
            )
            detail = validation.detail or "checkpoint validation failed"
            blocking.append(f"{stage}: {category} ({detail})")
            continue
        assert validation.loaded is not None
        results[stage] = validation.loaded.state
        manifests[stage] = validation.loaded.manifest

    if blocking:
        details = "; ".join(blocking)
        raise ValueError(
            "--from-step cannot start because required checkpoint(s) are invalid: "
            f"{details}. Use --resume to repair invalid stages or restart from an earlier stage."
        )


def _run_and_checkpoint_stage(
    stage: str,
    request,
    manager: CheckpointManager,
    config: PipelineConfig,
    output_dir: Path,
    results: Mapping[str, object],
    *,
    ncbi_client: NcbiClient | None,
    cdhit_runner: CdHitRunner | None,
    mafft_runner: MafftRunner | None,
    primer3_runner: Primer3Runner | None,
    refresh_online_input: bool,
) -> tuple[object, CheckpointManifest]:
    manager.begin(request)
    try:
        result = _run_stage(
            stage,
            config,
            output_dir,
            results,
            ncbi_client=ncbi_client,
            cdhit_runner=cdhit_runner,
            mafft_runner=mafft_runner,
            primer3_runner=primer3_runner,
            refresh_online_input=refresh_online_input,
        )
        manifest = manager.complete(
            request,
            result,
            STAGE_DEFINITIONS[stage].codec,
            stage_outputs(stage, result, output_dir),
        )
        return result, manifest
    except Exception as error:
        manager.fail(request, error)
        raise


def _run_stage(
    stage: str,
    config: PipelineConfig,
    output_dir: Path,
    results: Mapping[str, object],
    *,
    ncbi_client: NcbiClient | None,
    cdhit_runner: CdHitRunner | None,
    mafft_runner: MafftRunner | None,
    primer3_runner: Primer3Runner | None,
    refresh_online_input: bool,
) -> object:
    if stage == "input":
        return _run_input(
            config,
            output_dir,
            ncbi_client=ncbi_client,
            refresh_online=refresh_online_input,
        )

    records = results["input"]
    if stage == "qc":
        return _run_qc(records, config)

    qc_result = results["qc"]
    approved_records = _approved_records(records, qc_result)
    if stage == "clustering":
        return cluster_sequences(
            approved_records,
            qc_result.evaluation_set,
            config.clustering,
            output_dir,
            runner=cdhit_runner,
        )

    clustering = results["clustering"]
    discovery_records = _discovery_records(approved_records, clustering)
    if stage == "alignment":
        return align_discovery(
            discovery_records,
            clustering.discovery_set,
            config.alignment,
            output_dir,
            runner=mafft_runner,
        )

    alignment = results["alignment"]
    if stage == "conservation":
        return analyze_conservation(
            discovery_records,
            alignment,
            config.conservation,
            output_dir,
            target_name=config.target_name,
        )

    conservation = results["conservation"]
    if stage == "primer_design":
        return design_primers(
            conservation,
            config.primer_design,
            output_dir,
            runner=primer3_runner,
        )

    primer_design = results["primer_design"]
    if stage == "inclusivity":
        return evaluate_inclusivity(
            approved_records,
            qc_result.evaluation_set,
            primer_design,
            config.inclusivity,
            output_dir,
        )
    if stage == "specificity":
        return evaluate_specificity(
            primer_design,
            config.off_targets,
            config.specificity,
            output_dir,
        )

    if stage == "ranking":
        return evaluate_ranking(
            primer_design,
            results["inclusivity"],
            results["specificity"],
            config.ranking,
            output_dir,
            target_name=config.target_name,
        )
    raise ValueError(f"Unknown pipeline stage: {stage}")


def _run_input(
    config: PipelineConfig,
    output_dir: Path,
    *,
    ncbi_client: NcbiClient | None,
    refresh_online: bool,
) -> tuple[LocalSequenceRecord, ...]:
    selected_input = config.selected_input
    if isinstance(selected_input, NcbiInputConfig):
        if selected_input.frozen_dataset is not None:
            acquired = validate_frozen_dataset(selected_input.frozen_dataset)
        else:
            dataset_dir = output_dir / "ncbi_dataset"
            if refresh_online and dataset_dir.exists():
                shutil.rmtree(dataset_dir)
            acquired = acquire_ncbi_dataset(
                selected_input,
                dataset_dir,
                client=ncbi_client,
            )
        records = load_genbank(acquired.records_path)
        _copy_effective_manifest(
            acquired.manifest_path, output_dir / "ncbi_dataset_manifest.json"
        )
        return records

    input_path, input_format = selected_input
    return load_local_sequences(input_path, input_format)


def _run_qc(
    records: tuple[LocalSequenceRecord, ...], config: PipelineConfig
):
    return evaluate_sequences(
        records,
        min_length=config.qc.min_length,
        max_ambiguous_fraction=config.qc.max_ambiguous_fraction,
        expected_length=config.qc.expected_length,
        length_tolerance_fraction=config.qc.length_tolerance_fraction,
    )


def _approved_records(records, qc_result):
    approved_ids = set(qc_result.evaluation_set.sequence_ids)
    return tuple(
        record for record in records if record.sequence_id in approved_ids
    )


def _discovery_records(approved_records, clustering):
    discovery_ids = set(clustering.discovery_set.sequence_ids)
    return tuple(
        record
        for record in approved_records
        if record.sequence_id in discovery_ids
    )


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
