"""Checkpoint-aware orchestration for the Geison qPCR pipeline."""

from __future__ import annotations

import json
import os
import shutil
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Mapping

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
from qpcr_pipeline.contrastive_conservation import analyze_contrastive_conservation
from qpcr_pipeline.diagnostics import EnvironmentInspector
from qpcr_pipeline.execution import STAGE_ORDER, ExecutionPolicy
from qpcr_pipeline.inclusivity import evaluate_inclusivity
from qpcr_pipeline.local_input import LocalSequenceRecord, load_genbank, load_local_sequences
from qpcr_pipeline.ncbi import NcbiClient, acquire_ncbi_dataset, validate_frozen_dataset
from qpcr_pipeline.panel_manifest import (
    PanelResult,
    load_approved_panel_manifest,
    materialize_approved_panel,
    prepare_panel_preflight,
)
from qpcr_pipeline.planning import plan_pipeline
from qpcr_pipeline.provenance import (
    build_input_provenance,
    build_panel_provenance,
    build_reference_provenance,
    effective_config_payload,
)
from qpcr_pipeline.primer3 import Primer3Runner
from qpcr_pipeline.primer_design import design_primers
from qpcr_pipeline.qc import evaluate_sequences
from qpcr_pipeline.ranking_guard import evaluate_ranking_with_execution_guard
from qpcr_pipeline.researcher_report import generate_researcher_report
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
    action_required_code: str | None = None
    action_required_artifact: str | None = None


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
    environment_inspector: EnvironmentInspector | None = None,
    recorder_factory: Callable[[Path], RunRecorder] | None = None,
) -> RunSummary:
    policy = execution or ExecutionPolicy()
    selected_input = config.selected_input
    output_dir = Path(outdir)
    if (
        isinstance(selected_input, NcbiInputConfig)
        and selected_input.frozen_dataset is not None
    ):
        _reject_output_inside_frozen_dataset(output_dir, selected_input.frozen_dataset)
    output_dir.mkdir(parents=True, exist_ok=True)

    base_tool_provider = tool_identity_provider or SubprocessToolIdentityProvider()
    effective_tool_provider = _EffectiveToolIdentityProvider(
        base_tool_provider,
        cdhit_runner=cdhit_runner,
        mafft_runner=mafft_runner,
        primer3_runner=primer3_runner,
    )
    environment = (environment_inspector or EnvironmentInspector()).inspect(config)
    panel_preflight = prepare_panel_preflight(
        config.panel,
        output_dir,
        target_name=config.target_name,
    )
    recorder = (
        recorder_factory(output_dir)
        if recorder_factory is not None
        else RunRecorder(output_dir)
    )
    if panel_preflight.status == "ACTION_REQUIRED":
        proposal_path = panel_preflight.proposal_path
        assert proposal_path is not None
        synthetic_plan = [
            {
                "stage": "panel",
                "action": "ACTION_REQUIRED",
                "reason": "panel approval required before scientific execution",
            }
        ]
        recorder.begin_attempt(
            config.target_name,
            effective_config_payload(config),
            asdict(policy),
            environment,
            synthetic_plan,
        )
        (output_dir / "run_summary.json").unlink(missing_ok=True)
        (output_dir / "qc_report.json").unlink(missing_ok=True)
        recorder.action_required("PANEL_APPROVAL_REQUIRED", proposal_path)
        summary = RunSummary(
            status="ACTION_REQUIRED",
            target_name=config.target_name,
            sequence_count=0,
            sequence_ids=[],
            action_required_code="PANEL_APPROVAL_REQUIRED",
            action_required_artifact=str(proposal_path),
        )
        _write_json_atomic(output_dir / "run_summary.json", asdict(summary))
        _publish_researcher_report(output_dir)
        return summary
    plan = plan_pipeline(
        config,
        output_dir,
        execution=policy,
        tool_identity_provider=effective_tool_provider,
    )
    recorder.begin_attempt(
        config.target_name,
        effective_config_payload(config),
        asdict(policy),
        environment,
        [asdict(decision) for decision in plan.decisions],
    )

    (output_dir / "run_summary.json").unlink(missing_ok=True)
    (output_dir / "qc_report.json").unlink(missing_ok=True)

    manager = CheckpointManager(output_dir)
    results: dict[str, object] = dict(plan.reused_results)
    manifests: dict[str, CheckpointManifest] = dict(plan.reused_manifests)
    actions: list[StageActionSummary] = []
    current_stage: str | None = None

    try:
        for decision in plan.decisions:
            stage = decision.stage
            action = decision.action
            checkpoint_path = output_dir / ".checkpoints" / stage / "manifest.json"

            if action == "REUSE":
                if stage not in results or stage not in manifests:
                    raise RuntimeError(
                        f"Execution plan selected REUSE without loaded checkpoint for {stage!r}."
                    )
                recorder.stage_reused(stage, action, checkpoint_path)
                actions.append(StageActionSummary(stage, action))
                continue

            current_stage = stage
            recorder.stage_started(stage, action)
            request = stage_request(
                stage, config, manifests, results, effective_tool_provider
            )
            execution_missing_evidence: tuple[str, ...] = ()
            if stage == "ranking":
                execution_missing_evidence = _pre_ranking_completeness(results).missing_evidence

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
                execution_missing_evidence=execution_missing_evidence,
            )
            results[stage] = result
            manifests[stage] = manifest
            recorder.stage_completed(stage, action, checkpoint_path)
            actions.append(StageActionSummary(stage, action))
            current_stage = None

        qc_result = results["qc"]
        clustering = results["clustering"]
        alignment = results["alignment"]
        conservation = results["conservation"]
        contrastive = results["contrastive_conservation"]
        primer_design = results["primer_design"]
        inclusivity = results["inclusivity"]
        specificity = results["specificity"]
        ranking = results["ranking"]

        pre_ranking = assess_pre_ranking_completeness(
            evaluation_sequence_count=len(qc_result.evaluation_set.sequence_ids),
            assay_count=len(primer_design.assays),
            inclusivity_status=inclusivity.status,
            specificity_status=specificity.status,
        )
        final_completeness = assess_final_completeness(
            pre_ranking,
            ranking_status=ranking.status,
        )
        final_status = "COMPLETED" if final_completeness.complete else "PARTIAL"

        if not final_completeness.complete and any(
            assay.classification == "IN SILICO PASS" for assay in ranking.assays
        ):
            raise RuntimeError(
                "Incomplete run evidence cannot produce IN SILICO PASS."
            )

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
            "contrastive_conservation": {
                "status": contrastive.status,
                "window_count": len(contrastive.windows),
                "challenge_dataset_count": len(contrastive.challenge_datasets),
                "candidate_region_count": len(contrastive.candidates),
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
        recorder.complete(
            final_status,
            final_completeness,
            input_provenance=build_input_provenance(
                config, output_dir, qc_result, manifests["input"]
            ),
            reference=build_reference_provenance(alignment),
            panel_provenance=build_panel_provenance(results["panel"]),
        )
        _publish_researcher_report(output_dir)
        return summary
    except BaseException as error:
        try:
            recorder.fail(error, stage=current_stage)
        except BaseException as diagnostic_error:
            try:
                error.add_note(
                    "Geison could not persist run failure diagnostics: "
                    f"{type(diagnostic_error).__name__}: {diagnostic_error}"
                )
            except Exception:
                pass
        _publish_researcher_report(output_dir)
        raise


def _pre_ranking_completeness(results: Mapping[str, object]):
    qc_result = results["qc"]
    primer_design = results["primer_design"]
    inclusivity = results["inclusivity"]
    specificity = results["specificity"]
    return assess_pre_ranking_completeness(
        evaluation_sequence_count=len(qc_result.evaluation_set.sequence_ids),
        assay_count=len(primer_design.assays),
        inclusivity_status=inclusivity.status,
        specificity_status=specificity.status,
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
    execution_missing_evidence: tuple[str, ...] = (),
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
            execution_missing_evidence=execution_missing_evidence,
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
    execution_missing_evidence: tuple[str, ...] = (),
) -> object:
    if stage == "panel":
        if config.panel is None:
            return PanelResult(
                status="LEGACY",
                manifest_sha256=None,
                manifest_path=None,
                target_mode=None,
                non_target_count=0,
            )
        if config.panel.frozen_manifest is None:
            raise RuntimeError(
                "Panel proposal reached checkpoint execution before approval."
            )
        return materialize_approved_panel(config.panel.frozen_manifest, output_dir)

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
    if stage == "contrastive_conservation":
        approved_panel = None
        if config.contrastive_conservation.enabled:
            if config.panel is None or config.panel.frozen_manifest is None:
                raise ValueError(
                    "Enabled contrastive conservation requires an approved frozen panel."
                )
            approved_panel = load_approved_panel_manifest(config.panel.frozen_manifest)
        return analyze_contrastive_conservation(
            conservation,
            approved_panel,
            config.off_targets,
            config.contrastive_conservation,
            config.primer_design,
            output_dir,
        )

    contrastive = results["contrastive_conservation"]
    if stage == "primer_design":
        return design_primers(
            conservation,
            config.primer_design,
            output_dir,
            contrastive=contrastive,
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
        return evaluate_ranking_with_execution_guard(
            primer_design,
            results["inclusivity"],
            results["specificity"],
            config.ranking,
            output_dir,
            target_name=config.target_name,
            execution_missing_evidence=execution_missing_evidence,
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
    return tuple(record for record in records if record.sequence_id in approved_ids)


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


def _publish_researcher_report(output_dir: Path) -> None:
    manifest_path = output_dir / "run_manifest.json"
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return
    if not isinstance(manifest, dict):
        return
    effective_config = manifest.get("effective_config")
    if not isinstance(effective_config, dict):
        return
    primer_design = effective_config.get("primer_design")
    if not isinstance(primer_design, dict) or primer_design.get("enabled") is not True:
        return

    report_path = output_dir / "report.html"
    error_path = output_dir / "report_error.json"
    report_path.unlink(missing_ok=True)
    error_path.unlink(missing_ok=True)
    try:
        generate_researcher_report(output_dir)
    except Exception as error:
        report_path.unlink(missing_ok=True)
        try:
            _write_json_atomic(
                error_path,
                {
                    "error_type": type(error).__name__,
                    "message": str(error)[:1000],
                },
            )
        except Exception:
            pass


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
