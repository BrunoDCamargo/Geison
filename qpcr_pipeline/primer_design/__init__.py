"""Primer-design API backed by shared and contrastive region selection."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Literal

from . import legacy as _legacy
from qpcr_pipeline.region_selection import (
    CandidateRegion,
    RegionSelectionError,
    candidate_region_from_window,
    is_target_eligible,
    overlap_fraction,
    select_conservation_candidate_regions,
)

if TYPE_CHECKING:
    from qpcr_pipeline.contrastive_conservation import ContrastiveConservationResult
    from qpcr_pipeline.primer3 import Primer3Runner


for _name in dir(_legacy):
    if not _name.startswith("_") and _name not in {"CandidateRegion", "PrimerDesignResult"}:
        globals()[_name] = getattr(_legacy, _name)

PrimerDesignError = _legacy.PrimerDesignError
DesignedOligo = _legacy.DesignedOligo
AssayCandidate = _legacy.AssayCandidate


@dataclass(frozen=True, slots=True)
class PrimerDesignResult:
    status: Literal["SKIPPED", "COMPLETE"]
    reference_id: str | None
    candidates: tuple[CandidateRegion, ...]
    assays: tuple[AssayCandidate, ...]
    candidate_regions_path: Path | None
    assays_path: Path | None
    primer3_input_path: Path | None
    primer3_output_path: Path | None
    report_path: Path
    candidate_source: Literal["CONSERVATION_ONLY", "CONTRASTIVE_CONSERVATION"] = (
        "CONSERVATION_ONLY"
    )


def _select_candidate_regions(conservation, config):
    try:
        return select_conservation_candidate_regions(conservation, config)
    except RegionSelectionError as error:
        raise PrimerDesignError(str(error)) from error


def _candidates_and_source(
    conservation,
    config,
    contrastive: ContrastiveConservationResult | None,
) -> tuple[
    tuple[CandidateRegion, ...],
    Literal["CONSERVATION_ONLY", "CONTRASTIVE_CONSERVATION"],
]:
    if contrastive is not None and contrastive.status == "COMPLETE":
        return (
            tuple(item.region for item in contrastive.candidates),
            "CONTRASTIVE_CONSERVATION",
        )
    return _select_candidate_regions(conservation, config), "CONSERVATION_ONLY"


def primer3_required(
    conservation,
    config,
    *,
    contrastive: ContrastiveConservationResult | None = None,
) -> bool:
    """Return whether this exact primer-design request will invoke Primer3."""
    _legacy.validate_primer_design_config(config)
    if not config.enabled:
        return False
    candidates, _ = _candidates_and_source(conservation, config, contrastive)
    return bool(candidates)


def _report_with_source(
    *,
    status,
    config,
    reference_id,
    candidates,
    assays,
    primer3_details,
    artifacts,
    candidate_source,
):
    report = _legacy._report(
        status=status,
        config=config,
        reference_id=reference_id,
        candidates=candidates,
        assays=assays,
        primer3_details=primer3_details,
        artifacts=artifacts,
    )
    if status == "COMPLETE":
        report["candidate_source"] = candidate_source
    return report


def design_primers(
    conservation,
    config,
    output_dir: Path,
    *,
    contrastive: ContrastiveConservationResult | None = None,
    runner: Primer3Runner | None = None,
) -> PrimerDesignResult:
    """Design assays from contrastive candidates when available, else legacy regions."""
    _legacy.validate_primer_design_config(config)
    output_dir = Path(output_dir)
    paths = _legacy._artifact_paths(output_dir)

    if not config.enabled:
        candidate_source: Literal["CONSERVATION_ONLY", "CONTRASTIVE_CONSERVATION"] = (
            "CONSERVATION_ONLY"
        )
        report = _report_with_source(
            status="SKIPPED",
            config=config,
            reference_id=None,
            candidates=(),
            assays=(),
            primer3_details={},
            artifacts={
                "candidate_regions": None,
                "assays": None,
                "primer3_input": None,
                "primer3_output": None,
            },
            candidate_source=candidate_source,
        )
        paths["report"].parent.mkdir(parents=True, exist_ok=True)
        paths["report"].unlink(missing_ok=True)
        for key in ("candidates", "assays", "input", "output"):
            paths[key].unlink(missing_ok=True)
        _legacy._atomic_write_text(
            paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        return PrimerDesignResult(
            status="SKIPPED",
            reference_id=None,
            candidates=(),
            assays=(),
            candidate_regions_path=None,
            assays_path=None,
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=paths["report"],
            candidate_source=candidate_source,
        )

    candidates, candidate_source = _candidates_and_source(
        conservation, config, contrastive
    )
    if not candidates:
        report = _report_with_source(
            status="COMPLETE",
            config=config,
            reference_id=conservation.reference_id,
            candidates=(),
            assays=(),
            primer3_details={},
            artifacts={
                "candidate_regions": _legacy._relative_path(
                    paths["candidates"], output_dir
                ),
                "assays": _legacy._relative_path(paths["assays"], output_dir),
                "primer3_input": None,
                "primer3_output": None,
            },
            candidate_source=candidate_source,
        )
        paths["report"].parent.mkdir(parents=True, exist_ok=True)
        paths["report"].unlink(missing_ok=True)
        _legacy._atomic_write_text(paths["candidates"], _legacy._candidate_text(()))
        _legacy._atomic_write_text(paths["assays"], _legacy._assay_text(()))
        paths["input"].unlink(missing_ok=True)
        paths["output"].unlink(missing_ok=True)
        _legacy._atomic_write_text(
            paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
        )
        return PrimerDesignResult(
            status="COMPLETE",
            reference_id=conservation.reference_id,
            candidates=(),
            assays=(),
            candidate_regions_path=paths["candidates"],
            assays_path=paths["assays"],
            primer3_input_path=None,
            primer3_output_path=None,
            report_path=paths["report"],
            candidate_source=candidate_source,
        )

    from qpcr_pipeline.primer3 import (
        SubprocessPrimer3Runner,
        build_primer3_input,
        parse_primer3_output,
    )

    paths["report"].parent.mkdir(parents=True, exist_ok=True)
    paths["report"].unlink(missing_ok=True)
    input_text = build_primer3_input(
        conservation.major_consensus, candidates, config
    )
    if runner is None:
        runner = SubprocessPrimer3Runner()
    output_text = runner.run(input_text)
    assays, primer3_details = parse_primer3_output(
        output_text, candidates, conservation.major_consensus
    )
    report = _report_with_source(
        status="COMPLETE",
        config=config,
        reference_id=conservation.reference_id,
        candidates=candidates,
        assays=assays,
        primer3_details=primer3_details,
        artifacts={
            "candidate_regions": _legacy._relative_path(paths["candidates"], output_dir),
            "assays": _legacy._relative_path(paths["assays"], output_dir),
            "primer3_input": _legacy._relative_path(paths["input"], output_dir),
            "primer3_output": _legacy._relative_path(paths["output"], output_dir),
        },
        candidate_source=candidate_source,
    )
    _legacy._atomic_write_text(paths["candidates"], _legacy._candidate_text(candidates))
    _legacy._atomic_write_text(paths["assays"], _legacy._assay_text(assays))
    _legacy._atomic_write_text(paths["input"], input_text)
    _legacy._atomic_write_text(paths["output"], output_text)
    _legacy._atomic_write_text(
        paths["report"], json.dumps(report, indent=2, sort_keys=True) + "\n"
    )
    return PrimerDesignResult(
        status="COMPLETE",
        reference_id=conservation.reference_id,
        candidates=candidates,
        assays=assays,
        candidate_regions_path=paths["candidates"],
        assays_path=paths["assays"],
        primer3_input_path=paths["input"],
        primer3_output_path=paths["output"],
        report_path=paths["report"],
        candidate_source=candidate_source,
    )


# Preserve compatibility for code inside the legacy module that resolves these
# names dynamically while this package owns the public API.
_legacy.CandidateRegion = CandidateRegion
_legacy._select_candidate_regions = _select_candidate_regions

__all__ = [name for name in globals() if not name.startswith("_")]
