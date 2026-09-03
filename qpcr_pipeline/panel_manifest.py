from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import re
import tempfile
from typing import Literal

import yaml

from qpcr_pipeline.panel import (
    DiagnosticContext,
    PanelDefinition,
    PanelNonTarget,
    PanelTarget,
    SequenceSelectionProvenance,
    TargetGroup,
    validate_panel_definition,
)

_PROPOSAL_FIELDS = {"schema_version", "status", "definition"}
_APPROVED_FIELDS = {
    "schema_version",
    "status",
    "approved_by_user",
    "proposal_sha256",
    "definition",
}
_DEFINITION_FIELDS = {"target", "non_targets", "diagnostic_context"}
_TARGET_FIELDS = {"name", "taxid", "mode", "subtype", "groups"}
_TARGET_GROUP_FIELDS = {
    "name",
    "required",
    "dataset_roles",
    "reasons",
    "proposed_by",
    "sequence_selection",
}
_NON_TARGET_FIELDS = {
    "name",
    "taxid",
    "criticality",
    "dataset_roles",
    "reasons",
    "proposed_by",
    "sequence_selection",
}
_CONTEXT_FIELDS = {"syndrome", "geography", "sample_type", "vector"}
_SELECTION_FIELDS = {"dataset_role", "method", "source", "details"}
_SHA256_PATTERN = re.compile(r"^sha256:[0-9a-f]{64}$")


@dataclass(frozen=True, slots=True)
class PanelProposal:
    schema_version: int
    status: Literal["PROPOSED"]
    definition: PanelDefinition


@dataclass(frozen=True, slots=True)
class ApprovedPanelManifest:
    schema_version: int
    status: Literal["APPROVED"]
    approved_by_user: bool
    proposal_sha256: str
    definition: PanelDefinition


@dataclass(frozen=True, slots=True)
class PanelResult:
    status: Literal["APPROVED", "LEGACY"]
    manifest_sha256: str | None
    manifest_path: Path | None
    target_mode: str | None
    non_target_count: int


def _require_fields(raw: object, expected: set[str], label: str) -> dict:
    if not isinstance(raw, dict):
        raise ValueError(f"{label} must be a mapping.")
    if set(raw) != expected:
        rendered = ", ".join(sorted(expected))
        raise ValueError(f"{label} fields must be exactly: {rendered}.")
    return raw


def _require_list(raw: object, label: str) -> list:
    if not isinstance(raw, list):
        raise ValueError(f"{label} must be a list.")
    return raw


def _parse_sequence_selection(raw: object, label: str) -> SequenceSelectionProvenance:
    values = _require_fields(raw, _SELECTION_FIELDS, label)
    details = tuple(_require_list(values["details"], f"{label} details"))
    return SequenceSelectionProvenance(
        dataset_role=values["dataset_role"],
        method=values["method"],
        source=values["source"],
        details=details,
    )


def _parse_target_group(raw: object, label: str) -> TargetGroup:
    values = _require_fields(raw, _TARGET_GROUP_FIELDS, label)
    selections = tuple(
        _parse_sequence_selection(item, f"{label} sequence_selection entry {index}")
        for index, item in enumerate(
            _require_list(
                values["sequence_selection"],
                f"{label} sequence_selection",
            ),
            1,
        )
    )
    return TargetGroup(
        name=values["name"],
        required=values["required"],
        dataset_roles=tuple(
            _require_list(values["dataset_roles"], f"{label} dataset_roles")
        ),
        reasons=tuple(_require_list(values["reasons"], f"{label} reasons")),
        proposed_by=tuple(
            _require_list(values["proposed_by"], f"{label} proposed_by")
        ),
        sequence_selection=selections,
    )


def _parse_target(raw: object) -> PanelTarget:
    label = "panel target"
    values = _require_fields(raw, _TARGET_FIELDS, label)
    groups = tuple(
        _parse_target_group(item, f"{label} group {index}")
        for index, item in enumerate(
            _require_list(values["groups"], f"{label} groups"),
            1,
        )
    )
    return PanelTarget(
        name=values["name"],
        taxid=values["taxid"],
        mode=values["mode"],
        subtype=values["subtype"],
        groups=groups,
    )


def _parse_non_target(raw: object, label: str) -> PanelNonTarget:
    values = _require_fields(raw, _NON_TARGET_FIELDS, label)
    selections = tuple(
        _parse_sequence_selection(item, f"{label} sequence_selection entry {index}")
        for index, item in enumerate(
            _require_list(
                values["sequence_selection"],
                f"{label} sequence_selection",
            ),
            1,
        )
    )
    return PanelNonTarget(
        name=values["name"],
        taxid=values["taxid"],
        criticality=values["criticality"],
        dataset_roles=tuple(
            _require_list(values["dataset_roles"], f"{label} dataset_roles")
        ),
        reasons=tuple(_require_list(values["reasons"], f"{label} reasons")),
        proposed_by=tuple(
            _require_list(values["proposed_by"], f"{label} proposed_by")
        ),
        sequence_selection=selections,
    )


def _parse_context(raw: object) -> DiagnosticContext:
    label = "panel diagnostic_context"
    values = _require_fields(raw, _CONTEXT_FIELDS, label)
    return DiagnosticContext(
        syndrome=values["syndrome"],
        geography=values["geography"],
        sample_type=values["sample_type"],
        vector=values["vector"],
    )


def _parse_definition(raw: object) -> PanelDefinition:
    values = _require_fields(raw, _DEFINITION_FIELDS, "panel definition")
    non_targets = tuple(
        _parse_non_target(item, f"panel non-target {index}")
        for index, item in enumerate(
            _require_list(values["non_targets"], "panel non_targets"),
            1,
        )
    )
    definition = PanelDefinition(
        target=_parse_target(values["target"]),
        non_targets=non_targets,
        diagnostic_context=_parse_context(values["diagnostic_context"]),
    )
    validate_panel_definition(definition)
    return definition


def _parse_proposal(raw: object) -> PanelProposal:
    values = _require_fields(raw, _PROPOSAL_FIELDS, "panel proposal")
    if type(values["schema_version"]) is not int or values["schema_version"] != 1:
        raise ValueError("Panel proposal schema_version must be exactly 1.")
    if values["status"] != "PROPOSED":
        raise ValueError("Panel proposal status must be PROPOSED.")
    return PanelProposal(
        schema_version=1,
        status="PROPOSED",
        definition=_parse_definition(values["definition"]),
    )


def _parse_approved(raw: object) -> ApprovedPanelManifest:
    values = _require_fields(raw, _APPROVED_FIELDS, "approved panel manifest")
    if type(values["schema_version"]) is not int or values["schema_version"] != 1:
        raise ValueError("The approved panel schema_version must be exactly 1.")
    if values["status"] != "APPROVED":
        raise ValueError("The approved panel status must be APPROVED.")
    if values["approved_by_user"] is not True:
        raise ValueError("approved_by_user must be exactly true.")
    proposal_sha256 = values["proposal_sha256"]
    if not isinstance(proposal_sha256, str) or not _SHA256_PATTERN.fullmatch(
        proposal_sha256
    ):
        raise ValueError("Approved panel proposal_sha256 is invalid.")
    return ApprovedPanelManifest(
        schema_version=1,
        status="APPROVED",
        approved_by_user=True,
        proposal_sha256=proposal_sha256,
        definition=_parse_definition(values["definition"]),
    )


def _selection_payload(value: SequenceSelectionProvenance) -> dict[str, object]:
    return {
        "dataset_role": value.dataset_role,
        "method": value.method,
        "source": value.source,
        "details": list(value.details),
    }


def _definition_payload(definition: PanelDefinition) -> dict[str, object]:
    validate_panel_definition(definition)
    target = definition.target
    return {
        "target": {
            "name": target.name,
            "taxid": target.taxid,
            "mode": target.mode,
            "subtype": target.subtype,
            "groups": [
                {
                    "name": group.name,
                    "required": group.required,
                    "dataset_roles": list(group.dataset_roles),
                    "reasons": list(group.reasons),
                    "proposed_by": list(group.proposed_by),
                    "sequence_selection": [
                        _selection_payload(item)
                        for item in group.sequence_selection
                    ],
                }
                for group in target.groups
            ],
        },
        "non_targets": [
            {
                "name": item.name,
                "taxid": item.taxid,
                "criticality": item.criticality,
                "dataset_roles": list(item.dataset_roles),
                "reasons": list(item.reasons),
                "proposed_by": list(item.proposed_by),
                "sequence_selection": [
                    _selection_payload(selection)
                    for selection in item.sequence_selection
                ],
            }
            for item in definition.non_targets
        ],
        "diagnostic_context": {
            "syndrome": definition.diagnostic_context.syndrome,
            "geography": definition.diagnostic_context.geography,
            "sample_type": definition.diagnostic_context.sample_type,
            "vector": definition.diagnostic_context.vector,
        },
    }


def _proposal_payload(proposal: PanelProposal) -> dict[str, object]:
    if not isinstance(proposal, PanelProposal):
        raise ValueError("Panel proposal must be a PanelProposal.")
    if type(proposal.schema_version) is not int or proposal.schema_version != 1:
        raise ValueError("Panel proposal schema_version must be exactly 1.")
    if proposal.status != "PROPOSED":
        raise ValueError("Panel proposal status must be PROPOSED.")
    return {
        "schema_version": proposal.schema_version,
        "status": proposal.status,
        "definition": _definition_payload(proposal.definition),
    }


def _approved_payload(manifest: ApprovedPanelManifest) -> dict[str, object]:
    if not isinstance(manifest, ApprovedPanelManifest):
        raise ValueError("Approved panel manifest must be an ApprovedPanelManifest.")
    if type(manifest.schema_version) is not int or manifest.schema_version != 1:
        raise ValueError("The approved panel schema_version must be exactly 1.")
    if manifest.status != "APPROVED":
        raise ValueError("The approved panel status must be APPROVED.")
    if manifest.approved_by_user is not True:
        raise ValueError("approved_by_user must be exactly true.")
    if not isinstance(manifest.proposal_sha256, str) or not _SHA256_PATTERN.fullmatch(
        manifest.proposal_sha256
    ):
        raise ValueError("Approved panel proposal_sha256 is invalid.")
    return {
        "schema_version": manifest.schema_version,
        "status": manifest.status,
        "approved_by_user": manifest.approved_by_user,
        "proposal_sha256": manifest.proposal_sha256,
        "definition": _definition_payload(manifest.definition),
    }


def _canonical_json_text(payload: object) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ) + "\n"


def _approved_text(manifest: ApprovedPanelManifest) -> str:
    return json.dumps(
        _approved_payload(manifest),
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
    ) + "\n"


def _atomic_write_text(path: Path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.",
        suffix=".tmp",
        dir=path.parent,
        text=True,
    )
    temporary = Path(temporary_name)
    try:
        with open(
            descriptor,
            "w",
            encoding="utf-8",
            newline="\n",
            closefd=True,
        ) as handle:
            handle.write(text)
            handle.flush()
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def load_panel_proposal(path: Path) -> PanelProposal:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return _parse_proposal(raw)


def write_panel_proposal(proposal: PanelProposal, path: Path) -> None:
    payload = _proposal_payload(proposal)
    text = yaml.safe_dump(payload, sort_keys=False, allow_unicode=True)
    _atomic_write_text(path, text.rstrip("\n") + "\n")


def proposal_semantic_sha256(proposal: PanelProposal) -> str:
    text = _canonical_json_text(_proposal_payload(proposal))
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return "sha256:" + digest


def load_approved_panel_manifest(path: Path) -> ApprovedPanelManifest:
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    return _parse_approved(raw)


def write_approved_panel_manifest(
    manifest: ApprovedPanelManifest,
    path: Path,
) -> None:
    _atomic_write_text(path, _approved_text(manifest))


def approve_panel_proposal(
    proposal_path: Path,
    output_path: Path,
) -> ApprovedPanelManifest:
    proposal = load_panel_proposal(proposal_path)
    manifest = ApprovedPanelManifest(
        schema_version=1,
        status="APPROVED",
        approved_by_user=True,
        proposal_sha256=proposal_semantic_sha256(proposal),
        definition=proposal.definition,
    )
    write_approved_panel_manifest(manifest, output_path)
    return manifest


def materialize_approved_panel(source_path: Path, outdir: Path) -> PanelResult:
    manifest = load_approved_panel_manifest(source_path)
    destination = Path(outdir) / "panel" / "approved_panel.json"
    write_approved_panel_manifest(manifest, destination)
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    return PanelResult(
        status="APPROVED",
        manifest_sha256="sha256:" + digest,
        manifest_path=destination,
        target_mode=manifest.definition.target.mode,
        non_target_count=len(manifest.definition.non_targets),
    )
