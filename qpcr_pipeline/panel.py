from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

DatasetRole = Literal["DESIGN", "CHALLENGE"]
TargetMode = Literal["broad_detection", "subtype_specific"]
Criticality = Literal["CRITICAL", "IMPORTANT", "BACKGROUND"]

_DATASET_ROLES = frozenset({"DESIGN", "CHALLENGE"})
_TARGET_MODES = frozenset({"broad_detection", "subtype_specific"})
_CRITICALITIES = frozenset({"CRITICAL", "IMPORTANT", "BACKGROUND"})


@dataclass(frozen=True, slots=True)
class DiagnosticContext:
    syndrome: str | None = None
    geography: str | None = None
    sample_type: str | None = None
    vector: str | None = None


@dataclass(frozen=True, slots=True)
class SequenceSelectionProvenance:
    dataset_role: DatasetRole
    method: str
    source: str
    details: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class TargetGroup:
    name: str
    required: bool
    dataset_roles: tuple[DatasetRole, ...]
    reasons: tuple[str, ...]
    proposed_by: tuple[str, ...]
    sequence_selection: tuple[SequenceSelectionProvenance, ...] = ()


@dataclass(frozen=True, slots=True)
class PanelTarget:
    name: str
    taxid: int | None
    mode: TargetMode
    subtype: str | None
    groups: tuple[TargetGroup, ...]


@dataclass(frozen=True, slots=True)
class PanelNonTarget:
    name: str
    taxid: int | None
    criticality: Criticality
    dataset_roles: tuple[DatasetRole, ...]
    reasons: tuple[str, ...]
    proposed_by: tuple[str, ...]
    sequence_selection: tuple[SequenceSelectionProvenance, ...] = ()


@dataclass(frozen=True, slots=True)
class PanelDefinition:
    target: PanelTarget
    non_targets: tuple[PanelNonTarget, ...]
    diagnostic_context: DiagnosticContext = DiagnosticContext()


def _validate_non_blank(value: object, label: str) -> None:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-blank string.")


def _validate_taxid(value: object, label: str) -> None:
    if value is None:
        return
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ValueError(f"{label} must be a positive integer when configured.")


def _validate_string_tuple(value: object, label: str, *, allow_empty: bool) -> None:
    if not isinstance(value, tuple):
        raise ValueError(f"{label} must be a tuple.")
    if not allow_empty and not value:
        raise ValueError(f"{label} must not be empty.")
    if any(not isinstance(item, str) or not item.strip() for item in value):
        raise ValueError(f"{label} must contain only non-blank strings.")


def _validate_dataset_roles(value: object, label: str) -> None:
    _validate_string_tuple(value, label, allow_empty=False)
    assert isinstance(value, tuple)
    if any(item not in _DATASET_ROLES for item in value):
        raise ValueError(f"{label} contains an unsupported dataset role.")
    if len(set(value)) != len(value):
        raise ValueError(f"{label} values must be unique.")


def _validate_sequence_selection(
    value: object,
    *,
    label: str,
) -> None:
    if not isinstance(value, tuple):
        raise ValueError(f"{label} must be a tuple.")
    for index, item in enumerate(value, 1):
        if not isinstance(item, SequenceSelectionProvenance):
            raise ValueError(
                f"{label} entry {index} must be a SequenceSelectionProvenance."
            )
        if item.dataset_role not in _DATASET_ROLES:
            raise ValueError(
                f"{label} entry {index} dataset_role is unsupported."
            )
        _validate_non_blank(item.method, f"{label} entry {index} method")
        _validate_non_blank(item.source, f"{label} entry {index} source")
        _validate_string_tuple(
            item.details,
            f"{label} entry {index} details",
            allow_empty=True,
        )


def _validate_context(context: DiagnosticContext) -> None:
    if not isinstance(context, DiagnosticContext):
        raise ValueError("Panel diagnostic_context must be a DiagnosticContext.")
    for name in ("syndrome", "geography", "sample_type", "vector"):
        value = getattr(context, name)
        if value is not None:
            _validate_non_blank(value, f"Panel diagnostic_context.{name}")


def validate_panel_definition(definition: PanelDefinition) -> None:
    if not isinstance(definition, PanelDefinition):
        raise ValueError("Panel definition must be a PanelDefinition.")
    if not isinstance(definition.target, PanelTarget):
        raise ValueError("Panel target must be a PanelTarget.")
    if not isinstance(definition.non_targets, tuple):
        raise ValueError("Panel non_targets must be a tuple.")

    target = definition.target
    _validate_non_blank(target.name, "Panel target name")
    _validate_taxid(target.taxid, "Panel target TaxID")
    if target.mode not in _TARGET_MODES:
        raise ValueError("Panel target mode is unsupported.")
    if target.mode == "broad_detection" and target.subtype is not None:
        raise ValueError("Panel broad_detection target cannot specify a subtype.")
    if target.mode == "subtype_specific":
        _validate_non_blank(target.subtype, "Panel subtype_specific target subtype")
    if not isinstance(target.groups, tuple) or not target.groups:
        raise ValueError("Panel target groups must be a non-empty tuple.")

    group_names: set[str] = set()
    for index, group in enumerate(target.groups, 1):
        if not isinstance(group, TargetGroup):
            raise ValueError(f"Panel target group {index} must be a TargetGroup.")
        _validate_non_blank(group.name, f"Panel target group {index} name")
        folded = group.name.strip().casefold()
        if folded in group_names:
            raise ValueError("Panel target group names must be unique.")
        group_names.add(folded)
        if type(group.required) is not bool:
            raise ValueError(
                f"Panel target group {index} required must be a boolean."
            )
        _validate_dataset_roles(
            group.dataset_roles,
            f"Panel target group {index} dataset_roles",
        )
        _validate_string_tuple(
            group.reasons,
            f"Panel target group {index} reasons",
            allow_empty=False,
        )
        _validate_string_tuple(
            group.proposed_by,
            f"Panel target group {index} proposed_by",
            allow_empty=False,
        )
        _validate_sequence_selection(
            group.sequence_selection,
            label=f"Panel target group {index} sequence_selection",
        )

    non_target_names: set[str] = set()
    non_target_taxids: set[int] = set()
    target_name = target.name.strip().casefold()
    for index, item in enumerate(definition.non_targets, 1):
        if not isinstance(item, PanelNonTarget):
            raise ValueError(f"Panel non-target {index} must be a PanelNonTarget.")
        _validate_non_blank(item.name, f"Panel non-target {index} name")
        _validate_taxid(item.taxid, f"Panel non-target {index} TaxID")
        folded = item.name.strip().casefold()
        if folded == target_name:
            raise ValueError("Panel target and non-target names must be distinct.")
        if folded in non_target_names:
            raise ValueError("Panel non-target names must be unique.")
        non_target_names.add(folded)
        if item.taxid is not None:
            if target.taxid is not None and item.taxid == target.taxid:
                raise ValueError(
                    "Panel target and non-target TaxIDs must be distinct."
                )
            if item.taxid in non_target_taxids:
                raise ValueError("Panel non-target TaxIDs must be unique.")
            non_target_taxids.add(item.taxid)
        if item.criticality not in _CRITICALITIES:
            raise ValueError(f"Panel non-target {index} criticality is unsupported.")
        _validate_dataset_roles(
            item.dataset_roles,
            f"Panel non-target {index} dataset_roles",
        )
        _validate_string_tuple(
            item.reasons,
            f"Panel non-target {index} reasons",
            allow_empty=False,
        )
        _validate_string_tuple(
            item.proposed_by,
            f"Panel non-target {index} proposed_by",
            allow_empty=False,
        )
        _validate_sequence_selection(
            item.sequence_selection,
            label=f"Panel non-target {index} sequence_selection",
        )

    _validate_context(definition.diagnostic_context)
