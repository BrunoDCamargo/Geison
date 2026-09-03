# Panel Model and Manifests Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Geison's first panel-aware workflow: explicit target/non-target panel metadata, deterministic proposal and approval artifacts, a first-run human approval gate, a frozen approved manifest, and checkpoint/provenance integration.

**Architecture:** Add a pure panel domain module and a strict manifest I/O module. An inline panel proposal is a review request: `qpcr-pipeline run` writes `panel_proposal.yaml`, records `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED`, and stops before input acquisition. `qpcr-pipeline panel approve` converts a reviewed proposal into deterministic canonical JSON. A frozen approved manifest becomes a checkpointed `panel` stage before `input`; changing that manifest conservatively invalidates every downstream stage. Existing non-design workflows may omit a panel during this migration, while enabled primer design requires either a proposal or an approved manifest.

**Tech Stack:** Python 3.10+, stdlib `dataclasses`, `hashlib`, `json`, `pathlib`, `tempfile`, existing PyYAML dependency, pytest, existing Geison checkpoint/resume and run-recording infrastructure.

**Spec:** `docs/superpowers/specs/2026-09-03-decipher-hybrid-assay-architecture-design.md`

## Global Constraints

- Geison remains the Python orchestrator. This subproject adds no R or DECIPHER dependency.
- The panel model is pathogen-agnostic. The West Nile fixture is test data, not a production clinical knowledge base.
- A configured proposal must never silently enter scientific design.
- Approval is explicit and non-interactive. `qpcr-pipeline run` must never prompt for approval.
- Approved manifests are strict, deterministic, hashable scientific artifacts.
- Design/challenge roles are modeled now. Actual accession assignment and representative sequence partitioning remain subproject 4 work.
- Existing `off_targets` configuration remains unchanged in this subproject.
- `primer_design.enabled: true` requires a configured panel. QC/alignment/conservation-only legacy workflows may remain panel-free during migration.
- No universal biological threshold, automatic clinical inference, or live network lookup is added here.
- `run_manifest.json` stays at schema version 1 in this subproject. Older schema-1 manifests must be normalized with missing new optional fields during resume rather than rejected.
- The new `panel` stage makes old pre-panel checkpoints non-reusable across the new dependency graph. That one-time migration cost is acceptable.

---

## File Map

### New production files

- `qpcr_pipeline/panel.py`: immutable panel models and structural/scientific consistency validation; no file I/O.
- `qpcr_pipeline/panel_manifest.py`: strict proposal YAML, canonical approved JSON, semantic hashing, approval, runtime preflight, and approved-manifest materialization.

### Modified production files

- `qpcr_pipeline/config.py`: `PanelConfig`, strict YAML parsing, target-name consistency, primer-design panel requirement.
- `qpcr_pipeline/execution.py`: add `panel` as the first stage; `input` depends on it.
- `qpcr_pipeline/checkpoint_codecs.py`: add `PANEL_CODEC`.
- `qpcr_pipeline/checkpoint_stages.py`: add panel parameters, source identity, codec, and output identity.
- `qpcr_pipeline/pipeline.py`: panel preflight, action-required summary, panel stage execution, provenance call.
- `qpcr_pipeline/run_recording.py`: action-required lifecycle and panel provenance field.
- `qpcr_pipeline/provenance.py`: safe panel provenance projection.
- `qpcr_pipeline/dry_run.py`: non-mutating panel preflight reporting.
- `qpcr_pipeline/cli.py`: nested `panel approve` command and action-required exit code.
- `README.md`: proposal → review → approve → rerun workflow.

### New tests and fixture

- `tests/test_panel.py`
- `tests/test_panel_manifest.py`
- `tests/test_panel_config.py`
- `tests/test_pipeline_panel.py`
- `tests/test_panel_cli.py`
- `tests/fixtures/panels/west_nile_proposal.yaml`

### Existing tests that must be updated where their literal assumptions change

- `tests/test_execution.py`
- `tests/test_execution_plan.py`
- `tests/test_checkpoint_codecs.py`
- `tests/test_checkpoint_stages.py`
- `tests/test_pipeline_resume.py`
- `tests/test_dry_run.py`
- `tests/test_run_manifest.py`
- `tests/test_run_recording.py`
- `tests/test_cli.py`
- `tests/test_minimal_run.py`
- `tests/test_pipeline_ranking.py`
- `tests/test_pipeline_specificity.py`
- Any test discovered by `grep -R "PrimerDesignConfig(enabled=True\|primer_design:" -n tests` that constructs a full `PipelineConfig` with enabled primer design.

---

### Task 1: Add immutable panel models and validation

**Files:**
- Create: `qpcr_pipeline/panel.py`
- Create: `tests/test_panel.py`
- Create: `tests/fixtures/panels/west_nile_proposal.yaml`

**Interfaces:**
- Produces `DatasetRole`, `TargetMode`, `Criticality`, `DiagnosticContext`, `SequenceSelectionProvenance`, `TargetGroup`, `PanelTarget`, `PanelNonTarget`, `PanelDefinition`, and `validate_panel_definition()`.
- Consumes no new project interface.

- [ ] **Step 1: Write failing domain-model tests**

Create `tests/test_panel.py` with a small synthetic panel so the validation tests do not depend on disputed biological metadata:

```python
import pytest

from qpcr_pipeline.panel import (
    DiagnosticContext,
    PanelDefinition,
    PanelNonTarget,
    PanelTarget,
    SequenceSelectionProvenance,
    TargetGroup,
    validate_panel_definition,
)


def valid_panel() -> PanelDefinition:
    selection = SequenceSelectionProvenance(
        dataset_role="DESIGN",
        method="manual_fixture",
        source="unit-test",
        details=("representative seed",),
    )
    return PanelDefinition(
        target=PanelTarget(
            name="Target virus",
            taxid=1001,
            mode="broad_detection",
            subtype=None,
            groups=(
                TargetGroup(
                    name="group-a",
                    required=True,
                    dataset_roles=("DESIGN", "CHALLENGE"),
                    reasons=("target_diversity",),
                    proposed_by=("manual",),
                    sequence_selection=(selection,),
                ),
            ),
        ),
        non_targets=(
            PanelNonTarget(
                name="Neighbor virus",
                taxid=2001,
                criticality="CRITICAL",
                dataset_roles=("DESIGN", "CHALLENGE"),
                reasons=("phylogenetic_neighbor",),
                proposed_by=("manual",),
                sequence_selection=(),
            ),
        ),
        diagnostic_context=DiagnosticContext(
            syndrome="febrile illness",
            geography="test-region",
            sample_type="serum",
            vector="mosquito",
        ),
    )


def test_valid_panel_definition_is_accepted():
    validate_panel_definition(valid_panel())


def test_subtype_specific_requires_subtype():
    panel = valid_panel()
    invalid = PanelDefinition(
        target=PanelTarget(
            name=panel.target.name,
            taxid=panel.target.taxid,
            mode="subtype_specific",
            subtype=None,
            groups=panel.target.groups,
        ),
        non_targets=panel.non_targets,
        diagnostic_context=panel.diagnostic_context,
    )
    with pytest.raises(ValueError, match="subtype_specific.*subtype"):
        validate_panel_definition(invalid)


def test_broad_detection_rejects_subtype():
    panel = valid_panel()
    invalid = PanelDefinition(
        target=PanelTarget(
            name=panel.target.name,
            taxid=panel.target.taxid,
            mode="broad_detection",
            subtype="group-a",
            groups=panel.target.groups,
        ),
        non_targets=panel.non_targets,
        diagnostic_context=panel.diagnostic_context,
    )
    with pytest.raises(ValueError, match="broad_detection.*subtype"):
        validate_panel_definition(invalid)


def test_non_target_cannot_duplicate_target_taxid():
    panel = valid_panel()
    duplicate = PanelNonTarget(
        name="Target alias",
        taxid=1001,
        criticality="CRITICAL",
        dataset_roles=("DESIGN",),
        reasons=("test",),
        proposed_by=("manual",),
        sequence_selection=(),
    )
    invalid = PanelDefinition(
        target=panel.target,
        non_targets=(duplicate,),
        diagnostic_context=panel.diagnostic_context,
    )
    with pytest.raises(ValueError, match="target.*non-target.*TaxID"):
        validate_panel_definition(invalid)


def test_dataset_roles_must_be_unique():
    panel = valid_panel()
    invalid_group = TargetGroup(
        name="group-a",
        required=True,
        dataset_roles=("DESIGN", "DESIGN"),
        reasons=("target_diversity",),
        proposed_by=("manual",),
        sequence_selection=(),
    )
    invalid = PanelDefinition(
        target=PanelTarget(
            name=panel.target.name,
            taxid=panel.target.taxid,
            mode=panel.target.mode,
            subtype=None,
            groups=(invalid_group,),
        ),
        non_targets=panel.non_targets,
        diagnostic_context=panel.diagnostic_context,
    )
    with pytest.raises(ValueError, match="dataset_roles.*unique"):
        validate_panel_definition(invalid)
```

Add focused tests for blank names, non-positive TaxID, duplicate target-group names using casefolding, duplicate non-target names using casefolding, duplicate non-null non-target TaxIDs, invalid criticality, invalid dataset role, empty dataset-role tuple, blank reason/proposer values, invalid optional diagnostic-context strings, and invalid `SequenceSelectionProvenance` fields.

- [ ] **Step 2: Run the test and verify the expected RED state**

```bash
python -m pytest tests/test_panel.py -v
```

Expected: test collection fails with `ModuleNotFoundError: No module named 'qpcr_pipeline.panel'`.

- [ ] **Step 3: Implement `qpcr_pipeline/panel.py`**

Create the file with this complete public model and validation implementation:

```python
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
            raise ValueError(f"Panel target group {index} required must be a boolean.")
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
                raise ValueError("Panel target and non-target TaxIDs must be distinct.")
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
```

- [ ] **Step 4: Add the manual West Nile review fixture**

Create `tests/fixtures/panels/west_nile_proposal.yaml`:

```yaml
schema_version: 1
status: PROPOSED
definition:
  target:
    name: West Nile virus
    taxid: null
    mode: broad_detection
    subtype: null
    groups:
      - name: lineage_1
        required: true
        dataset_roles: [DESIGN, CHALLENGE]
        reasons: [target_diversity]
        proposed_by: [manual]
        sequence_selection: []
      - name: lineage_2
        required: true
        dataset_roles: [DESIGN, CHALLENGE]
        reasons: [target_diversity]
        proposed_by: [manual]
        sequence_selection: []
  non_targets:
    - name: Usutu virus
      taxid: null
      criticality: CRITICAL
      dataset_roles: [DESIGN, CHALLENGE]
      reasons: [phylogenetic_neighbor]
      proposed_by: [manual]
      sequence_selection: []
    - name: Japanese encephalitis virus
      taxid: null
      criticality: CRITICAL
      dataset_roles: [DESIGN, CHALLENGE]
      reasons: [phylogenetic_neighbor]
      proposed_by: [manual]
      sequence_selection: []
    - name: Dengue virus
      taxid: null
      criticality: IMPORTANT
      dataset_roles: [CHALLENGE]
      reasons: [clinical_differential]
      proposed_by: [manual]
      sequence_selection: []
  diagnostic_context:
    syndrome: arboviral febrile disease
    geography: Brazil
    sample_type: human serum
    vector: mosquito
```

The null TaxIDs are deliberate in this test-only fixture: subproject 1 validates the identifier field but does not yet implement taxonomy resolution.

- [ ] **Step 5: Run the domain tests**

```bash
python -m pytest tests/test_panel.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit the domain layer**

```bash
git add qpcr_pipeline/panel.py tests/test_panel.py tests/fixtures/panels/west_nile_proposal.yaml
git commit -m "feat: add panel domain model"
```

---

### Task 2: Add strict proposal and approved-manifest serialization

**Files:**
- Create: `qpcr_pipeline/panel_manifest.py`
- Create: `tests/test_panel_manifest.py`

**Interfaces:**
- Consumes `PanelDefinition` and `validate_panel_definition()` from Task 1.
- Produces `PanelProposal`, `ApprovedPanelManifest`, `PanelResult`, `load_panel_proposal()`, `write_panel_proposal()`, `approve_panel_proposal()`, `load_approved_panel_manifest()`, `write_approved_panel_manifest()`, `proposal_semantic_sha256()`, and `materialize_approved_panel()`.

- [ ] **Step 1: Write failing manifest tests**

Create `tests/test_panel_manifest.py`:

```python
import json
from pathlib import Path

import pytest

from qpcr_pipeline.panel_manifest import (
    approve_panel_proposal,
    load_approved_panel_manifest,
    load_panel_proposal,
    proposal_semantic_sha256,
)

FIXTURE = Path("tests/fixtures/panels/west_nile_proposal.yaml")


def test_load_proposal_validates_schema_and_domain():
    proposal = load_panel_proposal(FIXTURE)
    assert proposal.schema_version == 1
    assert proposal.status == "PROPOSED"
    assert proposal.definition.target.name == "West Nile virus"
    assert proposal.definition.non_targets[0].criticality == "CRITICAL"


def test_approval_is_byte_deterministic(tmp_path):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    approve_panel_proposal(FIXTURE, first)
    approve_panel_proposal(FIXTURE, second)
    assert first.read_bytes() == second.read_bytes()


def test_semantic_hash_ignores_yaml_comments_and_formatting(tmp_path):
    original = load_panel_proposal(FIXTURE)
    reformatted = tmp_path / "reformatted.yaml"
    reformatted.write_text(
        "# review comment\n" + FIXTURE.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    assert proposal_semantic_sha256(original) == proposal_semantic_sha256(
        load_panel_proposal(reformatted)
    )


def test_approved_manifest_rejects_extra_fields(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["unexpected"] = True
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved panel manifest fields"):
        load_approved_panel_manifest(approved)


def test_approved_manifest_rejects_invalid_proposal_hash(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["proposal_sha256"] = "sha256:not-a-digest"
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="proposal_sha256"):
        load_approved_panel_manifest(approved)
```

Add tests for wrong schema version, wrong status, `approved_by_user` not exactly `True`, missing fields, invalid nested panel content, and atomic-write preservation when serialization raises before replacement. Do not add duplicate-YAML-key detection in this subproject; PyYAML's loader behavior is outside this manifest contract.

- [ ] **Step 2: Run the test and verify RED**

```bash
python -m pytest tests/test_panel_manifest.py -v
```

Expected: collection fails because `qpcr_pipeline.panel_manifest` does not exist.

- [ ] **Step 3: Implement strict payload conversion and manifest types**

Create `qpcr_pipeline/panel_manifest.py`. The public dataclasses are:

```python
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
```

Implement a strict `_require_fields(raw, expected, label)` helper that requires `raw` to be a dict and `set(raw) == expected`. Implement explicit parsers for every nested mapping; reuse the dataclass names from Task 1 and call `validate_panel_definition()` after construction. Use these exact allowed-field sets:

```python
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
```

Convert YAML/JSON lists to tuples only in the explicit parser. Reject non-list values instead of accepting arbitrary iterables.

- [ ] **Step 4: Implement canonical hashing and deterministic output**

Use this exact canonicalization policy:

```python
_SHA256_PATTERN = re.compile(r"^sha256:[0-9a-f]{64}$")


def _canonical_json_text(payload: object) -> str:
    return json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
    ) + "\n"


def proposal_semantic_sha256(proposal: PanelProposal) -> str:
    text = _canonical_json_text(_proposal_payload(proposal))
    digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return "sha256:" + digest
```

`_proposal_payload()` and `_definition_payload()` must emit ordinary dict/list/scalar structures in the same field semantics accepted by the strict parser. Approved files use:

```python
def _approved_text(manifest: ApprovedPanelManifest) -> str:
    return json.dumps(
        _approved_payload(manifest),
        indent=2,
        sort_keys=True,
        ensure_ascii=False,
    ) + "\n"
```

No timestamp is stored in the approved manifest.

- [ ] **Step 5: Implement atomic writes, approval, and materialization**

Use a shared helper:

```python
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
```

Approval logic is exactly:

```python
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
```

`write_panel_proposal()` normalizes YAML with `yaml.safe_dump(payload, sort_keys=False, allow_unicode=True)` and one final newline. `materialize_approved_panel()` must:

1. load and validate the source manifest;
2. write canonical approved JSON to `outdir / "panel" / "approved_panel.json"`;
3. compute SHA-256 over the materialized file bytes;
4. return:

```python
PanelResult(
    status="APPROVED",
    manifest_sha256="sha256:" + digest,
    manifest_path=destination,
    target_mode=manifest.definition.target.mode,
    non_target_count=len(manifest.definition.non_targets),
)
```

- [ ] **Step 6: Run manifest tests**

```bash
python -m pytest tests/test_panel.py tests/test_panel_manifest.py -v
```

Expected: PASS.

- [ ] **Step 7: Commit manifest I/O**

```bash
git add qpcr_pipeline/panel_manifest.py tests/test_panel_manifest.py
git commit -m "feat: add deterministic panel manifests"
```

---

### Task 3: Add panel configuration and the primer-design requirement

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Create: `tests/test_panel_config.py`
- Modify relevant existing config/pipeline fixtures discovered by the grep command in this task.

**Interfaces:**
- Consumes `PanelDefinition` from Task 1.
- Produces `PanelConfig` and `PipelineConfig.panel`.

- [ ] **Step 1: Write failing config tests**

Create `tests/test_panel_config.py`. Cover inline proposal and frozen-manifest modes:

```python
from pathlib import Path

import pytest

from qpcr_pipeline.config import load_config


def test_loads_inline_panel_proposal(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        f"""
target:
  name: West Nile virus
input:
  fasta: {fasta.as_posix()}
panel:
  proposal:
    target:
      name: West Nile virus
      taxid: null
      mode: broad_detection
      subtype: null
      groups:
        - name: lineage_1
          required: true
          dataset_roles: [DESIGN, CHALLENGE]
          reasons: [target_diversity]
          proposed_by: [manual]
          sequence_selection: []
    non_targets:
      - name: Usutu virus
        taxid: null
        criticality: CRITICAL
        dataset_roles: [DESIGN, CHALLENGE]
        reasons: [phylogenetic_neighbor]
        proposed_by: [manual]
        sequence_selection: []
    diagnostic_context:
      syndrome: arboviral febrile disease
      geography: Brazil
      sample_type: human serum
      vector: mosquito
""",
        encoding="utf-8",
    )
    config = load_config(config_path)
    assert config.panel is not None
    assert config.panel.proposal is not None
    assert config.panel.frozen_manifest is None
    assert config.panel.proposal.target.mode == "broad_detection"


def test_primer_design_requires_panel(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        f"""
target:
  name: target
input:
  fasta: {fasta.as_posix()}
alignment:
  enabled: true
conservation:
  enabled: true
primer_design:
  enabled: true
""",
        encoding="utf-8",
    )
    with pytest.raises(
        ValueError,
        match="Enabled primer design requires a panel proposal or frozen manifest",
    ):
        load_config(config_path)
```

Also test: `panel.frozen_manifest` parses to `Path`; panel config rejects both proposal and frozen manifest together; rejects neither when a panel section is present; rejects unknown fields at every nested level; inline target name must equal top-level target name after `.strip()`; panel-free primer-design-disabled config remains valid.

- [ ] **Step 2: Run and verify RED**

```bash
python -m pytest tests/test_panel_config.py -v
```

Expected: FAIL because `PanelConfig` and the parser do not exist.

- [ ] **Step 3: Add `PanelConfig` and `PipelineConfig.panel`**

In `qpcr_pipeline/config.py`, import the panel model types and define:

```python
@dataclass(frozen=True, slots=True)
class PanelConfig:
    proposal: PanelDefinition | None = None
    frozen_manifest: Path | None = None
```

Insert `panel: PanelConfig | None = None` after the three input fields in `PipelineConfig`. In `load_config()`, parse `panel_config = _parse_panel_config(raw.get("panel"))` and pass `panel=panel_config` to `PipelineConfig`.

- [ ] **Step 4: Implement strict panel config parsing**

Add these complete parser rules:

```python
def _parse_panel_config(raw: Any) -> PanelConfig | None:
    if raw is None:
        return None
    if not isinstance(raw, dict):
        raise ValueError("Configuration section 'panel' must be a mapping.")
    allowed = {"proposal", "frozen_manifest"}
    unknown = set(raw) - allowed
    if unknown:
        rendered = ", ".join(sorted(str(item) for item in unknown))
        raise ValueError(
            f"Configuration section 'panel' fields {rendered} are unrecognized."
        )
    proposal = (
        _parse_panel_definition(raw["proposal"])
        if "proposal" in raw
        else None
    )
    frozen = _optional_path(raw, "frozen_manifest")
    config = PanelConfig(proposal=proposal, frozen_manifest=frozen)
    if sum(value is not None for value in (proposal, frozen)) != 1:
        raise ValueError(
            "Panel configuration must specify exactly one proposal or frozen_manifest."
        )
    return config
```

Implement `_parse_panel_definition`, `_parse_panel_target`, `_parse_target_group`, `_parse_panel_non_target`, `_parse_diagnostic_context`, and `_parse_sequence_selection` using the same exact allowed-field sets from Task 2. Each list-valued YAML field must be checked as a `list`, then converted to tuple. After constructing `PanelDefinition`, call `validate_panel_definition()`.

Do not import `panel_manifest.py` into `config.py`; frozen-manifest file validation belongs to runtime preflight, which avoids a config/manifest import cycle.

- [ ] **Step 5: Add config-level consistency validation**

Add:

```python
def validate_panel_config(config: PanelConfig, *, target_name: str) -> None:
    if not isinstance(config, PanelConfig):
        raise ValueError("Panel configuration must be a PanelConfig.")
    if sum(
        value is not None
        for value in (config.proposal, config.frozen_manifest)
    ) != 1:
        raise ValueError(
            "Panel configuration must specify exactly one proposal or frozen_manifest."
        )
    if config.proposal is not None:
        validate_panel_definition(config.proposal)
        if config.proposal.target.name.strip() != target_name.strip():
            raise ValueError(
                "Inline panel target name must match pipeline target_name."
            )
    if config.frozen_manifest is not None and not isinstance(
        config.frozen_manifest,
        Path,
    ):
        raise ValueError("Panel frozen_manifest must be a Path when configured.")
```

In `validate_pipeline_config()` call it when `config.panel` is not `None`, then add:

```python
if config.primer_design.enabled and config.panel is None:
    raise ValueError(
        "Enabled primer design requires a panel proposal or frozen manifest."
    )
```

- [ ] **Step 6: Repair existing enabled-primer pipeline/config fixtures**

Run:

```bash
grep -R "PrimerDesignConfig(enabled=True\|primer_design:" -n tests
```

For full `PipelineConfig` tests whose purpose is not panel validation, add an explicit minimal `PanelConfig(proposal=...)` helper. For tests that exercise `primer_design.py` directly without `PipelineConfig`, make no change. Do not weaken the new validation.

- [ ] **Step 7: Run config-focused tests**

```bash
python -m pytest \
  tests/test_panel_config.py \
  tests/test_config.py \
  tests/test_specificity_config.py \
  tests/test_ranking_config.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit config integration**

```bash
git add qpcr_pipeline/config.py tests/test_panel_config.py tests/test_config.py tests/test_minimal_run.py tests/test_pipeline_ranking.py tests/test_pipeline_specificity.py
git commit -m "feat: configure assay panels"
```

Only stage the existing test files from that command if they actually changed.

---

### Task 4: Add proposal preflight and `ACTION_REQUIRED` lifecycle

**Files:**
- Modify: `qpcr_pipeline/panel_manifest.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_pipeline_panel.py`
- Modify: `tests/test_run_manifest.py`
- Modify: `tests/test_run_logging.py`

**Interfaces:**
- Produces `PanelPreflight`, `prepare_panel_preflight(panel, outdir, *, target_name)`, `RunRecorder.action_required()`, `RunSummary.action_required_code`, and `RunSummary.action_required_artifact`.

- [ ] **Step 1: Write failing first-run gate test**

In `tests/test_pipeline_panel.py`, define a helper that uses a small local FASTA and the parsed fixture definition, then test:

```python
def test_panel_proposal_stops_before_input_and_writes_review_artifact(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    summary = run_pipeline(config, outdir)
    assert summary.status == "ACTION_REQUIRED"
    assert summary.action_required_code == "PANEL_APPROVAL_REQUIRED"
    assert summary.sequence_count == 0
    assert summary.sequence_ids == []
    assert (outdir / "panel_proposal.yaml").is_file()
    assert not (outdir / ".checkpoints" / "input" / "manifest.json").exists()
    assert not (outdir / "qc_report.json").exists()
```

Add run-manifest assertions that status is `ACTION_REQUIRED`, the action code is `PANEL_APPROVAL_REQUIRED`, the artifact path equals the proposal path, and the first attempt is also `ACTION_REQUIRED` rather than `FAILED`.

- [ ] **Step 2: Run and verify RED**

```bash
python -m pytest tests/test_pipeline_panel.py::test_panel_proposal_stops_before_input_and_writes_review_artifact -v
```

Expected: FAIL because proposal preflight/action-required lifecycle is absent.

- [ ] **Step 3: Add `PanelPreflight` with one consistent signature**

In `panel_manifest.py` add:

```python
@dataclass(frozen=True, slots=True)
class PanelPreflight:
    status: Literal["READY", "ACTION_REQUIRED", "LEGACY"]
    proposal_path: Path | None = None
    approved_manifest_path: Path | None = None


def prepare_panel_preflight(
    panel: PanelConfig | None,
    outdir: Path,
    *,
    target_name: str,
) -> PanelPreflight:
    if panel is None:
        return PanelPreflight(status="LEGACY")
    if panel.proposal is not None:
        if panel.proposal.target.name.strip() != target_name.strip():
            raise ValueError(
                f"Panel target {panel.proposal.target.name!r} does not match "
                f"pipeline target {target_name!r}."
            )
        proposal_path = Path(outdir) / "panel_proposal.yaml"
        write_panel_proposal(panel.proposal, proposal_path)
        return PanelPreflight(
            status="ACTION_REQUIRED",
            proposal_path=proposal_path,
        )
    if panel.frozen_manifest is None:
        raise ValueError("Panel preflight requires a proposal or frozen manifest.")
    manifest = load_approved_panel_manifest(panel.frozen_manifest)
    if manifest.definition.target.name.strip() != target_name.strip():
        raise ValueError(
            f"Approved panel target {manifest.definition.target.name!r} does not "
            f"match pipeline target {target_name!r}."
        )
    return PanelPreflight(
        status="READY",
        approved_manifest_path=panel.frozen_manifest,
    )
```

Use `TYPE_CHECKING` for the `PanelConfig` annotation and a quoted annotation if needed, so `config.py` never imports `panel_manifest.py` and no runtime cycle is introduced.

- [ ] **Step 4: Extend run recording without bumping schema version**

In `EVENT_FIELDS`, add:

```python
"run_action_required": frozenset({"status", "code", "artifact"}),
```

When creating a new run manifest, add:

```python
"action_required": None,
"panel_provenance": {},
```

When loading an existing schema-1 manifest in `begin_attempt()`, normalize old manifests before any access:

```python
payload.setdefault("action_required", None)
payload.setdefault("panel_provenance", {})
payload["action_required"] = None
```

Add this method:

```python
def action_required(self, code: str, artifact: Path) -> None:
    if not isinstance(code, str) or not code:
        raise ValueError("Action-required code must be non-empty.")
    if self._payload is None:
        raise RuntimeError("No active run manifest.")
    now = self.clock()
    action = {
        "code": code,
        "artifact": sanitize_diagnostic(artifact),
    }
    attempt = self._active_attempt()
    attempt["status"] = "ACTION_REQUIRED"
    attempt["finished_at"] = now
    attempt["failure"] = None
    self._payload["status"] = "ACTION_REQUIRED"
    self._payload["updated_at"] = now
    self._payload["action_required"] = action
    self._payload["failure"] = None
    self._write()
    assert self._logger is not None
    self._logger.emit(
        "run_action_required",
        status="ACTION_REQUIRED",
        code=code,
        artifact=artifact,
    )
```

- [ ] **Step 5: Add action-required preflight before normal planning**

Extend `RunSummary` in `pipeline.py`:

```python
@dataclass(frozen=True, slots=True)
class RunSummary:
    status: str
    target_name: str
    sequence_count: int
    sequence_ids: list[str]
    stage_actions: tuple[StageActionSummary, ...] = ()
    action_required_code: str | None = None
    action_required_artifact: str | None = None
```

Keep the existing environment inspection. Immediately afterward call:

```python
panel_preflight = prepare_panel_preflight(
    config.panel,
    output_dir,
    target_name=config.target_name,
)
```

If `panel_preflight.status == "ACTION_REQUIRED"`, do not call `plan_pipeline()`. Create the recorder and call `begin_attempt()` with a synthetic plan:

```python
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
return summary
```

For `READY` or `LEGACY`, continue into existing `plan_pipeline()` and run logic.

- [ ] **Step 6: Test repeatability and old-manifest compatibility**

Add tests that:

- Two runs with the same inline proposal write byte-identical proposal YAML.
- A stale prior `run_summary.json` is replaced by the action-required summary.
- A schema-1 run manifest created before these fields existed can be resumed because `setdefault()` adds `action_required` and `panel_provenance`.
- `run.log.jsonl` includes exactly one `run_action_required` event for the action-required attempt.

- [ ] **Step 7: Run lifecycle tests**

```bash
python -m pytest \
  tests/test_pipeline_panel.py \
  tests/test_run_manifest.py \
  tests/test_run_recording.py \
  tests/test_run_logging.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit the approval gate**

```bash
git add qpcr_pipeline/panel_manifest.py qpcr_pipeline/pipeline.py qpcr_pipeline/run_recording.py tests/test_pipeline_panel.py tests/test_run_manifest.py tests/test_run_recording.py tests/test_run_logging.py
git commit -m "feat: gate assay design on panel approval"
```

---

### Task 5: Add the approved panel as the first checkpoint stage

**Files:**
- Modify: `qpcr_pipeline/execution.py`
- Modify: `qpcr_pipeline/checkpoint_codecs.py`
- Modify: `qpcr_pipeline/checkpoint_stages.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_execution.py`
- Modify: `tests/test_execution_plan.py`
- Modify: `tests/test_checkpoint_codecs.py`
- Modify: `tests/test_checkpoint_stages.py`
- Modify: `tests/test_pipeline_resume.py`

**Interfaces:**
- Consumes `PanelResult` and `materialize_approved_panel()`.
- Produces stage `panel`, always before `input`.

- [ ] **Step 1: Write failing execution-graph tests**

Update `tests/test_execution.py`:

```python
from qpcr_pipeline.execution import (
    STAGE_DEPENDENCIES,
    STAGE_ORDER,
    ExecutionPolicy,
    plan_from_validity,
    required_reuse_boundary,
    transitive_descendants,
)


def test_panel_is_first_stage_and_input_depends_on_it():
    assert STAGE_ORDER[0] == "panel"
    assert STAGE_DEPENDENCIES["panel"] == ()
    assert STAGE_DEPENDENCIES["input"] == ("panel",)


def test_resume_invalid_panel_forces_every_downstream_stage():
    reusable = {stage: True for stage in STAGE_ORDER}
    reusable["panel"] = False
    actions = {
        item.stage: item.action
        for item in plan_from_validity(ExecutionPolicy(resume=True), reusable)
    }
    assert actions["panel"] == "RUN"
    assert all(actions[stage] == "FORCED" for stage in STAGE_ORDER[1:])
```

Update literal action dictionaries in this file and `tests/test_execution_plan.py` to include `"panel": "REUSE"`, `"RUN"`, or `"FORCED"` according to the dependency under test.

- [ ] **Step 2: Run execution tests and verify RED**

```bash
python -m pytest tests/test_execution.py tests/test_execution_plan.py -v
```

Expected: FAIL because `panel` is absent from the graph.

- [ ] **Step 3: Modify `execution.py`**

Use exactly this order and dependency prefix:

```python
StageName = Literal[
    "panel",
    "input",
    "qc",
    "clustering",
    "alignment",
    "conservation",
    "primer_design",
    "inclusivity",
    "specificity",
    "ranking",
]

STAGE_ORDER: tuple[StageName, ...] = (
    "panel",
    "input",
    "qc",
    "clustering",
    "alignment",
    "conservation",
    "primer_design",
    "inclusivity",
    "specificity",
    "ranking",
)

STAGE_DEPENDENCIES: dict[StageName, tuple[StageName, ...]] = {
    "panel": (),
    "input": ("panel",),
    "qc": ("input",),
    "clustering": ("qc",),
    "alignment": ("clustering",),
    "conservation": ("alignment",),
    "primer_design": ("conservation",),
    "inclusivity": ("primer_design", "qc"),
    "specificity": ("primer_design",),
    "ranking": ("primer_design", "inclusivity", "specificity"),
}
```

- [ ] **Step 4: Add `PANEL_CODEC` and round-trip tests**

In `checkpoint_codecs.py` import `PanelResult` and add:

```python
PANEL_CODEC = _StructuralCodec[PanelResult](PanelResult)
```

Add tests for:

```python
legacy = PanelResult(
    status="LEGACY",
    manifest_sha256=None,
    manifest_path=None,
    target_mode=None,
    non_target_count=0,
)
```

and an approved `PanelResult` whose manifest path is inside the test output directory. Verify encode/decode equality.

- [ ] **Step 5: Add panel stage definition and checkpoint identity**

In `checkpoint_stages.py`:

1. import `PANEL_CODEC`;
2. add `"panel": PANEL_CODEC` to `_CODECS`;
3. add to `stage_parameters()` before `input`:

```python
if stage == "panel":
    if config.panel is None:
        return {"mode": "LEGACY"}
    return {"mode": "APPROVED_MANIFEST"}
```

4. add to `stage_input_identities()` before `input`:

```python
if stage == "panel":
    if config.panel is None:
        return {}
    if config.panel.frozen_manifest is None:
        raise ValueError("Panel proposal cannot enter checkpoint planning.")
    load_approved_panel_manifest(config.panel.frozen_manifest)
    return {
        "approved_panel": {
            "sha256": file_sha256(config.panel.frozen_manifest),
        }
    }
```

5. add to `stage_outputs()`:

```python
if stage == "panel":
    values = _paths(getattr(result, "manifest_path", None))
elif stage == "input":
    candidates = (
        outdir / "ncbi_dataset" / "records.gb",
        outdir / "ncbi_dataset" / "dataset_manifest.json",
        outdir / "ncbi_dataset_manifest.json",
    )
    values = tuple(path for path in candidates if path.is_file())
```

Leave panel tool identities empty.

- [ ] **Step 6: Execute the panel stage in `pipeline.py`**

At the top of `_run_stage()` add:

```python
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
    return materialize_approved_panel(
        config.panel.frozen_manifest,
        output_dir,
    )
```

Only after this branch should `_run_stage()` read `results["input"]`.

- [ ] **Step 7: Test checkpoint invalidation**

In `tests/test_pipeline_resume.py`, build two approved panel files from two proposal files that differ only in Usutu criticality (`CRITICAL` vs `IMPORTANT`). Use a local input and keep primer design disabled so this checkpoint test needs no Primer3 binary.

Assert:

- first run creates `.checkpoints/panel/manifest.json`;
- unchanged `--resume` returns `REUSE` for all stages;
- changing `frozen_manifest` to the scientifically changed approved file returns `RUN` for panel and `FORCED` for every stage after panel.

This conservative invalidation is intentional until panel-driven multi-dataset acquisition exists.

- [ ] **Step 8: Run checkpoint/resume suites**

```bash
python -m pytest \
  tests/test_execution.py \
  tests/test_execution_plan.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py -v
```

Expected: PASS.

- [ ] **Step 9: Commit checkpoint integration**

```bash
git add qpcr_pipeline/execution.py qpcr_pipeline/checkpoint_codecs.py qpcr_pipeline/checkpoint_stages.py qpcr_pipeline/pipeline.py tests/test_execution.py tests/test_execution_plan.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py tests/test_pipeline_resume.py
git commit -m "feat: checkpoint approved assay panels"
```

---

### Task 6: Add panel provenance and the explicit approval CLI

**Files:**
- Modify: `qpcr_pipeline/provenance.py`
- Modify: `qpcr_pipeline/run_recording.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `qpcr_pipeline/cli.py`
- Create: `tests/test_panel_cli.py`
- Modify: `tests/test_cli.py`
- Modify: `tests/test_run_manifest.py`
- Modify: `tests/test_run_recording.py`

**Interfaces:**
- Produces `build_panel_provenance()` and `qpcr-pipeline panel approve`.
- CLI exit code 3 denotes expected `ACTION_REQUIRED`, not technical failure.

- [ ] **Step 1: Write failing CLI tests**

Create `tests/test_panel_cli.py`:

```python
import json
import sys
from pathlib import Path

from qpcr_pipeline.cli import main


def test_panel_approve_command_writes_frozen_manifest(
    tmp_path,
    monkeypatch,
    capsys,
):
    proposal = Path("tests/fixtures/panels/west_nile_proposal.yaml")
    approved = tmp_path / "approved.json"
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "qpcr-pipeline",
            "panel",
            "approve",
            str(proposal),
            "--output",
            str(approved),
        ],
    )
    assert main() == 0
    payload = json.loads(approved.read_text(encoding="utf-8"))
    assert payload["status"] == "APPROVED"
    assert payload["approved_by_user"] is True
    assert payload["proposal_sha256"].startswith("sha256:")
    assert "Approved panel manifest:" in capsys.readouterr().out
```

Add a CLI-run test that monkeypatches `run_pipeline()` to return an action-required `RunSummary`, then asserts `main() == 3` and output contains both `PANEL_APPROVAL_REQUIRED` and the proposal path.

- [ ] **Step 2: Run and verify RED**

```bash
python -m pytest tests/test_panel_cli.py -v
```

Expected: FAIL because the nested command is absent.

- [ ] **Step 3: Add nested `panel approve` command**

In `build_parser()` add:

```python
panel_parser = subparsers.add_parser(
    "panel",
    help="Review and freeze assay panels",
)
panel_subparsers = panel_parser.add_subparsers(
    dest="panel_command",
    required=True,
)
approve_parser = panel_subparsers.add_parser(
    "approve",
    help="Approve a reviewed panel proposal and write a frozen manifest",
)
approve_parser.add_argument("proposal", type=Path)
approve_parser.add_argument("--output", type=Path, required=True)
```

Before the `run` branch in `main()` add:

```python
if args.command == "panel" and args.panel_command == "approve":
    manifest = approve_panel_proposal(args.proposal, args.output)
    print(f"Approved panel manifest: {args.output}")
    print(f"Proposal identity: {manifest.proposal_sha256}")
    return 0
```

After `summary = run_pipeline(...)`, before the normal status print, add:

```python
if summary.status == "ACTION_REQUIRED":
    print(
        f"ACTION_REQUIRED: {summary.action_required_code} "
        f"({summary.action_required_artifact})"
    )
    return 3
```

- [ ] **Step 4: Add safe panel provenance**

In `provenance.py` add:

```python
def build_panel_provenance(panel_result: object) -> dict[str, object]:
    status = getattr(panel_result, "status", None)
    if status == "LEGACY":
        return {"mode": "legacy_unconfigured"}
    if status != "APPROVED":
        raise ValueError(
            "Panel provenance requires an approved or legacy panel result."
        )
    manifest_sha256 = _sha256_identity(
        getattr(panel_result, "manifest_sha256", None)
    )
    if manifest_sha256 is None:
        raise ValueError(
            "Approved panel provenance is missing a SHA-256 identity."
        )
    target_mode = getattr(panel_result, "target_mode", None)
    if target_mode not in {"broad_detection", "subtype_specific"}:
        raise ValueError("Approved panel provenance target_mode is invalid.")
    non_target_count = getattr(panel_result, "non_target_count", None)
    if (
        isinstance(non_target_count, bool)
        or not isinstance(non_target_count, int)
        or non_target_count < 0
    ):
        raise ValueError("Approved panel provenance non_target_count is invalid.")
    return {
        "mode": "approved_manifest",
        "manifest_sha256": manifest_sha256,
        "target_mode": target_mode,
        "non_target_count": non_target_count,
    }
```

Do not embed the whole panel in `run_manifest.json`; the canonical frozen manifest is the detailed artifact.

- [ ] **Step 5: Persist panel provenance on successful execution**

Change `RunRecorder.complete()` to use keyword-only scientific provenance arguments so future additions do not create positional ambiguity:

```python
def complete(
    self,
    status: str,
    scientific_completeness: ScientificCompleteness | dict[str, object],
    *,
    input_provenance: object,
    reference: object,
    panel_provenance: object,
) -> None:
```

Inside it set:

```python
self._payload["input_provenance"] = sanitize_diagnostic(input_provenance)
self._payload["reference"] = sanitize_diagnostic(reference)
self._payload["panel_provenance"] = sanitize_diagnostic(panel_provenance)
self._payload["action_required"] = None
```

Update the pipeline call:

```python
recorder.complete(
    final_status,
    final_completeness,
    input_provenance=build_input_provenance(
        config,
        output_dir,
        qc_result,
        manifests["input"],
    ),
    reference=build_reference_provenance(alignment),
    panel_provenance=build_panel_provenance(results["panel"]),
)
```

Update direct `RunRecorder.complete()` tests to use these keywords.

- [ ] **Step 6: Add provenance assertions**

For an approved-panel run in `tests/test_run_manifest.py`:

```python
provenance = manifest["panel_provenance"]
assert provenance["mode"] == "approved_manifest"
assert provenance["manifest_sha256"].startswith("sha256:")
assert provenance["target_mode"] == "broad_detection"
assert provenance["non_target_count"] == 3
```

For a legacy panel-free non-design run:

```python
assert manifest["panel_provenance"] == {"mode": "legacy_unconfigured"}
```

Also assert no sequence string or complete panel definition is copied into the panel provenance summary.

- [ ] **Step 7: Run CLI/provenance tests**

```bash
python -m pytest \
  tests/test_panel_cli.py \
  tests/test_cli.py \
  tests/test_run_manifest.py \
  tests/test_run_recording.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit CLI and provenance**

```bash
git add qpcr_pipeline/provenance.py qpcr_pipeline/run_recording.py qpcr_pipeline/pipeline.py qpcr_pipeline/cli.py tests/test_panel_cli.py tests/test_cli.py tests/test_run_manifest.py tests/test_run_recording.py
git commit -m "feat: approve and record assay panels"
```

---

### Task 7: Make dry-run non-mutating and verify the complete panel workflow

**Files:**
- Modify: `qpcr_pipeline/dry_run.py`
- Modify: `qpcr_pipeline/cli.py`
- Modify: `tests/test_dry_run.py`
- Modify: `tests/test_pipeline_panel.py`
- Modify: `README.md`
- Modify any remaining enabled-primer pipeline fixtures identified by the full suite.

**Interfaces:**
- Consumes all prior tasks.
- Produces a complete offline proposal → approve → frozen rerun workflow.

- [ ] **Step 1: Write dry-run tests before modifying dry-run implementation**

Add tests in `tests/test_dry_run.py`:

```python
def test_dry_run_proposal_reports_action_required_without_writing(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    report = dry_run_pipeline(config, outdir)
    assert report.panel_action_required is True
    assert report.panel_proposal_would_be_written == str(
        outdir / "panel_proposal.yaml"
    )
    assert not outdir.exists()


def test_dry_run_with_frozen_panel_starts_with_panel_stage(tmp_path):
    config = frozen_panel_pipeline_config(tmp_path)
    report = dry_run_pipeline(config, tmp_path / "run")
    assert report.panel_action_required is False
    assert report.decisions[0].stage == "panel"
```

- [ ] **Step 2: Run and verify RED**

```bash
python -m pytest tests/test_dry_run.py -v
```

Expected: FAIL because dry-run has no panel preflight fields.

- [ ] **Step 3: Extend `DryRunReport` and preflight before planning**

In `dry_run.py`, extend the report dataclass with:

```python
panel_action_required: bool = False
panel_proposal_would_be_written: str | None = None
```

Do not call the mutating `prepare_panel_preflight()` for inline proposals. Instead inspect validated config directly:

```python
if config.panel is not None and config.panel.proposal is not None:
    proposal_path = None if outdir is None else str(Path(outdir) / "panel_proposal.yaml")
    return DryRunReport(
        target_name=config.target_name,
        environment=environment,
        decisions=(),
        panel_action_required=True,
        panel_proposal_would_be_written=proposal_path,
    )
```

For frozen/legacy panel state, call `plan_pipeline()` normally. Frozen manifests are read-only validated through checkpoint request creation; no output file may be written in dry-run.

Update `_render_dry_run()` in `cli.py` so action-required reports include:

```text
Panel approval required before scientific execution.
Proposal would be written to: <path>
```

Only print the second line when a path is available.

- [ ] **Step 4: Add an offline proposal → approve → rerun integration test**

In `tests/test_pipeline_panel.py`:

```python
def test_proposal_approve_rerun_reaches_checkpointed_pipeline(tmp_path):
    proposal_config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"
    first = run_pipeline(proposal_config, outdir)
    assert first.status == "ACTION_REQUIRED"

    approved = tmp_path / "approved_panel.json"
    approve_panel_proposal(outdir / "panel_proposal.yaml", approved)
    approved_config = replace(
        proposal_config,
        panel=PanelConfig(frozen_manifest=approved),
        primer_design=PrimerDesignConfig(enabled=False),
    )
    second = run_pipeline(
        approved_config,
        outdir,
        execution=ExecutionPolicy(resume=True),
    )

    assert second.status in {"PARTIAL", "COMPLETED"}
    assert (outdir / "panel" / "approved_panel.json").is_file()
    assert (outdir / ".checkpoints" / "panel" / "manifest.json").is_file()
    assert (outdir / ".checkpoints" / "input" / "manifest.json").is_file()
```

The helper's original proposal run may have primer design disabled; an explicitly configured proposal still requires review. Keep a separate Task-3 config test proving enabled primer design cannot omit a panel. This keeps the end-to-end approval test offline and independent of Primer3.

- [ ] **Step 5: Add target-mismatch protection**

Add a test that approves the WNV fixture, then calls `run_pipeline()` with top-level `target_name="Zika virus"` and that frozen manifest. Assert:

```python
with pytest.raises(
    ValueError,
    match="Approved panel target 'West Nile virus'.*pipeline target 'Zika virus'",
):
    run_pipeline(config, tmp_path / "run")
```

The implementation is already in `prepare_panel_preflight()` from Task 4; this test proves the frozen-manifest path, not just inline config validation.

- [ ] **Step 6: Document the exact user workflow**

Add a concise README section named `Panel approval workflow` with these commands and config transitions:

```text
1. Configure `panel.proposal` and run with `--outdir`.
2. Geison returns ACTION_REQUIRED and writes `panel_proposal.yaml`.
3. Review/edit that proposal as a scientific input.
4. Freeze the reviewed proposal:
   qpcr-pipeline panel approve panel_proposal.yaml --output approved_panel.json
5. Replace `panel.proposal` with:
   panel:
     frozen_manifest: approved_panel.json
6. Rerun with `--resume`.
```

State that the frozen manifest is scientific provenance and should be versioned or archived with run inputs. State that automatic panel construction is not part of this subproject.

- [ ] **Step 7: Run the focused panel suite**

```bash
python -m pytest \
  tests/test_panel.py \
  tests/test_panel_manifest.py \
  tests/test_panel_config.py \
  tests/test_pipeline_panel.py \
  tests/test_panel_cli.py \
  tests/test_execution.py \
  tests/test_execution_plan.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py \
  tests/test_dry_run.py \
  tests/test_run_manifest.py -v
```

Expected: PASS.

- [ ] **Step 8: Run the complete regression suite**

```bash
python -m pytest -q
```

Expected: all tests pass. If an existing test fails because it constructs an enabled-primer full pipeline without a panel, make that test's panel state explicit. Do not remove or relax the new validation.

- [ ] **Step 9: Verify CLI behavior from the installed entry point**

```bash
qpcr-pipeline panel approve --help
qpcr-pipeline panel approve tests/fixtures/panels/west_nile_proposal.yaml --output /tmp/geison-west-nile-approved.json
python -c "import json; p=json.load(open('/tmp/geison-west-nile-approved.json')); assert p['status']=='APPROVED' and p['approved_by_user'] is True and p['proposal_sha256'].startswith('sha256:')"
```

Expected: all three commands exit 0.

- [ ] **Step 10: Commit end-to-end coverage and docs**

```bash
git add qpcr_pipeline/dry_run.py qpcr_pipeline/cli.py tests/test_dry_run.py tests/test_pipeline_panel.py README.md
git commit -m "test: cover panel approval workflow"
```

---

## Acceptance Checklist

Before declaring subproject 1 complete, obtain fresh evidence for every item:

- Panel models represent target mode, required target groups, non-target criticality, design/challenge roles, reasons, proposal provenance, sequence-selection provenance placeholders, and diagnostic context.
- Invalid or contradictory panel models fail deterministically.
- A proposal execution returns `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` before input acquisition/checkpointing.
- `qpcr-pipeline run` never asks an interactive approval question.
- `qpcr-pipeline panel approve` produces byte-deterministic canonical JSON for unchanged semantic proposal content.
- Frozen approved manifests reject extra/missing/malformed fields and malformed proposal hashes.
- Approved panels are checkpointed in a first `panel` stage before `input`.
- Changing the approved manifest invalidates `panel` and conservatively forces every downstream stage.
- `run_manifest.json` records only safe panel provenance summary and artifact identity, not sequence content or the full panel definition.
- Schema-1 run manifests created before this feature can resume after optional-field normalization.
- Legacy panel-free non-design workflows still run.
- Enabled primer design cannot run without a configured panel.
- Dry-run reports approval requirements without writing the proposal or any checkpoint.
- The West Nile fixture is offline and explicitly test-only.
- `python -m pytest -q` passes.

## Deferred to Later Subprojects

Do not pull these into subproject 1:

- R/DECIPHER detection or execution.
- DECIPHER-compatible target/non-target design alignment.
- Common `PrimerPairCandidate` engine abstraction.
- Separate hydrolysis-probe design.
- Automatic taxonomy-based panel construction.
- Production arbovirus clinical knowledge base.
- Representative accession acquisition/stratification.
- Assignment of actual accessions to design/challenge partitions.
- Panel-driven target/non-target multi-dataset orchestration.
- CRITICAL/IMPORTANT/BACKGROUND assay scoring behavior.
- New scientific outcomes such as `NO_VALID_ASSAY`.
