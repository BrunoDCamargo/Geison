# Panel Model and Manifests Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the first panel-aware Geison workflow: explicit target/non-target panel metadata, deterministic proposal/approval artifacts, a first-time human approval gate, frozen approved manifests, and checkpoint/provenance integration.

**Architecture:** Introduce a focused panel domain module plus a strict manifest I/O module. A configured panel proposal stops an executing run before scientific design and writes `panel_proposal.yaml`; a separate non-interactive CLI command converts that proposal into a deterministic approved JSON manifest. Approved manifests become a checkpointed `panel` stage placed before `input`, so a manifest identity change conservatively invalidates every downstream stage. Existing non-design workflows may omit `panel`; any run with `primer_design.enabled: true` must configure either a proposal or frozen approved manifest.

**Tech Stack:** Python 3.10+, stdlib dataclasses/hashlib/json/pathlib/tempfile, existing PyYAML dependency, pytest, existing Geison checkpoint and run-recording infrastructure.

**Spec:** `docs/superpowers/specs/2026-09-03-decipher-hybrid-assay-architecture-design.md`

## Global Constraints

- Geison remains the Python orchestrator; this subproject adds no R/DECIPHER dependency.
- The panel scientific model is pathogen-agnostic; the West Nile fixture is test data, not a hard-coded production knowledge base.
- A newly configured panel proposal must not silently proceed into assay design. It must produce `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` and a review artifact.
- Approval must be non-interactive and explicit. Do not prompt inside `qpcr-pipeline run`.
- Approved manifests must be deterministic, strict, immutable artifacts suitable for hashing and reproducible reruns.
- Challenge/design roles are modeled now, but actual representative sequence selection and enforcement against sequence membership remain subproject 4 work.
- Existing `off_targets` configuration remains in place during this subproject; do not migrate specificity datasets yet.
- `primer_design.enabled: true` requires a configured panel. QC/alignment/conservation-only legacy workflows may continue without a panel during this incremental migration.
- Do not introduce hidden biological thresholds or automatic clinical inference in this subproject.
- No live NCBI or network access is required by panel unit/integration tests.

---

## File Structure

### New production files

- `qpcr_pipeline/panel.py`
  - Immutable panel domain models.
  - Validation of target mode, target groups, non-target criticality, design/challenge roles, reasons, provenance labels, and diagnostic context.
  - No file I/O.

- `qpcr_pipeline/panel_manifest.py`
  - Strict YAML proposal serialization/parsing.
  - Deterministic canonical JSON approval serialization/parsing.
  - Semantic proposal SHA-256 calculation.
  - Runtime proposal preflight and approved-manifest materialization.

### Modified production files

- `qpcr_pipeline/config.py`
  - `PanelConfig` and YAML parsing.
  - Add optional `PipelineConfig.panel`.
  - Require a panel when primer design is enabled.

- `qpcr_pipeline/execution.py`
  - Add `panel` as the first stage and make `input` depend on it.

- `qpcr_pipeline/checkpoint_codecs.py`
  - Add a structural codec for `PanelResult`.

- `qpcr_pipeline/checkpoint_stages.py`
  - Add panel stage parameters, approved-manifest file identity, output artifact, and dependency integration.

- `qpcr_pipeline/pipeline.py`
  - Execute panel preflight before normal checkpoint planning.
  - Return an action-required summary for proposals.
  - Execute/materialize approved panel as the first checkpointed stage.

- `qpcr_pipeline/run_recording.py`
  - Add explicit action-required lifecycle state and panel provenance storage.

- `qpcr_pipeline/provenance.py`
  - Project safe panel provenance into the run manifest.

- `qpcr_pipeline/cli.py`
  - Add `qpcr-pipeline panel approve`.
  - Return a distinct exit code for action-required runs.

### New tests/fixtures

- `tests/test_panel.py`
- `tests/test_panel_manifest.py`
- `tests/test_panel_config.py`
- `tests/test_pipeline_panel.py`
- `tests/test_panel_cli.py`
- `tests/fixtures/panels/west_nile_proposal.yaml`

### Modified tests

- `tests/test_execution.py`
- `tests/test_checkpoint_codecs.py`
- `tests/test_checkpoint_stages.py`
- `tests/test_pipeline_resume.py`
- `tests/test_run_manifest.py`
- Any assertion that enumerates the complete `STAGE_ORDER` literally.

---

### Task 1: Add immutable panel domain models and validation

**Files:**
- Create: `qpcr_pipeline/panel.py`
- Create: `tests/test_panel.py`
- Create: `tests/fixtures/panels/west_nile_proposal.yaml`

**Interfaces:**
- Produces:
  - `DatasetRole = Literal["DESIGN", "CHALLENGE"]`
  - `TargetMode = Literal["broad_detection", "subtype_specific"]`
  - `Criticality = Literal["CRITICAL", "IMPORTANT", "BACKGROUND"]`
  - `DiagnosticContext`
  - `SequenceSelectionProvenance`
  - `TargetGroup`
  - `PanelTarget`
  - `PanelNonTarget`
  - `PanelDefinition`
  - `validate_panel_definition(definition: PanelDefinition) -> None`
- Consumes: no new project interfaces.

- [ ] **Step 1: Write the failing domain-model tests**

Create `tests/test_panel.py` with direct-construction tests that define the contract before any parser exists:

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


def _west_nile_panel() -> PanelDefinition:
    selection = SequenceSelectionProvenance(
        dataset_role="DESIGN",
        method="manual_fixture",
        source="tests/fixtures/panels/west_nile_proposal.yaml",
        details=("lineage-aware seed",),
    )
    return PanelDefinition(
        target=PanelTarget(
            name="West Nile virus",
            taxid=11082,
            mode="broad_detection",
            subtype=None,
            groups=(
                TargetGroup(
                    name="lineage_1",
                    required=True,
                    dataset_roles=("DESIGN", "CHALLENGE"),
                    reasons=("target_diversity",),
                    proposed_by=("manual",),
                    sequence_selection=(selection,),
                ),
                TargetGroup(
                    name="lineage_2",
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
                name="Usutu virus",
                taxid=64286,
                criticality="CRITICAL",
                dataset_roles=("DESIGN", "CHALLENGE"),
                reasons=("phylogenetic_neighbor",),
                proposed_by=("manual",),
                sequence_selection=(selection,),
            ),
        ),
        diagnostic_context=DiagnosticContext(
            syndrome="arboviral febrile disease",
            geography="Brazil",
            sample_type="human serum",
            vector="mosquito",
        ),
    )


def test_valid_panel_definition_is_accepted():
    validate_panel_definition(_west_nile_panel())


def test_subtype_specific_requires_subtype():
    panel = _west_nile_panel()
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


def test_non_target_cannot_duplicate_target_taxid():
    panel = _west_nile_panel()
    duplicate = PanelNonTarget(
        name="West Nile alias",
        taxid=11082,
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


def test_dataset_roles_must_be_unique_and_non_empty():
    panel = _west_nile_panel()
    invalid_group = TargetGroup(
        name="lineage_1",
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

Also cover: blank names, non-positive TaxID, duplicate target group names, duplicate non-target names case-insensitively, duplicate non-target TaxIDs, invalid criticality, blank `reasons`/`proposed_by` members, `broad_detection` with a non-null subtype, invalid optional diagnostic-context strings, and invalid `SequenceSelectionProvenance.dataset_role`.

- [ ] **Step 2: Run the tests and verify the module is missing**

Run:

```bash
python -m pytest tests/test_panel.py -v
```

Expected: collection fails because `qpcr_pipeline.panel` does not exist.

- [ ] **Step 3: Implement the panel dataclasses and strict validators**

Create `qpcr_pipeline/panel.py` with these exact public shapes:

```python
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

DatasetRole = Literal["DESIGN", "CHALLENGE"]
TargetMode = Literal["broad_detection", "subtype_specific"]
Criticality = Literal["CRITICAL", "IMPORTANT", "BACKGROUND"]


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


def validate_panel_definition(definition: PanelDefinition) -> None:
    ...
```

Validation rules are deterministic and structural:

- `PanelDefinition`, `PanelTarget`, `TargetGroup`, `PanelNonTarget`, `DiagnosticContext`, and `SequenceSelectionProvenance` must be their expected dataclass types.
- Explicit panel definitions require at least one target group.
- TaxIDs are `None` or positive `int`, excluding `bool`.
- `broad_detection` requires `subtype is None`; `subtype_specific` requires a non-blank subtype.
- `TargetGroup.required` must be a real `bool`.
- `dataset_roles` must be a non-empty tuple, contain only `DESIGN`/`CHALLENGE`, and contain no duplicates.
- `reasons` and `proposed_by` must each be non-empty tuples of non-blank strings.
- `sequence_selection` may be empty in this subproject but, when present, each item must validate.
- Target group names are unique using `.casefold()`.
- Non-target names are unique using `.casefold()`; non-null non-target TaxIDs are unique.
- A non-target may not duplicate the target by non-null TaxID or casefolded name.
- Diagnostic-context fields are `None` or non-blank strings.

- [ ] **Step 4: Add the manual West Nile proposal fixture**

Create `tests/fixtures/panels/west_nile_proposal.yaml` with this normalized scientific seed:

```yaml
schema_version: 1
status: PROPOSED
definition:
  target:
    name: West Nile virus
    taxid: 11082
    mode: broad_detection
    subtype: null
    groups:
      - name: lineage_1
        required: true
        dataset_roles: [DESIGN, CHALLENGE]
        reasons: [target_diversity]
        proposed_by: [manual]
        sequence_selection:
          - dataset_role: DESIGN
            method: manual_fixture
            source: tests/fixtures/panels/west_nile_proposal.yaml
            details: [lineage-aware seed]
      - name: lineage_2
        required: true
        dataset_roles: [DESIGN, CHALLENGE]
        reasons: [target_diversity]
        proposed_by: [manual]
        sequence_selection:
          - dataset_role: DESIGN
            method: manual_fixture
            source: tests/fixtures/panels/west_nile_proposal.yaml
            details: [lineage-aware seed]
  non_targets:
    - name: Usutu virus
      taxid: 64286
      criticality: CRITICAL
      dataset_roles: [DESIGN, CHALLENGE]
      reasons: [phylogenetic_neighbor]
      proposed_by: [manual]
      sequence_selection: []
    - name: Japanese encephalitis virus
      taxid: 11072
      criticality: CRITICAL
      dataset_roles: [DESIGN, CHALLENGE]
      reasons: [phylogenetic_neighbor]
      proposed_by: [manual]
      sequence_selection: []
    - name: Dengue virus
      taxid: 12637
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

Do not treat this fixture as production clinical knowledge.

- [ ] **Step 5: Run the domain tests**

Run:

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
- Consumes: `PanelDefinition`, `validate_panel_definition` from Task 1.
- Produces:
  - `PanelProposal`
  - `ApprovedPanelManifest`
  - `PanelResult`
  - `load_panel_proposal(path: Path) -> PanelProposal`
  - `write_panel_proposal(definition: PanelDefinition, path: Path) -> PanelProposal`
  - `approve_panel_proposal(proposal_path: Path, output_path: Path) -> ApprovedPanelManifest`
  - `load_approved_panel_manifest(path: Path) -> ApprovedPanelManifest`
  - `materialize_approved_panel(manifest_path: Path, outdir: Path) -> PanelResult`
  - `proposal_semantic_sha256(proposal: PanelProposal) -> str`

- [ ] **Step 1: Write failing manifest round-trip and determinism tests**

Create `tests/test_panel_manifest.py` around the West Nile fixture:

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


def test_approval_is_deterministic(tmp_path):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"

    approve_panel_proposal(FIXTURE, first)
    approve_panel_proposal(FIXTURE, second)

    assert first.read_bytes() == second.read_bytes()
    payload = json.loads(first.read_text(encoding="utf-8"))
    assert payload["status"] == "APPROVED"
    assert payload["approved_by_user"] is True
    assert payload["proposal_sha256"].startswith("sha256:")


def test_semantic_hash_ignores_yaml_formatting(tmp_path):
    original = load_panel_proposal(FIXTURE)
    reformatted = tmp_path / "reformatted.yaml"
    reformatted.write_text(
        "# formatting-only change\n" + FIXTURE.read_text(encoding="utf-8"),
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
```

Also test wrong `schema_version`, wrong `status`, `approved_by_user: false`, malformed proposal SHA, invalid nested domain content, duplicate keys represented after YAML load, and atomic output behavior when approval serialization fails.

- [ ] **Step 2: Run the tests and verify failure**

```bash
python -m pytest tests/test_panel_manifest.py -v
```

Expected: import failure for `qpcr_pipeline.panel_manifest`.

- [ ] **Step 3: Implement strict manifest dataclasses and canonical payload conversion**

Create `qpcr_pipeline/panel_manifest.py` with:

```python
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import os
from pathlib import Path
import tempfile
from typing import Literal

import yaml

from qpcr_pipeline.panel import PanelDefinition, validate_panel_definition


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

Use explicit payload-to-dataclass parsing helpers. Do not instantiate nested dataclasses with unchecked `**payload` recursion. Every mapping level must reject unknown and missing fields with clear `ValueError` messages.

Canonical proposal identity is computed from semantic content, not source YAML bytes:

```python
def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        + "\n"
    ).encode("utf-8")


def proposal_semantic_sha256(proposal: PanelProposal) -> str:
    payload = _proposal_payload(proposal)
    digest = hashlib.sha256(_canonical_json_bytes(payload)).hexdigest()
    return "sha256:" + digest
```

Approved JSON must be deterministic pretty JSON (`indent=2`, `sort_keys=True`, one final newline). Approval contains no timestamp so approving an unchanged proposal twice yields identical bytes.

- [ ] **Step 4: Implement atomic writes and runtime materialization**

Use a small `_atomic_write_text(path: Path, text: str) -> None` helper based on `tempfile.mkstemp(..., dir=path.parent)` and `Path.replace()`.

`write_panel_proposal()` writes normalized YAML with `yaml.safe_dump(..., sort_keys=False, allow_unicode=True)`.

`approve_panel_proposal()` performs:

```python
proposal = load_panel_proposal(proposal_path)
approved = ApprovedPanelManifest(
    schema_version=1,
    status="APPROVED",
    approved_by_user=True,
    proposal_sha256=proposal_semantic_sha256(proposal),
    definition=proposal.definition,
)
write_approved_panel_manifest(approved, output_path)
return approved
```

`materialize_approved_panel()` validates the source manifest, writes a canonical copy to `outdir / "panel" / "approved_panel.json"`, calculates the copied file SHA-256 as `sha256:<hex>`, and returns `PanelResult(status="APPROVED", ...)`.

- [ ] **Step 5: Run manifest tests**

```bash
python -m pytest tests/test_panel.py tests/test_panel_manifest.py -v
```

Expected: PASS.

- [ ] **Step 6: Commit manifest I/O**

```bash
git add qpcr_pipeline/panel_manifest.py tests/test_panel_manifest.py
git commit -m "feat: add deterministic panel manifests"
```

---

### Task 3: Add panel configuration and require it for primer design

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Create: `tests/test_panel_config.py`
- Modify: `tests/test_config.py` only where existing primer-design fixtures need an explicit panel.

**Interfaces:**
- Consumes: `PanelDefinition` and parser helpers implemented in Tasks 1-2.
- Produces:
  - `PanelConfig`
  - `PipelineConfig.panel: PanelConfig | None`
  - YAML support for `panel.proposal` and `panel.frozen_manifest`.

- [ ] **Step 1: Write failing config parsing tests**

Create `tests/test_panel_config.py` with helpers that write a minimal input FASTA and config YAML. Cover both configuration modes:

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
      taxid: 11082
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
        taxid: 64286
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
```

Add tests for:

- `panel.frozen_manifest` path parsing.
- Exactly one of `proposal` or `frozen_manifest` must be configured when `panel` exists.
- Unknown fields rejected at `panel`, `target`, target-group, non-target, context, and sequence-selection levels.
- Top-level `target.name` must equal inline proposal target name exactly after `.strip()` normalization.
- `primer_design.enabled: true` with no panel raises `Enabled primer design requires a panel proposal or frozen manifest.`
- Primer design disabled with no panel remains valid.

- [ ] **Step 2: Run the tests and verify failure**

```bash
python -m pytest tests/test_panel_config.py -v
```

Expected: FAIL because `PipelineConfig` has no panel field/parser.

- [ ] **Step 3: Add `PanelConfig` and parser integration**

In `qpcr_pipeline/config.py`, import the panel domain models and define:

```python
@dataclass(frozen=True, slots=True)
class PanelConfig:
    proposal: PanelDefinition | None = None
    frozen_manifest: Path | None = None
```

Add:

```python
@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path | None = None
    input_genbank: Path | None = None
    input_ncbi: NcbiInputConfig | None = None
    panel: PanelConfig | None = None
    ...
```

`load_config()` reads `panel = _parse_panel_config(raw.get("panel"))`; absence returns `None`.

Implement strict helpers in `config.py`:

```python
def _parse_panel_config(raw: Any) -> PanelConfig | None: ...
def _parse_panel_definition(raw: Any) -> PanelDefinition: ...
def _parse_panel_target(raw: Any) -> PanelTarget: ...
def _parse_target_group(raw: Any, *, index: int) -> TargetGroup: ...
def _parse_panel_non_target(raw: Any, *, index: int) -> PanelNonTarget: ...
def _parse_diagnostic_context(raw: Any) -> DiagnosticContext: ...
def _parse_sequence_selection(raw: Any, *, section: str) -> tuple[SequenceSelectionProvenance, ...]: ...
def validate_panel_config(config: PanelConfig, *, target_name: str) -> None: ...
```

Use the same strict-unknown-field style already used by clustering/off-target parsing.

- [ ] **Step 4: Enforce the incremental primer-design gate**

In `validate_pipeline_config()` add:

```python
if config.panel is not None:
    if not isinstance(config.panel, PanelConfig):
        raise ValueError("Pipeline panel must be a PanelConfig when configured.")
    validate_panel_config(config.panel, target_name=config.target_name)

if config.primer_design.enabled and config.panel is None:
    raise ValueError(
        "Enabled primer design requires a panel proposal or frozen manifest."
    )
```

Do not require a panel for target-only QC/alignment/conservation workflows yet.

- [ ] **Step 5: Repair existing primer-design config fixtures intentionally**

Search:

```bash
grep -R "primer_design" -n tests | head -50
```

For tests whose purpose is primer design rather than panel validation, construct a temporary approved panel helper in later pipeline tests or use a minimal direct `PanelConfig(proposal=...)` where no execution occurs. Do not weaken the new validation merely to preserve old tests.

- [ ] **Step 6: Run config suites**

```bash
python -m pytest tests/test_panel_config.py tests/test_config.py tests/test_specificity_config.py tests/test_ranking_config.py -v
```

Expected: PASS.

- [ ] **Step 7: Commit config integration**

```bash
git add qpcr_pipeline/config.py tests/test_panel_config.py tests/test_config.py
git commit -m "feat: configure assay panels"
```

---

### Task 4: Add non-interactive action-required preflight and approval lifecycle

**Files:**
- Modify: `qpcr_pipeline/panel_manifest.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `qpcr_pipeline/run_recording.py`
- Create: `tests/test_pipeline_panel.py`
- Modify: `tests/test_run_manifest.py`

**Interfaces:**
- Consumes: `PanelConfig` from Task 3; manifest I/O from Task 2.
- Produces:
  - `PanelPreflight`
  - `prepare_panel_preflight(panel: PanelConfig | None, outdir: Path) -> PanelPreflight`
  - `RunRecorder.action_required(code: str, artifact: Path) -> None`
  - `RunSummary.action_required_code`
  - `RunSummary.action_required_artifact`

- [ ] **Step 1: Write the first-run action-required test**

In `tests/test_pipeline_panel.py`, construct a config with `primer_design.enabled=True`, conservation/alignment enabled, a tiny local FASTA, and `PanelConfig(proposal=west_nile_definition)`.

```python
def test_panel_proposal_stops_before_input_and_writes_review_artifact(tmp_path):
    config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"

    summary = run_pipeline(config, outdir)

    assert summary.status == "ACTION_REQUIRED"
    assert summary.action_required_code == "PANEL_APPROVAL_REQUIRED"
    assert summary.sequence_count == 0
    assert (outdir / "panel_proposal.yaml").is_file()
    assert not (outdir / ".checkpoints" / "input" / "manifest.json").exists()
    assert not (outdir / "qc_report.json").exists()
```

Add a run-manifest assertion:

```python
manifest = json.loads((outdir / "run_manifest.json").read_text(encoding="utf-8"))
assert manifest["status"] == "ACTION_REQUIRED"
assert manifest["action_required"] == {
    "code": "PANEL_APPROVAL_REQUIRED",
    "artifact": str(outdir / "panel_proposal.yaml"),
}
```

- [ ] **Step 2: Run the test and verify failure**

```bash
python -m pytest tests/test_pipeline_panel.py::test_panel_proposal_stops_before_input_and_writes_review_artifact -v
```

Expected: FAIL because action-required lifecycle is not implemented.

- [ ] **Step 3: Add `PanelPreflight` and deterministic proposal output**

In `panel_manifest.py` add:

```python
@dataclass(frozen=True, slots=True)
class PanelPreflight:
    status: Literal["READY", "ACTION_REQUIRED", "LEGACY"]
    proposal_path: Path | None = None
    approved_manifest_path: Path | None = None


def prepare_panel_preflight(panel: PanelConfig | None, outdir: Path) -> PanelPreflight:
    if panel is None:
        return PanelPreflight(status="LEGACY")
    if panel.proposal is not None:
        proposal_path = Path(outdir) / "panel_proposal.yaml"
        write_panel_proposal(panel.proposal, proposal_path)
        return PanelPreflight(
            status="ACTION_REQUIRED",
            proposal_path=proposal_path,
        )
    assert panel.frozen_manifest is not None
    load_approved_panel_manifest(panel.frozen_manifest)
    return PanelPreflight(
        status="READY",
        approved_manifest_path=panel.frozen_manifest,
    )
```

Import `PanelConfig` lazily or move the type under `TYPE_CHECKING` to avoid a runtime import cycle (`config.py` imports panel domain models; `panel_manifest.py` may import `PanelConfig`).

- [ ] **Step 4: Extend run recording with `ACTION_REQUIRED`**

In `run_recording.py`:

- Add `"run_action_required": frozenset({"status", "code", "artifact"})` to `EVENT_FIELDS`.
- Initialize top-level `"action_required": None` and `"panel_provenance": {}` in a new run manifest.
- Clear stale `action_required` when beginning a new attempt.
- Add:

```python
def action_required(self, code: str, artifact: Path) -> None:
    if not code:
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

- [ ] **Step 5: Add panel preflight to `run_pipeline()` before checkpoint planning**

Extend `RunSummary`:

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

After output-dir creation and environment inspection, but before `plan_pipeline()`, call `prepare_panel_preflight()`.

For `ACTION_REQUIRED`:

1. Build a minimal recorded plan containing one row: `{"stage": "panel", "action": "ACTION_REQUIRED", "reason": "panel approval required before assay design"}`.
2. `recorder.begin_attempt(...)`.
3. Remove stale `run_summary.json` and `qc_report.json`.
4. `recorder.action_required("PANEL_APPROVAL_REQUIRED", proposal_path)`.
5. Write `run_summary.json` atomically with the action-required summary.
6. Return without creating any scientific checkpoints.

Do not raise an exception for normal action-required state.

- [ ] **Step 6: Test repeatability and stale-state cleanup**

Add tests that:

- Running the same proposal twice rewrites identical `panel_proposal.yaml` content.
- An earlier `run_summary.json` is replaced by an `ACTION_REQUIRED` summary.
- An action-required attempt is preserved when a later approved run reuses the same output directory.

- [ ] **Step 7: Run panel lifecycle and run-manifest tests**

```bash
python -m pytest tests/test_pipeline_panel.py tests/test_run_manifest.py tests/test_run_recording.py tests/test_run_logging.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit action-required lifecycle**

```bash
git add qpcr_pipeline/panel_manifest.py qpcr_pipeline/pipeline.py qpcr_pipeline/run_recording.py tests/test_pipeline_panel.py tests/test_run_manifest.py
git commit -m "feat: gate assay design on panel approval"
```

---

### Task 5: Add approved panel as a checkpointed first stage

**Files:**
- Modify: `qpcr_pipeline/execution.py`
- Modify: `qpcr_pipeline/checkpoint_codecs.py`
- Modify: `qpcr_pipeline/checkpoint_stages.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_execution.py`
- Modify: `tests/test_checkpoint_codecs.py`
- Modify: `tests/test_checkpoint_stages.py`
- Modify: `tests/test_pipeline_resume.py`

**Interfaces:**
- Consumes: `PanelResult`, `materialize_approved_panel()`.
- Produces: checkpointed stage named `panel`, always first in `STAGE_ORDER`; `input` depends on `panel`.

- [ ] **Step 1: Write failing execution-graph tests**

Update/add in `tests/test_execution.py`:

```python
def test_panel_is_first_stage_and_input_depends_on_it():
    assert STAGE_ORDER[0] == "panel"
    assert STAGE_DEPENDENCIES["panel"] == ()
    assert STAGE_DEPENDENCIES["input"] == ("panel",)


def test_invalid_panel_forces_entire_pipeline_on_resume():
    reusable = _all_valid()
    reusable["panel"] = False
    actions = _actions(plan_from_validity(ExecutionPolicy(resume=True), reusable))
    assert actions["panel"] == "RUN"
    assert all(actions[stage] == "FORCED" for stage in STAGE_ORDER[1:])
```

Import `STAGE_DEPENDENCIES` in this test.

- [ ] **Step 2: Run execution tests and verify failure**

```bash
python -m pytest tests/test_execution.py -v
```

Expected: FAIL because `panel` is not a stage.

- [ ] **Step 3: Add the stage to `execution.py`**

Change the stage declarations to:

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

STAGE_ORDER = (
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

STAGE_DEPENDENCIES = {
    "panel": (),
    "input": ("panel",),
    ...
}
```

No other dependency needs changing in this subproject because transitive invalidation now flows through `input`.

- [ ] **Step 4: Add `PANEL_CODEC`**

In `checkpoint_codecs.py` import `PanelResult` and add:

```python
PANEL_CODEC = _StructuralCodec[PanelResult](PanelResult)
```

Add a round-trip test using both:

```python
PanelResult(
    status="LEGACY",
    manifest_sha256=None,
    manifest_path=None,
    target_mode=None,
    non_target_count=0,
)
```

and an approved result whose `manifest_path` is inside the test output directory.

- [ ] **Step 5: Add checkpoint request identity for the panel stage**

In `checkpoint_stages.py`:

- Import `PANEL_CODEC`, `PanelResult`, and `load_approved_panel_manifest` only where needed.
- Add `"panel": PANEL_CODEC` to `_CODECS`.
- `stage_parameters("panel", config)` returns:

```python
if config.panel is None:
    return {"mode": "LEGACY"}
return {"mode": "APPROVED_MANIFEST"}
```

A proposal must never reach checkpoint planning because preflight returns action-required first.

- `stage_input_identities("panel", config)` returns `{}` for legacy mode; for frozen manifest:

```python
return {
    "approved_panel": {
        "sha256": file_sha256(config.panel.frozen_manifest),
    }
}
```

Validate the manifest before returning the identity so corrupted approved manifests fail before execution.

- `stage_outputs("panel", result, outdir)` returns `result.manifest_path` when present.

- [ ] **Step 6: Execute the stage in `pipeline.py`**

In `_run_stage()` handle panel before reading `results["input"]`:

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
        raise RuntimeError("Panel proposal reached checkpoint execution before approval.")
    return materialize_approved_panel(config.panel.frozen_manifest, output_dir)
```

The existing `input` branch follows this and therefore only runs after a completed panel checkpoint.

- [ ] **Step 7: Update resume tests for conservative panel invalidation**

Add an integration test in `tests/test_pipeline_resume.py`:

1. Approve a West Nile proposal into `approved-a.json`.
2. Run a minimal pipeline with `PanelConfig(frozen_manifest=approved-a.json)`.
3. Resume unchanged and assert all checkpoints, including `panel`, are `REUSE`.
4. Produce `approved-b.json` from a scientifically changed proposal (for example change Usutu `criticality` from `CRITICAL` to `IMPORTANT`).
5. Resume with `approved-b.json` and assert `panel` is `RUN` and every downstream stage is `FORCED`.

This subproject deliberately uses conservative full invalidation. More selective dependency boundaries can be introduced only after the panel starts controlling separate target/non-target sequence acquisition.

- [ ] **Step 8: Run checkpoint/resume suites**

```bash
python -m pytest \
  tests/test_execution.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_execution_plan.py \
  tests/test_pipeline_resume.py -v
```

Expected: PASS.

- [ ] **Step 9: Commit checkpoint integration**

```bash
git add qpcr_pipeline/execution.py qpcr_pipeline/checkpoint_codecs.py qpcr_pipeline/checkpoint_stages.py qpcr_pipeline/pipeline.py tests/test_execution.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py tests/test_pipeline_resume.py
git commit -m "feat: checkpoint approved assay panels"
```

---

### Task 6: Add panel provenance and explicit CLI approval command

**Files:**
- Modify: `qpcr_pipeline/provenance.py`
- Modify: `qpcr_pipeline/run_recording.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `qpcr_pipeline/cli.py`
- Create: `tests/test_panel_cli.py`
- Modify: `tests/test_run_manifest.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Consumes: panel manifest and `PanelResult` from earlier tasks.
- Produces:
  - `build_panel_provenance(panel_result: object) -> dict[str, object]`
  - `qpcr-pipeline panel approve <proposal> --output <manifest>`
  - CLI exit code `3` for normal `ACTION_REQUIRED` run results.

- [ ] **Step 1: Write failing CLI parser and command tests**

Create `tests/test_panel_cli.py` using `monkeypatch` on `sys.argv` and `capsys`, matching the existing CLI test style:

```python
def test_panel_approve_command_writes_frozen_manifest(tmp_path, monkeypatch, capsys):
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
    assert "Approved panel manifest:" in capsys.readouterr().out
```

Add a test that `qpcr-pipeline run ... --outdir ...` returns `3` when the pipeline summary status is `ACTION_REQUIRED` and prints the proposal path.

- [ ] **Step 2: Run CLI tests and verify failure**

```bash
python -m pytest tests/test_panel_cli.py -v
```

Expected: FAIL because the nested `panel approve` command does not exist.

- [ ] **Step 3: Add nested CLI parser**

In `build_parser()`:

```python
panel_parser = subparsers.add_parser("panel", help="Review and freeze assay panels")
panel_subparsers = panel_parser.add_subparsers(dest="panel_command", required=True)
approve_parser = panel_subparsers.add_parser(
    "approve",
    help="Approve a reviewed panel proposal and write a frozen manifest",
)
approve_parser.add_argument("proposal", type=Path)
approve_parser.add_argument("--output", type=Path, required=True)
```

In `main()`:

```python
if args.command == "panel" and args.panel_command == "approve":
    manifest = approve_panel_proposal(args.proposal, args.output)
    print(f"Approved panel manifest: {args.output}")
    print(f"Proposal identity: {manifest.proposal_sha256}")
    return 0
```

For `run`:

```python
if summary.status == "ACTION_REQUIRED":
    print(
        f"ACTION_REQUIRED: {summary.action_required_code} "
        f"({summary.action_required_artifact})"
    )
    return 3
```

Keep technical exceptions as exceptions/current failure behavior; do not map them to action-required.

- [ ] **Step 4: Add safe panel provenance**

In `provenance.py`:

```python
def build_panel_provenance(panel_result: object) -> dict[str, object]:
    status = getattr(panel_result, "status", None)
    if status == "LEGACY":
        return {"mode": "legacy_unconfigured"}
    if status != "APPROVED":
        raise ValueError("Panel provenance requires an approved or legacy panel result.")
    manifest_sha256 = getattr(panel_result, "manifest_sha256", None)
    if _sha256_identity(manifest_sha256) is None:
        raise ValueError("Approved panel provenance is missing a SHA-256 identity.")
    return {
        "mode": "approved_manifest",
        "manifest_sha256": _sha256_identity(manifest_sha256),
        "target_mode": getattr(panel_result, "target_mode", None),
        "non_target_count": getattr(panel_result, "non_target_count", 0),
    }
```

Do not copy the full panel or sequence content into `run_manifest.json`; the frozen manifest artifact is the detailed record.

- [ ] **Step 5: Persist panel provenance on successful runs**

Extend `RunRecorder.complete()` with a `panel_provenance: object` argument and store sanitized output in top-level `panel_provenance`.

Update the only pipeline call site:

```python
recorder.complete(
    final_status,
    final_completeness,
    build_input_provenance(...),
    build_reference_provenance(alignment),
    build_panel_provenance(results["panel"]),
)
```

Keep argument order explicit/keyworded if that improves readability; update direct unit-test calls accordingly.

- [ ] **Step 6: Add run-manifest provenance tests**

In `tests/test_run_manifest.py` assert for an approved-panel run:

```python
assert manifest["panel_provenance"] == {
    "mode": "approved_manifest",
    "manifest_sha256": pytest.helpers.any_sha256 if such a helper exists else manifest["panel_provenance"]["manifest_sha256"],
    "target_mode": "broad_detection",
    "non_target_count": 3,
}
assert manifest["panel_provenance"]["manifest_sha256"].startswith("sha256:")
```

Do not introduce a new global pytest helper solely for this assertion; if none exists, assert fields individually.

For a non-design legacy run, assert:

```python
assert manifest["panel_provenance"] == {"mode": "legacy_unconfigured"}
```

- [ ] **Step 7: Run CLI/provenance suites**

```bash
python -m pytest tests/test_panel_cli.py tests/test_cli.py tests/test_run_manifest.py tests/test_run_recording.py -v
```

Expected: PASS.

- [ ] **Step 8: Commit CLI and provenance**

```bash
git add qpcr_pipeline/provenance.py qpcr_pipeline/run_recording.py qpcr_pipeline/pipeline.py qpcr_pipeline/cli.py tests/test_panel_cli.py tests/test_cli.py tests/test_run_manifest.py tests/test_run_recording.py
git commit -m "feat: approve and record assay panels"
```

---

### Task 7: Add end-to-end panel workflow coverage and regression verification

**Files:**
- Modify: `tests/test_pipeline_panel.py`
- Modify: `tests/test_minimal_run.py` only where primer-design-enabled cases require approved panel setup.
- Modify: `tests/test_dry_run.py` for the added `panel` stage.
- Modify: `README.md` with the minimal proposal → approve → rerun workflow.

**Interfaces:**
- Consumes all earlier tasks.
- Produces a working, offline, reviewable subproject-1 workflow.

- [ ] **Step 1: Add an end-to-end proposal → approval → rerun test**

In `tests/test_pipeline_panel.py` implement one complete offline scenario:

```python
def test_proposal_approve_rerun_reaches_scientific_pipeline(tmp_path):
    proposal_config = panel_proposal_pipeline_config(tmp_path)
    outdir = tmp_path / "run"

    first = run_pipeline(proposal_config, outdir)
    assert first.status == "ACTION_REQUIRED"

    approved = tmp_path / "approved_panel.json"
    approve_panel_proposal(outdir / "panel_proposal.yaml", approved)

    approved_config = replace(
        proposal_config,
        panel=PanelConfig(frozen_manifest=approved),
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

The exact final scientific status is intentionally not the assertion target here; this subproject validates the panel gate and checkpoint transition.

- [ ] **Step 2: Add a target mismatch failure test for frozen manifests**

Approve the WNV fixture, then run a config with `target_name="Zika virus"` and that manifest. Runtime preflight or panel materialization must raise a clear error such as:

```text
Approved panel target 'West Nile virus' does not match pipeline target 'Zika virus'.
```

Implement this target consistency check in `prepare_panel_preflight()` after loading the approved manifest. Pass the pipeline target name into preflight:

```python
def prepare_panel_preflight(
    panel: PanelConfig | None,
    outdir: Path,
    *,
    target_name: str,
) -> PanelPreflight:
    ...
```

Use the same check for inline proposals in `validate_panel_config()`.

- [ ] **Step 3: Update dry-run expectations**

Dry run with a frozen approved manifest should include `panel` as the first stage. Dry run with an inline proposal must not mutate the filesystem; it should still report that execution cannot proceed without approval.

If the current `dry_run_pipeline()` cannot represent action-required, extend its report with:

```python
panel_action_required: bool = False
panel_proposal_would_be_written: str | None = None
```

Do **not** write `panel_proposal.yaml` during `--dry-run`.

- [ ] **Step 4: Document the exact user workflow**

Add a concise `README.md` section:

```text
1. Configure `panel.proposal` and run with `--outdir`.
2. Geison exits with ACTION_REQUIRED and writes `panel_proposal.yaml`.
3. Review/edit the proposal outside the running pipeline.
4. Freeze it explicitly:
   qpcr-pipeline panel approve panel_proposal.yaml --output approved_panel.json
5. Replace `panel.proposal` with `panel.frozen_manifest: approved_panel.json`.
6. Rerun with `--resume`.
```

State that the approved manifest is scientific provenance and should be versioned with the project/run inputs.

- [ ] **Step 5: Run the focused panel suite**

```bash
python -m pytest \
  tests/test_panel.py \
  tests/test_panel_manifest.py \
  tests/test_panel_config.py \
  tests/test_pipeline_panel.py \
  tests/test_panel_cli.py \
  tests/test_execution.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py \
  tests/test_dry_run.py \
  tests/test_run_manifest.py -v
```

Expected: PASS.

- [ ] **Step 6: Run the complete test suite**

```bash
python -m pytest -q
```

Expected: all tests pass. Any failure caused by the new primer-design panel requirement must be fixed by making the relevant test's panel state explicit, not by weakening validation.

- [ ] **Step 7: Verify CLI help and a real offline approval command**

Run:

```bash
qpcr-pipeline panel approve --help
qpcr-pipeline panel approve tests/fixtures/panels/west_nile_proposal.yaml --output /tmp/geison-west-nile-approved.json
python -c "import json; p=json.load(open('/tmp/geison-west-nile-approved.json')); assert p['status']=='APPROVED' and p['approved_by_user'] is True"
```

Expected: help exits 0; approval exits 0; JSON assertion exits 0.

- [ ] **Step 8: Commit end-to-end coverage and docs**

```bash
git add tests/test_pipeline_panel.py tests/test_minimal_run.py tests/test_dry_run.py README.md
git commit -m "test: cover panel approval workflow"
```

---

## Acceptance Checklist

Before declaring subproject 1 complete, verify all of the following with fresh evidence:

- A valid explicit panel can represent target mode, required target groups, non-target criticality, design/challenge roles, reasons, proposal provenance, sequence-selection provenance placeholders, and diagnostic context.
- Invalid/contradictory panel models fail deterministically.
- A proposal run writes `panel_proposal.yaml` and returns `ACTION_REQUIRED / PANEL_APPROVAL_REQUIRED` before any scientific input checkpoint is created.
- `qpcr-pipeline panel approve` is the only approval action added; `run` never asks interactively.
- Approval of identical semantic proposal content is byte-for-byte deterministic.
- Frozen approved manifests are strict and reject extra/missing/malformed fields.
- An approved panel becomes a first-class `panel` checkpoint before `input`.
- Changing the approved manifest invalidates `panel` and conservatively forces all descendants.
- `run_manifest.json` records only safe panel provenance summary plus the frozen panel artifact identity, not sequence content.
- Legacy non-design workflows can omit a panel, but primer design cannot.
- Dry-run remains non-mutating.
- The West Nile fixture is offline and explicitly test-only.
- Full `python -m pytest -q` passes.

## Deferred to Later Subprojects

The following are intentionally excluded from this plan and must not be pulled into implementation opportunistically:

- DECIPHER/R environment detection or adapter code.
- DECIPHER-compatible homologous non-target alignment.
- `PrimerPairCandidate` engine abstraction.
- Separate hydrolysis-probe design.
- Automatic taxonomy/clinical panel construction.
- Arbovirus production knowledge base.
- Representative sequence acquisition/stratification.
- Actual assignment of sequence accessions into design/challenge partitions.
- Full target/non-target dataset orchestration from the panel.
- New CRITICAL/IMPORTANT/BACKGROUND specificity scoring behavior.
- New scientific outcome model such as `NO_VALID_ASSAY`.
