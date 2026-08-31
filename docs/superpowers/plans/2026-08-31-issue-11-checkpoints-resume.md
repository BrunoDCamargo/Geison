# Issue #11 Checkpoints, Resume, and Selective Invalidation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make every Geison pipeline stage checkpointable so interrupted runs can resume safely and configuration changes recompute only the affected dependency subgraph.

**Architecture:** Add a generic `CheckpointManager` for atomic manifests, typed state snapshots, hashes, and validation; an explicit execution DAG/policy planner; stage-specific codecs and fingerprint projections; then route the existing `run_pipeline()` orchestration through those components. Scientific stage functions remain responsible for science and artifact publication, while checkpoint logic stays in focused infrastructure modules.

**Tech Stack:** Python >=3.10, stdlib `dataclasses`, `hashlib`, `json`, `importlib.metadata`, `subprocess`, `pathlib`; existing Biopython/PyYAML application dependencies; pytest/unittest test suite; existing CircleCI Python 3.12 job.

**Spec:** `docs/superpowers/specs/2026-08-31-issue-11-checkpoints-resume-design.md`

## Global Constraints

- Checkpoint schema version is `1`.
- Checkpoint statuses are exactly `RUNNING`, `COMPLETE`, and `FAILED`; only `COMPLETE` is reusable.
- Checkpoints live only under `<outdir>/.checkpoints`; no global/shared cache.
- No `pickle`, workflow-engine dependency, or cross-project cache discovery.
- Normal execution without resume flags always recomputes every defined stage and refreshes completed online NCBI acquisition.
- `--resume` may reuse only checkpoints whose causal fingerprint, result fingerprint, state hash, output hashes, schema, stage identity, and typed state all validate.
- `--from-step` is strict: all required boundary checkpoints are validated before any scientific write.
- `--force-step` requires `--resume` and forces the selected stage plus transitive dependents.
- Relevant configuration is projected per stage; never hash the entire `PipelineConfig` for every stage.
- Changing Geison version invalidates every checkpoint. External tool identity invalidates only the stage that actually uses that tool and its dependents.
- Existing frozen-NCBI output-directory guard remains in force before writes.
- Existing scientific artifact formats remain user-facing and authoritative.
- Conservation never treats the shared root `report.html` alias as checkpoint-critical; ranking owns/hashes it only when ranking actually publishes it.
- Scientifically disabled stages may return `SKIPPED` while their checkpoint status is `COMPLETE`.
- Existing `run_summary.json` fields remain compatible and gain deterministic stage actions.
- Existing CircleCI policy remains: normal unit tests on `develop`; real external-tool integration tests on `main`.

---

## File Structure

**Create**

- `qpcr_pipeline/checkpointing.py`: canonical JSON hashing, manifest model, state codec protocol, atomic checkpoint lifecycle, validation, and typed checkpoint loading.
- `qpcr_pipeline/execution.py`: stage names/order/dependencies, execution policy, graph traversal, strict `--from-step` boundary calculation, and RUN/REUSE/FORCED planning.
- `qpcr_pipeline/checkpoint_codecs.py`: strict JSON codecs for input records and every current stage result.
- `qpcr_pipeline/checkpoint_stages.py`: stage-specific relevant-config projections, external-input identities, durable output ownership, Geison/tool identities, and checkpoint request construction.
- `tests/test_checkpointing.py`: manifest/fingerprint/hash/lifecycle tests.
- `tests/test_execution.py`: DAG and policy tests.
- `tests/test_checkpoint_codecs.py`: round-trip and malformed-state rejection tests.
- `tests/test_checkpoint_stages.py`: selective projection, external input, output ownership, and tool-identity tests.
- `tests/test_pipeline_resume.py`: end-to-end interruption/resume and selective invalidation tests with fake scientific runners.

**Modify**

- `qpcr_pipeline/pipeline.py`: execute/reuse stages through the checkpoint plan; preserve scientific summaries and add `stage_actions`.
- `qpcr_pipeline/cli.py`: expose and validate `--resume`, `--from-step`, and `--force-step`.
- `qpcr_pipeline/primer_design.py`: expose a pure `primer3_required()` decision using the same candidate-selection logic as `design_primers`, so tool identity is included only when Primer3 can actually run.
- `tests/test_cli.py`: parser/CLI contract tests.
- `tests/test_minimal_run.py`: regression assertions for normal full recomputation and unchanged scientific outputs.
- `README.md`: user-facing checkpoint/resume examples and semantics.

---

### Task 1: Build the generic checkpoint persistence layer

**Files:**
- Create: `qpcr_pipeline/checkpointing.py`
- Create: `tests/test_checkpointing.py`

**Interfaces:**
- Consumes: a stage name, canonical checkpoint request metadata, typed state codec, output paths, and `<outdir>`.
- Produces:
  - `CHECKPOINT_SCHEMA_VERSION = 1`
  - `CheckpointStatus = Literal["RUNNING", "COMPLETE", "FAILED"]`
  - `CheckpointInvalidity` enum/string values for `MISSING_MANIFEST`, `INVALID_MANIFEST`, `NON_COMPLETE_STATUS`, `FINGERPRINT_MISMATCH`, `MISSING_STATE`, `STATE_HASH_MISMATCH`, `INVALID_STATE`, `MISSING_OUTPUT`, `OUTPUT_HASH_MISMATCH`, `RESULT_FINGERPRINT_MISMATCH`.
  - `StateCodec[T]` protocol with `encode(value: T, outdir: Path) -> object` and `decode(payload: object, outdir: Path) -> T`.
  - `CheckpointRequest(stage, dependencies, inputs, parameters, software, tools)`.
  - `OutputIdentity(path: str, sha256: str)`.
  - `CheckpointManifest(...)` with causal `fingerprint` and `result_fingerprint`.
  - `CheckpointLoad[T](state, manifest)` and `CheckpointValidation[T](valid, invalidity, detail, loaded)`.
  - `canonical_sha256(value: object) -> str`.
  - `file_sha256(path: Path) -> str`.
  - `causal_fingerprint(request: CheckpointRequest) -> str`.
  - `result_fingerprint(causal: str, state_sha256: str, outputs: tuple[OutputIdentity, ...]) -> str`.
  - `CheckpointManager.begin(request)`, `complete(request, state, codec, outputs)`, `fail(request, error)`, and `validate(request, codec)`.

- [ ] **Step 1: Write failing canonical-hash and fingerprint tests**

Create tests proving key ordering does not change the digest, dependency result fingerprints do change the causal digest, and result bytes change only `result_fingerprint`:

```python
from pathlib import Path

from qpcr_pipeline.checkpointing import (
    CheckpointRequest,
    OutputIdentity,
    canonical_sha256,
    causal_fingerprint,
    result_fingerprint,
)


def test_canonical_sha256_is_mapping_order_independent():
    assert canonical_sha256({"b": 2, "a": 1}) == canonical_sha256({"a": 1, "b": 2})


def test_causal_fingerprint_uses_dependency_result_identity():
    base = dict(
        stage="alignment",
        inputs={},
        parameters={"enabled": True},
        software={"geison": "0.1.0.dev0"},
        tools={"mafft": {"version": "7.526"}},
    )
    first = CheckpointRequest(dependencies={"clustering": "sha256:a"}, **base)
    second = CheckpointRequest(dependencies={"clustering": "sha256:b"}, **base)
    assert causal_fingerprint(first) != causal_fingerprint(second)


def test_result_fingerprint_uses_state_and_sorted_outputs():
    first = result_fingerprint(
        "sha256:causal",
        "sha256:state-a",
        (OutputIdentity("a.tsv", "sha256:a"), OutputIdentity("b.json", "sha256:b")),
    )
    reordered = result_fingerprint(
        "sha256:causal",
        "sha256:state-a",
        (OutputIdentity("b.json", "sha256:b"), OutputIdentity("a.tsv", "sha256:a")),
    )
    changed = result_fingerprint(
        "sha256:causal",
        "sha256:state-b",
        (OutputIdentity("a.tsv", "sha256:a"), OutputIdentity("b.json", "sha256:b")),
    )
    assert first == reordered
    assert first != changed
```

- [ ] **Step 2: Run the focused tests and verify RED**

Run: `python -m pytest -q tests/test_checkpointing.py`

Expected: import failure because `qpcr_pipeline.checkpointing` does not exist.

- [ ] **Step 3: Implement canonical hashing and manifest dataclasses**

Use deterministic JSON with `sort_keys=True`, compact separators, UTF-8, no timestamps in fingerprints, and `sha256:` prefixes:

```python
def _canonical_bytes(value: object) -> bytes:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")


def canonical_sha256(value: object) -> str:
    return "sha256:" + hashlib.sha256(_canonical_bytes(value)).hexdigest()
```

`CheckpointRequest.fingerprint_payload()` must contain exactly schema version, stage, dependencies, inputs, parameters, software, and tools. `result_fingerprint()` must sort outputs by relative path before hashing.

- [ ] **Step 4: Add failing lifecycle/integrity tests**

Use a tiny test codec:

```python
class DictCodec:
    def encode(self, value, outdir):
        del outdir
        return value

    def decode(self, payload, outdir):
        del outdir
        if not isinstance(payload, dict) or "value" not in payload:
            raise ValueError("invalid test state")
        return payload
```

Cover:

```python
def test_complete_checkpoint_round_trips_and_validates(tmp_path): ...
def test_running_checkpoint_is_not_reusable(tmp_path): ...
def test_failed_checkpoint_is_not_reusable(tmp_path): ...
def test_missing_output_is_invalid(tmp_path): ...
def test_modified_output_hash_is_invalid(tmp_path): ...
def test_modified_state_hash_is_invalid(tmp_path): ...
def test_malformed_state_is_invalid(tmp_path): ...
def test_modified_result_fingerprint_is_invalid(tmp_path): ...
def test_output_outside_outdir_is_rejected(tmp_path): ...
```

- [ ] **Step 5: Run the lifecycle tests and verify RED**

Run: `python -m pytest -q tests/test_checkpointing.py`

Expected: failures for missing `CheckpointManager` behavior.

- [ ] **Step 6: Implement atomic lifecycle and validation**

Use `<outdir>/.checkpoints/<stage>/manifest.json` and `state.json`. `begin()` atomically writes `RUNNING`. `complete()` atomically writes state, hashes declared outputs, computes result fingerprint, then atomically replaces the manifest with `COMPLETE`. `fail()` best-effort writes `FAILED`. `validate()` catches JSON/schema/codec errors and returns an invalidity result instead of guessing.

Atomic text publication must follow the existing temp-file-and-replace pattern used in `pipeline.py`; output paths must resolve inside `<outdir>` and be stored relative to it.

- [ ] **Step 7: Run focused and regression tests**

Run:

```bash
python -m pytest -q tests/test_checkpointing.py
python -m pytest -q tests/test_minimal_run.py
```

Expected: PASS.

- [ ] **Step 8: Commit**

```bash
git add qpcr_pipeline/checkpointing.py tests/test_checkpointing.py
git commit -m "feat: add trustworthy stage checkpoints"
```

---

### Task 2: Add the explicit stage DAG and execution-policy planner

**Files:**
- Create: `qpcr_pipeline/execution.py`
- Create: `tests/test_execution.py`

**Interfaces:**
- Consumes: `ExecutionPolicy` plus read-only checkpoint probe results supplied by the runtime.
- Produces:
  - `STAGE_ORDER = ("input", "qc", "clustering", "alignment", "conservation", "primer_design", "inclusivity", "specificity", "ranking")`.
  - `STAGE_DEPENDENCIES` exactly matching the approved spec, including `inclusivity -> (primer_design, qc)` and `ranking -> (primer_design, inclusivity, specificity)`.
  - `StageAction = Literal["RUN", "REUSE", "FORCED"]`.
  - `ExecutionPolicy(resume=False, from_step=None, force_step=None)` with validation.
  - `StageDecision(stage, action, reason)`.
  - `transitive_descendants(stage)`.
  - `required_reuse_boundary(stage)` for strict `--from-step`.
  - `plan_from_validity(policy, reusable: Mapping[str, bool]) -> tuple[StageDecision, ...]` for deterministic graph semantics; the pipeline runtime is responsible for obtaining validity in topological order.

- [ ] **Step 1: Write failing DAG tests**

```python
from qpcr_pipeline.execution import required_reuse_boundary, transitive_descendants


def test_specificity_descendants_only_include_ranking():
    assert transitive_descendants("specificity") == ("ranking",)


def test_inclusivity_force_does_not_force_specificity():
    assert set(transitive_descendants("inclusivity")) == {"ranking"}


def test_from_specificity_requires_inclusivity_branch():
    required = set(required_reuse_boundary("specificity"))
    assert "inclusivity" in required
    assert "primer_design" in required
    assert "qc" in required
    assert "specificity" not in required
    assert "ranking" not in required
```

- [ ] **Step 2: Run and verify RED**

Run: `python -m pytest -q tests/test_execution.py`

Expected: missing module/functions.

- [ ] **Step 3: Implement DAG traversal and policy validation**

`ExecutionPolicy` rejects:

```text
resume + from_step
from_step + force_step
force_step without resume
unknown stage name
```

Traversal must always return stages in `STAGE_ORDER`, never set iteration order.

- [ ] **Step 4: Add failing planner tests**

```python
def test_normal_run_runs_every_stage(): ...
def test_resume_reuses_all_valid_stages(): ...
def test_resume_invalid_specificity_forces_only_specificity_and_ranking(): ...
def test_resume_invalid_clustering_forces_entire_dependent_chain(): ...
def test_resume_invalid_ranking_runs_only_ranking(): ...
def test_force_inclusivity_reuses_specificity_and_forces_ranking(): ...
def test_from_step_specificity_runs_specificity_and_ranking_only(): ...
```

For specificity invalidation, assert the exact actions:

```python
assert actions == {
    "input": "REUSE",
    "qc": "REUSE",
    "clustering": "REUSE",
    "alignment": "REUSE",
    "conservation": "REUSE",
    "primer_design": "REUSE",
    "inclusivity": "REUSE",
    "specificity": "RUN",
    "ranking": "FORCED",
}
```

- [ ] **Step 5: Implement deterministic planner semantics**

For `--resume`, an invalid stage gets `RUN`; all transitive dependents are promoted to `FORCED` without consulting their stale checkpoints. Independent branches remain eligible for `REUSE`. For `--force-step`, the selected stage and descendants are `FORCED`; other stages follow resume validity. For `--from-step`, selected stage/descendants are `FORCED` and all required boundary nodes are `REUSE` only.

- [ ] **Step 6: Run focused tests**

Run: `python -m pytest -q tests/test_execution.py`

Expected: PASS.

- [ ] **Step 7: Commit**

```bash
git add qpcr_pipeline/execution.py tests/test_execution.py
git commit -m "feat: plan selective pipeline execution"
```

---

### Task 3: Add strict typed state codecs for every stage

**Files:**
- Create: `qpcr_pipeline/checkpoint_codecs.py`
- Create: `tests/test_checkpoint_codecs.py`

**Interfaces:**
- Consumes current domain/result dataclasses from `local_input`, `qc`, `clustering`, `alignment`, `conservation`, `primer_design`, `inclusivity`, `specificity`, and `ranking`.
- Produces codec singletons/functions:
  - `INPUT_CODEC`
  - `QC_CODEC`
  - `CLUSTERING_CODEC`
  - `ALIGNMENT_CODEC`
  - `CONSERVATION_CODEC`
  - `PRIMER_DESIGN_CODEC`
  - `INCLUSIVITY_CODEC`
  - `SPECIFICITY_CODEC`
  - `RANKING_CODEC`
- All generated paths are encoded relative to `<outdir>` and decoded back under that same output root.

- [ ] **Step 1: Write failing path and QC round-trip tests**

Create representative `LocalSequenceRecord`, `QCRecord`, `QCResult`, `TargetSequenceSet`, and `EvaluationSet` values. Assert `decode(encode(value)) == value` for QC and that an encoded generated path never contains the absolute temp directory.

- [ ] **Step 2: Run and verify RED**

Run: `python -m pytest -q tests/test_checkpoint_codecs.py`

Expected: missing codecs.

- [ ] **Step 3: Implement strict primitive helpers and input-record codec**

Use helpers that explicitly require JSON shapes/types rather than permissive `dict.get()` defaults. For `LocalSequenceRecord`, persist `sequence_id`, `sequence`, and only metadata required by current downstream behavior. Preserve Biopython feature annotations as JSON:

```json
{
  "type": "CDS",
  "parts": [{"start": 0, "end": 100, "strand": 1, "ref": null}],
  "qualifiers": {"gene": ["example"], "locus_tag": [], "product": []}
}
```

Decode to `SeqFeature` with `SimpleLocation` or `CompoundLocation`, preserving the fields consumed by conservation annotation extraction. Do not attempt to serialize arbitrary Biopython/Python objects.

- [ ] **Step 4: Add failing round-trip tests for clustering/alignment/conservation**

Construct small valid values for:

```text
ClusteringResult + SequenceCluster + ClusterMember
AlignmentResult + AlignedSequence + AlignmentCoordinate
ConservationResult + PositionConservation + WindowConservation + ReferenceAnnotation
```

Include both `COMPLETE` and `SKIPPED` path shapes.

- [ ] **Step 5: Implement clustering/alignment/conservation codecs**

Every constructor field must be decoded explicitly. Tuple fields remain tuples after decode; optional paths remain `None`; result paths are resolved against `<outdir>`.

- [ ] **Step 6: Add failing round-trip tests for primer/inclusivity/specificity/ranking**

Cover all nested public dataclasses required by downstream stages:

```text
CandidateRegion, DesignedOligo, AssayCandidate, PrimerDesignResult
OligoMatch, ProposedOligoCompatibility, AssayInclusivity, OligoVariation, DegeneracyProposal, InclusivityResult
OffTargetHit, HitRetentionSummary, PlausibleAmplicon, SpecificityResult
RankingReason, ScoreComponents, RankedAssay, RankingResult
```

For `RankingReason.evidence`, verify primitive values survive JSON round-trip and tuple ordering remains deterministic.

- [ ] **Step 7: Implement the remaining codecs and malformed-state rejection**

Reject unknown/missing required top-level fields, wrong collection types, invalid enum/literal values, and paths that escape `<outdir>`.

- [ ] **Step 8: Run focused codec tests**

Run: `python -m pytest -q tests/test_checkpoint_codecs.py`

Expected: PASS.

- [ ] **Step 9: Commit**

```bash
git add qpcr_pipeline/checkpoint_codecs.py tests/test_checkpoint_codecs.py
git commit -m "feat: serialize checkpoint stage state"
```

---

### Task 4: Define stage checkpoint metadata, input identities, and tool identities

**Files:**
- Create: `qpcr_pipeline/checkpoint_stages.py`
- Create: `tests/test_checkpoint_stages.py`
- Modify: `qpcr_pipeline/primer_design.py`

**Interfaces:**
- Consumes: `PipelineConfig`, current upstream stage results, output root, and injectable tool-identity provider.
- Produces:
  - `StageCheckpointDefinition(name, dependencies, codec, parameters, inputs, tools, outputs)` registry for all nine stages.
  - `geison_version() -> str` using `importlib.metadata.version("geison-qpcr")` behind an injectable/patchable function.
  - `SubprocessToolIdentityProvider.identity(tool_name) -> Mapping[str, str]`.
  - `stage_request(stage, config, dependency_manifests, stage_context, tool_provider) -> CheckpointRequest`.
  - `stage_outputs(stage, result, outdir) -> tuple[Path, ...]`.
  - `primer3_required(conservation, config) -> bool` in `primer_design.py`.

- [ ] **Step 1: Write failing relevant-config isolation tests**

Build two `PipelineConfig` objects differing only in specificity configuration. Assert the generated parameter payloads for `input`, `qc`, `clustering`, `alignment`, `conservation`, `primer_design`, and `inclusivity` are equal, while `specificity` differs. Repeat with a clustering-only change and assert clustering differs.

Also assert ranking weights affect only ranking parameters.

- [ ] **Step 2: Run and verify RED**

Run: `python -m pytest -q tests/test_checkpoint_stages.py`

Expected: missing registry/projection implementation.

- [ ] **Step 3: Implement explicit per-stage parameter projections**

Use dataclass-to-primitive conversion only for the stage's own config object. Include `target_name` only where current stage outputs embed it (`conservation` and `ranking`). Do not include `outdir`.

For specificity, separate off-target file contents from parameters: dataset names/order/source mode belong in parameters; actual bytes/manifests belong in stage-local input identities.

- [ ] **Step 4: Add failing local/frozen/off-target identity tests**

Assert:

```text
same local path + changed bytes -> input identity changes
same bytes at a different local path -> scientific input identity remains equal
frozen NCBI records/manifest change -> input identity changes
same off-target config + changed FASTA bytes -> specificity input identity changes
```

The frozen input resolver may call existing `validate_frozen_dataset`; specificity may use existing offline dataset loading/provenance and must never access the network.

- [ ] **Step 5: Implement external input identities**

Local target identity uses input mode/format plus SHA-256 of source bytes. Frozen NCBI identity uses validated frozen `records.gb` plus manifest identity. Online NCBI causal identity uses only acquisition request/configuration; materialized dataset bytes enter the input stage's completed result fingerprint through checkpoint-declared outputs/state.

- [ ] **Step 6: Add failing tool-identity and `primer3_required` tests**

Use a fake provider returning stable values:

```python
class FakeToolIdentityProvider:
    def __init__(self, versions):
        self.versions = versions
        self.calls = []

    def identity(self, tool_name):
        self.calls.append(tool_name)
        return {"name": tool_name, "version": self.versions[tool_name]}
```

Assert:

```text
disabled clustering does not request CD-HIT identity
enabled clustering with non-empty Evaluation Set requests it
alignment with one discovery sequence does not request MAFFT identity
alignment with >1 sequence requests it
disabled primer design does not request Primer3 identity
enabled primer design with no candidate region does not request Primer3 identity
enabled primer design with candidates requests it
changing one tool identity changes only that stage request and descendants later via result fingerprints
```

- [ ] **Step 7: Refactor primer candidate-selection decision minimally**

Add a pure public helper:

```python
def primer3_required(conservation: ConservationResult, config: PrimerDesignConfig) -> bool:
    validate_primer_design_config(config)
    if not config.enabled:
        return False
    return bool(_select_candidate_regions(conservation, config))
```

Keep `design_primers()` on the same `_select_candidate_regions()` implementation so checkpoint identity cannot drift from actual execution behavior.

- [ ] **Step 8: Implement subprocess tool identity provider**

Resolve binaries with `shutil.which`. Capture normalized version output with bounded size and no shell:

```text
cd-hit-est: `cd-hit-est -h`, extract `CD-HIT version ...`
mafft: `mafft --version`
primer3_core: `primer3_core --version`
```

Return a stable mapping containing tool name and normalized version. Absolute executable paths may be retained for diagnostics but must not be part of the canonical fingerprint payload. Missing required binaries raise the existing stage-specific error path before execution.

- [ ] **Step 9: Add failing durable-output ownership tests**

Assert result-owned outputs are declared and that:

```text
conservation never declares root report.html
ranking SKIPPED does not declare root report.html
ranking COMPLETE declares root report.html when html_report_path is not None
None optional paths are omitted
all declared outputs resolve inside outdir
```

- [ ] **Step 10: Implement stage output ownership**

Use result dataclass path fields instead of duplicating artifact filenames where possible. QC/local input may legitimately have no separate scientific artifact beyond checkpoint state. Online NCBI input declares materialized consolidated records/manifest plus the effective copied NCBI manifest; local input does not copy the user's source into checkpoint outputs.

- [ ] **Step 11: Run focused and primer regression tests**

Run:

```bash
python -m pytest -q tests/test_checkpoint_stages.py
python -m pytest -q tests/test_primer_design.py
```

Expected: PASS.

- [ ] **Step 12: Commit**

```bash
git add qpcr_pipeline/checkpoint_stages.py qpcr_pipeline/primer_design.py tests/test_checkpoint_stages.py
git commit -m "feat: define checkpoint stage identities"
```

---

### Task 5: Route `run_pipeline()` through checkpoints and the execution DAG

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Create: `tests/test_pipeline_resume.py`
- Modify: `tests/test_minimal_run.py`

**Interfaces:**
- Consumes: existing `PipelineConfig`, `outdir`, injected NCBI/scientific runners, plus new `ExecutionPolicy` and optional tool identity provider.
- Produces:
  - `run_pipeline(..., execution: ExecutionPolicy | None = None, tool_identity_provider=None) -> RunSummary`.
  - `RunSummary.stage_actions: tuple[StageActionSummary, ...]` or equivalent JSON-compatible deterministic field while preserving existing summary fields.
  - Internal topological execution that returns the same scientific result objects whether a stage is freshly run or decoded from a checkpoint.

- [ ] **Step 1: Write a failing normal-run checkpoint test**

With a minimal local FASTA and disabled/default expensive stages, call `run_pipeline()` normally and assert all nine checkpoint manifests exist and are `COMPLETE`. Assert normal run actions are all `RUN`, even if a previous complete checkpoint exists.

- [ ] **Step 2: Run and verify RED**

Run: `python -m pytest -q tests/test_pipeline_resume.py::test_normal_run_writes_all_stage_checkpoints_and_does_not_reuse`

Expected: missing execution/checkpoint integration.

- [ ] **Step 3: Split existing orchestration into stage-local helpers without changing science**

Inside `pipeline.py`, extract the current sequential bodies into focused private helpers with these conceptual outputs:

```text
_run_input -> tuple[LocalSequenceRecord, ...] + input artifact context
_run_qc -> QCResult
_run_clustering -> ClusteringResult
_run_alignment -> AlignmentResult
_run_conservation -> ConservationResult
_run_primer_design -> PrimerDesignResult
_run_inclusivity -> InclusivityResult
_run_specificity -> SpecificityResult
_run_ranking -> RankingResult
```

The helper bodies must call the same existing scientific functions with the same arguments as today's `run_pipeline()`.

For online NCBI:

- normal full execution and explicit `--resume --force-step input` intentionally refresh by removing/recreating only `<outdir>/ncbi_dataset` before acquisition;
- an input stage rerun caused by interrupted/invalid checkpoint during plain `--resume` preserves `<outdir>/ncbi_dataset` so existing batch-level NCBI acquisition can continue its partial download;
- frozen external datasets are never deleted or modified.

- [ ] **Step 4: Implement stage execute/publish/reuse loop**

For each planned stage:

```text
REUSE  -> CheckpointManager.validate + load typed state
RUN    -> begin -> scientific helper -> complete
FORCED -> begin -> scientific helper -> complete
failure -> best-effort fail -> re-raise original exception
```

For newly executed stages, direct dependency references use the newly completed manifest `result_fingerprint`. For reused stages they use the validated existing manifest `result_fingerprint`.

- [ ] **Step 5: Add failing all-valid resume test**

Run once normally, then run `execution=ExecutionPolicy(resume=True)` with fake runners whose `run()` methods raise if invoked. Assert all stage actions are `REUSE` and the final scientific summary equals the first run.

- [ ] **Step 6: Implement read-only resume probing and reuse**

Probe stages topologically. When one checkpoint is invalid, mark it `RUN` and promote only its transitive descendants to `FORCED`; do not inspect/reuse stale descendants. Preserve independent branches.

- [ ] **Step 7: Add failing strict `--from-step` preflight tests**

Create a complete run, corrupt/delete the inclusivity state, then request `from_step="specificity"`. Assert the call raises before specificity/ranking artifacts change. Error text must include `inclusivity`, invalidity category, and guidance to use `--resume` or restart earlier.

Also test a valid `from_step="specificity"` executes specificity+ranking while reusing required boundary stages.

- [ ] **Step 8: Implement strict from-step preflight**

Validate and decode every stage in `required_reuse_boundary(from_step)` before calling `CheckpointManager.begin()` for the selected stage. Aggregate blocking stages into one deterministic error ordered by `STAGE_ORDER`.

- [ ] **Step 9: Add failing run-summary diagnostics test**

Assert `run_summary.json` contains existing fields plus ordered stage actions:

```json
"stage_actions": [
  {"stage": "input", "action": "REUSE"},
  {"stage": "specificity", "action": "RUN"},
  {"stage": "ranking", "action": "FORCED"}
]
```

The actual list contains all nine stages, in stage order.

- [ ] **Step 10: Preserve aggregate QC report behavior**

After all stage results are either run or reused, construct `qc_report.json` from the reconstructed/current result objects exactly as today. This guarantees resume does not change the report contract.

- [ ] **Step 11: Run focused and full minimal pipeline regression tests**

Run:

```bash
python -m pytest -q tests/test_pipeline_resume.py
python -m pytest -q tests/test_minimal_run.py
```

Expected: PASS.

- [ ] **Step 12: Commit**

```bash
git add qpcr_pipeline/pipeline.py tests/test_pipeline_resume.py tests/test_minimal_run.py
git commit -m "feat: resume checkpointed pipeline runs"
```

---

### Task 6: Expose the approved CLI contract

**Files:**
- Modify: `qpcr_pipeline/cli.py`
- Modify: `tests/test_cli.py`

**Interfaces:**
- Consumes CLI options.
- Produces validated `ExecutionPolicy` passed to `run_pipeline()`.

- [ ] **Step 1: Write failing parser/contract tests**

Add subprocess or direct-parser tests for:

```text
--resume accepted with --outdir
--from-step specificity accepted
--resume --force-step specificity accepted
--resume + --from-step rejected
--from-step + --force-step rejected
--force-step without --resume rejected
unknown stage rejected by argparse choices
resume/from/force rejected when --outdir is absent because no run is performed
```

- [ ] **Step 2: Run and verify RED**

Run: `python -m pytest -q tests/test_cli.py`

Expected: unrecognized flags / missing validation.

- [ ] **Step 3: Implement flags with registry-backed choices**

Add:

```python
run_parser.add_argument("--resume", action="store_true")
run_parser.add_argument("--from-step", choices=STAGE_ORDER)
run_parser.add_argument("--force-step", choices=STAGE_ORDER)
```

Construct `ExecutionPolicy` only for actual runs with `--outdir`; report policy validation through `parser.error(...)` so invalid usage exits non-zero with actionable text.

- [ ] **Step 4: Run CLI tests**

Run: `python -m pytest -q tests/test_cli.py`

Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/cli.py tests/test_cli.py
git commit -m "feat: expose pipeline resume controls"
```

---

### Task 7: Prove interruption recovery and selective invalidation end to end

**Files:**
- Modify: `tests/test_pipeline_resume.py`

**Interfaces:**
- Consumes implemented checkpointed pipeline.
- Produces acceptance-level tests corresponding one-to-one with issue #11 criteria.

- [ ] **Step 1: Write the interruption/resume test**

Use fake runners and a deterministic failing stage boundary. A first run must complete upstream checkpoints and fail during a later stage. A second run with `ExecutionPolicy(resume=True)` must assert upstream fake runner call counts do not increase, the interrupted stage runs again, downstream stages complete, and final `run_summary.json` is `COMPLETED`.

Use an injected failure that does not require real CD-HIT/MAFFT/Primer3 binaries.

- [ ] **Step 2: Run and verify the new test exposes any missing behavior**

Run: `python -m pytest -q tests/test_pipeline_resume.py -k interruption`

Expected before final fixes: FAIL if interrupted stage/status/reuse handling is incomplete.

- [ ] **Step 3: Fix only the interruption/resume defects revealed by the test**

Ensure a catchable exception leaves the failed stage non-reusable while already completed upstream manifests remain valid. Do not delete valid upstream state during resume.

- [ ] **Step 4: Add specificity-only invalidation acceptance test**

Run once with specificity config A, then resume with config B differing only in a specificity parameter. Assert:

```text
input REUSE
qc REUSE
clustering REUSE
alignment REUSE
conservation REUSE
primer_design REUSE
inclusivity REUSE
specificity RUN
ranking FORCED
```

Use fake runner call counters to independently prove alignment/conservation were not recalculated.

- [ ] **Step 5: Add clustering-chain invalidation acceptance test**

Run once, alter only `ClusteringConfig.identity`, then resume. Assert `input` and `qc` are reused while clustering and every transitive dependent are `RUN`/`FORCED` as appropriate.

- [ ] **Step 6: Add force-branch acceptance test**

For `ExecutionPolicy(resume=True, force_step="inclusivity")`, assert specificity is reused while inclusivity and ranking are forced.

- [ ] **Step 7: Add output-integrity acceptance test**

Corrupt `alignment/discovery_alignment.fasta`, resume, and assert alignment plus transitive dependents rerun while input/qc/clustering remain reused.

- [ ] **Step 8: Add software/tool-version invalidation tests**

Patch Geison version from `A` to `B` and assert all checkpoints invalidate. Change only fake MAFFT identity and assert alignment plus dependents invalidate, while upstream remains reusable. Verify disabled/unused tool identities are not queried.

- [ ] **Step 9: Add online acquisition result-fingerprint regression**

With a fake NCBI client, perform an acquisition request, record its input checkpoint result fingerprint, explicitly force/refresh input with the same request but different returned record bytes, and assert the new input `result_fingerprint` differs and downstream stages are forced.

- [ ] **Step 10: Run all issue #11 focused tests**

Run:

```bash
python -m pytest -q \
  tests/test_checkpointing.py \
  tests/test_execution.py \
  tests/test_checkpoint_codecs.py \
  tests/test_checkpoint_stages.py \
  tests/test_pipeline_resume.py \
  tests/test_cli.py
```

Expected: PASS.

- [ ] **Step 11: Commit**

```bash
git add tests/test_pipeline_resume.py qpcr_pipeline
git commit -m "test: verify checkpoint resume semantics"
```

---

### Task 8: Document the workflow and run final regression verification

**Files:**
- Modify: `README.md`

**Interfaces:**
- Consumes the completed CLI/runtime contract.
- Produces concise user documentation and final verification evidence.

- [ ] **Step 1: Add README examples and semantics**

Document exactly:

```bash
qpcr-pipeline run config.yaml --outdir run1
qpcr-pipeline run config.yaml --outdir run1 --resume
qpcr-pipeline run config.yaml --outdir run1 --from-step conservation
qpcr-pipeline run config.yaml --outdir run1 --resume --force-step specificity
```

Explain normal full recomputation, strict `--from-step`, force semantics, `.checkpoints/` being internal, hash validation, and that copying a complete output directory preserves its local resume state.

- [ ] **Step 2: Run the complete normal test suite**

Run: `python -m pytest -q`

Expected: PASS.

- [ ] **Step 3: Run syntax/import validation explicitly**

Run:

```bash
python -m compileall -q qpcr_pipeline tests
python -c "from qpcr_pipeline.cli import build_parser; from qpcr_pipeline.pipeline import run_pipeline"
```

Expected: exit 0.

- [ ] **Step 4: Review the branch diff against the spec**

Check that implementation changes are limited to issue #11, all seven GitHub acceptance criteria have direct tests, no shared root report ownership regression was introduced, and `.circleci/config.yml` has no permanent branch-filter change.

- [ ] **Step 5: Commit documentation**

```bash
git add README.md
git commit -m "docs: explain checkpoint resume workflow"
```

- [ ] **Step 6: Verify the exact final branch SHA**

Run `python -m pytest -q` on the exact candidate merge SHA. If this environment cannot execute the repository locally, obtain equivalent RED/GREEN/final evidence from the existing CircleCI test job without leaving a feature-branch filter change in the final diff.

Expected: all normal tests PASS on the exact candidate merge SHA.

---

## Acceptance Criteria Traceability

| GitHub #11 criterion | Primary implementation task | Primary verification |
| --- | --- | --- |
| Manifest per stage with inputs, relevant params, tool version, outputs, status | Tasks 1, 4, 5 | `test_checkpointing.py`, `test_checkpoint_stages.py`, normal pipeline checkpoint test |
| `--resume` reuses valid stages | Tasks 2, 5 | all-valid resume + interruption tests |
| `--from-step` continues from chosen stage | Tasks 2, 5, 6 | strict boundary/preflight tests |
| `--force-step` recalculates chosen stage and invalidates dependents | Tasks 2, 5, 6 | force-inclusivity branch test |
| Specificity parameter does not recalc alignment/conservation | Tasks 4, 7 | specificity-only invalidation test |
| Clustering parameter invalidates full dependent chain | Tasks 2, 4, 7 | clustering-chain invalidation test |
| Tests simulate interruption and resume | Tasks 5, 7 | interruption/resume acceptance test |
