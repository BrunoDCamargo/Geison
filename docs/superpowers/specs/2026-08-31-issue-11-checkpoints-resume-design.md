# Issue #11 design: checkpoints, resume and selective invalidation

## Status

Approved design for GitHub issue #11, "Adicionar checkpoints, resume e invalidação seletiva".

## Problem

The pipeline currently executes the full scientific chain in one sequential `run_pipeline()` call. Long runs cannot resume safely after interruption, and changing a late-stage parameter can force expensive upstream recomputation.

Issue #11 requires stage manifests, dependency hashes, `--resume`, `--from-step`, `--force-step`, and selective invalidation.

## Goals

1. Persist a trustworthy checkpoint after every stage.
2. Reuse only checkpoints proven valid for the current inputs, relevant configuration, software identity and dependency state.
3. Resume an interrupted run without recomputing valid stages.
4. Allow explicit restart from one stage.
5. Allow explicit recomputation of one stage plus only its transitive dependents.
6. Preserve scientific stage isolation: checkpoint infrastructure must not be embedded throughout the scientific algorithms.
7. Keep checkpoint state local to the selected output directory.
8. Make invalidation deterministic and testable.

## Non-goals

- Global/shared cache across output directories.
- Cross-project cache discovery.
- Workflow-engine dependency.
- `pickle` or other Python-runtime-specific checkpoint serialization.
- Multiple forced stages in one command.
- Compatibility promises across Geison versions.
- Reproduction metadata beyond what is required to validate checkpoints; broader execution diagnostics belong to issue #12.

## Stage graph

The checkpoint graph is explicit rather than inferred from call order.

```text
input
  |
  v
 qc ---------------------------+
  |                             |
  v                             |
clustering                       |
  |                             |
  v                             |
alignment                         |
  |                             |
  v                             |
conservation                      |
  |                             |
  v                             |
primer_design                     |
  |\                             |
  | \                            |
  v  v                           v
inclusivity <--------------------+   specificity
  \                                  /
   \                                /
    +-------------> ranking <-------+
```

Direct dependencies:

| Stage | Direct dependencies |
| --- | --- |
| `input` | none |
| `qc` | `input` |
| `clustering` | `qc` |
| `alignment` | `clustering` |
| `conservation` | `alignment` |
| `primer_design` | `conservation` |
| `inclusivity` | `primer_design`, `qc` |
| `specificity` | `primer_design` |
| `ranking` | `primer_design`, `inclusivity`, `specificity` |

`inclusivity` depends on both the approved evaluation sequence set produced by QC and the assay definitions produced by primer design. Declaring both dependencies prevents a QC change from leaving inclusivity incorrectly reusable through a superficially unchanged primer-design result.

External off-target inputs are stage-local inputs of `specificity`, not graph stages.

## Architecture

### Pipeline orchestrator

`run_pipeline()` remains the high-level scientific entry point but delegates checkpoint decisions to an execution planner and `CheckpointManager`.

The orchestrator owns:

- stage order;
- dependency graph;
- execution versus reuse decisions;
- reconstruction of stage results from checkpoint snapshots;
- publication of final run diagnostics.

Scientific functions such as clustering, alignment, conservation and specificity remain focused on their scientific work and do not implement cache policy.

### Stage definitions

A small stage registry declares, for each stage:

- name;
- direct dependencies;
- relevant configuration projection;
- external-input identity resolver, when applicable;
- external tool identity resolver, when applicable;
- typed JSON state codec;
- durable stage-owned outputs that participate in checkpoint validation.

The registry is the single source of truth for dependency traversal and selective invalidation.

### CheckpointManager

`CheckpointManager` is responsible for:

- deterministic fingerprint calculation;
- manifest paths;
- atomic manifest and state publication;
- SHA-256 calculation;
- checkpoint validation;
- state loading;
- state saving;
- invalidation reasons.

It must not know scientific domain rules beyond data supplied by the stage definitions.

### Execution planner

The planner converts CLI intent plus checkpoint validity into one deterministic action per stage.

Actions are:

- `RUN`: execute the stage and write a new checkpoint;
- `REUSE`: load a valid checkpoint;
- `FORCED`: execute because the user explicitly forced the stage or it is a transitive dependent of a forced/invalid stage.

`FORCED` is operationally a run but remains distinct in diagnostics.

For `--from-step`, the planner computes the complete plan and validates every required reusable boundary checkpoint before scientific execution, so an invalid prerequisite fails before any stage output is modified.

For `--resume`, once a stage is invalid, that stage and all transitive dependents are scheduled to run even if an old downstream checkpoint happens to remain on disk. The planner does not wait to see whether the rerun produces byte-identical results before deciding to reuse descendants.

## Checkpoint layout

Checkpoints are scoped to the current output directory.

```text
<outdir>/
  .checkpoints/
    input/
      manifest.json
      state.json
    qc/
      manifest.json
      state.json
    clustering/
      manifest.json
      state.json
    alignment/
      manifest.json
      state.json
    conservation/
      manifest.json
      state.json
    primer_design/
      manifest.json
      state.json
    inclusivity/
      manifest.json
      state.json
    specificity/
      manifest.json
      state.json
    ranking/
      manifest.json
      state.json
```

There is no global cache in issue #11.

## Typed state snapshots

Every stage stores a deterministic JSON snapshot sufficient to reconstruct the Python result required by downstream stages.

Each stage has an explicit encode/decode codec. The codec is versioned through the checkpoint schema and must reject malformed or incomplete state rather than guessing defaults.

`pickle` is forbidden because it is opaque, unsafe for untrusted data, and tightly coupled to Python implementation details.

Published scientific artifacts remain in their existing user-facing locations. `state.json` is checkpoint infrastructure and is not a replacement for those reports.

Paths stored in checkpoint state should be represented relative to `<outdir>` whenever they refer to generated files inside the run directory.

## Manifest model

A completed manifest has this conceptual structure:

```json
{
  "schema_version": 1,
  "stage": "specificity",
  "status": "COMPLETE",
  "fingerprint": "sha256:causal-inputs-and-policy",
  "result_fingerprint": "sha256:causal-fingerprint-plus-published-state",
  "dependencies": {
    "primer_design": "sha256:dependency-result-fingerprint"
  },
  "inputs": {
    "off_targets": "sha256:..."
  },
  "parameters": {
    "...": "only parameters relevant to this stage"
  },
  "software": {
    "geison": "0.1.0.dev0"
  },
  "tools": {},
  "outputs": [
    {
      "path": "specificity/specificity_report.json",
      "sha256": "..."
    }
  ],
  "state": {
    "path": ".checkpoints/specificity/state.json",
    "sha256": "..."
  }
}
```

Manifest statuses are:

- `RUNNING`: stage execution started but is not reusable;
- `COMPLETE`: stage and checkpoint publication finished successfully and may be reusable;
- `FAILED`: stage failed and is not reusable.

Only `COMPLETE` is eligible for reuse.

A killed process may leave `RUNNING`; this is intentionally invalid.

## Atomic publication

Before a stage begins, its manifest is atomically replaced with a `RUNNING` manifest. This intentionally supersedes any previous reusable manifest for that stage.

On success:

1. the scientific stage completes its own output publication;
2. the typed `state.json` is written atomically;
3. hashes of all checkpoint-declared outputs and state are calculated;
4. `result_fingerprint` is calculated from the causal `fingerprint` plus state/output identities;
5. the final `COMPLETE` manifest is written atomically.

On a catchable stage exception, the manager attempts to publish `FAILED` and then re-raises the original failure. A failure to publish the diagnostic manifest must not hide the original scientific exception.

If the process dies between steps, no `COMPLETE` manifest is created and the stage cannot be reused.

## Fingerprints

Two deterministic hashes have distinct roles.

### Causal fingerprint

`fingerprint` answers: "Would the same stage be requested under the current scientific inputs and policy?"

It is SHA-256 over canonical deterministic JSON containing:

- checkpoint schema version;
- stage name;
- Geison version;
- relevant stage parameters only;
- direct dependency `result_fingerprint` values;
- identities of stage-local external inputs that are knowable without performing the stage;
- external tool identity when that tool participates in the enabled stage.

A checkpoint can be reused only when its stored causal fingerprint matches the currently expected causal fingerprint.

### Result fingerprint

`result_fingerprint` answers: "Which concrete completed result did this stage produce?"

It is SHA-256 over canonical JSON containing:

- the causal `fingerprint`;
- the checkpoint-state SHA-256;
- the sorted durable output path/SHA-256 pairs.

Downstream dependency fingerprints always reference the dependency's `result_fingerprint`, not merely its causal fingerprint. Therefore a stage that legitimately produces different bytes under the same request, such as a refreshed online acquisition, necessarily invalidates downstream checkpoints.

Canonical serialization uses stable key ordering and excludes timestamps, runtime action (`RUN`/`REUSE`), absolute output-directory location and other incidental metadata.

### Relevant configuration projections

Each stage hashes only configuration that can affect that stage's result.

Examples:

- specificity parameters affect `specificity` and therefore `ranking`, but not alignment or conservation;
- clustering identity or other clustering parameters affect `clustering` and all transitive dependents;
- ranking weights affect only `ranking`;
- output-directory location does not by itself invalidate scientific stages.

The projections are explicit and covered by tests. The implementation must not hash the entire `PipelineConfig` for every stage because that would defeat selective invalidation.

## Input identity

### Local FASTA/GenBank

The `input` causal fingerprint includes the selected input mode, relevant parsing configuration and SHA-256 identity of the source bytes used by the run.

A changed local input file therefore invalidates `input` and every dependent stage.

The `input` state snapshot contains the normalized records required to resume downstream processing without relying on in-memory objects from the original process.

### Frozen NCBI dataset

The `input` causal fingerprint includes the frozen dataset identity/manifest and relevant acquisition configuration. Existing frozen-dataset immutability and output-directory safety rules remain authoritative.

### Online NCBI acquisition

Online acquisition needs separate request and result identity because NCBI may return a different concrete dataset for the same query over time.

For reuse, the `input` causal fingerprint includes the acquisition request/configuration identity and the existing materialized checkpoint is validated by its saved output/state hashes. `--resume` does not query NCBI merely to check freshness.

When online acquisition is actually run, the newly materialized dataset hashes participate in `result_fingerprint`. If refreshed acquisition returns different data while the request is unchanged, the new `result_fingerprint` changes and all dependent stages are invalidated.

To intentionally refresh acquisition, the user can run normally or use `--resume --force-step input`.

## External software identity

Changing Geison version invalidates all checkpoints. Issue #11 does not attempt cross-version compatibility analysis.

External tools invalidate only the stage that invokes them and its dependents:

- CD-HIT-EST identity belongs to `clustering` when clustering is enabled;
- MAFFT identity belongs to `alignment` when alignment is enabled;
- Primer3 identity belongs to `primer_design` when primer design invokes Primer3.

Tool identity acquisition must be behind an injectable/testable boundary so unit tests do not require installed binaries. Real subprocess-backed runners must provide the detected executable identity/version used for the run.

A disabled stage must not become invalid merely because an unused external tool version changed.

## Output integrity

A `COMPLETE` checkpoint is reusable only when:

1. its schema and stage name are valid;
2. its expected causal fingerprint matches the current expected fingerprint;
3. every declared output exists;
4. every declared output SHA-256 matches;
5. `state.json` exists;
6. the state SHA-256 matches;
7. the stored `result_fingerprint` matches a recomputation from the validated state/outputs;
8. the stage codec can decode and validate the state.

Any failure makes the checkpoint invalid. Partial reuse is forbidden.

Checkpoint output lists contain durable outputs owned by that stage. Shared or conditionally overwritten presentation aliases must not create false invalidation.

In particular, the root `report.html` has conditional ownership between conservation/ranking behavior from issue #10. Conservation must not declare that shared root alias as a checkpoint-critical output if a later valid ranking stage may replace it. Ranking may declare the root report only when ranking is the active owner for that configuration. Stable scientific/report artifacts with exclusive ownership remain checkpoint-declared and hash-validated.

## CLI contract

The existing command gains:

```text
qpcr-pipeline run CONFIG --outdir OUTDIR [--resume]
qpcr-pipeline run CONFIG --outdir OUTDIR --from-step STAGE
qpcr-pipeline run CONFIG --outdir OUTDIR --resume --force-step STAGE
```

Valid stage names are the names in the stage registry.

### Normal execution

A run without resume flags:

- recalculates every defined stage;
- does not silently reuse prior checkpoints;
- writes or replaces checkpoints as stages complete.

A scientifically disabled stage still executes its configured disabled behavior and can publish a valid checkpoint.

Any successful normal run is automatically resumable later.

### `--resume`

`--resume` reuses every valid checkpoint and runs only invalid stages plus their transitive dependents.

Example after a specificity-only parameter change:

```text
input          REUSE
qc             REUSE
clustering     REUSE
alignment      REUSE
conservation   REUSE
primer_design  REUSE
inclusivity    REUSE
specificity    RUN
ranking        FORCED
```

`ranking` is marked `FORCED` because a dependency is being recomputed. It is not reused even if its old manifest remains physically present.

### `--from-step STAGE`

`--from-step` is strict and predictable:

- the selected stage is always run;
- all transitive dependents of that stage are run;
- every stage outside that forced subgraph whose result is required by the selected stage or one of its forced descendants must have a valid reusable checkpoint;
- if any required checkpoint is missing or invalid, the command fails before scientific execution begins and identifies the blocking stage(s).

This rule handles the branch after `primer_design`. For example, `--from-step specificity` runs specificity and ranking but ranking still needs a valid inclusivity checkpoint. It does not silently recalculate inclusivity.

### `--force-step STAGE`

`--force-step` requires `--resume`.

The selected stage and all of its transitive dependents are forced to run even if their fingerprints would otherwise match. Other branches may still be reused when valid.

Example:

```text
--resume --force-step inclusivity

input          REUSE
qc             REUSE
clustering     REUSE
alignment      REUSE
conservation   REUSE
primer_design  REUSE
inclusivity    FORCED
specificity    REUSE
ranking        FORCED
```

A forced upstream stage deliberately creates a new causal execution chain for its descendants; descendants are not reused merely because their old fingerprint happens to be equal.

### Flag combinations

- `--resume` and `--from-step` are mutually exclusive.
- `--from-step` and `--force-step` are mutually exclusive.
- `--force-step` without `--resume` is rejected.
- issue #11 supports one `--from-step` or one `--force-step`, not repeated values.

## Disabled scientific stages

A scientifically disabled stage still participates in the execution graph. Its scientific result may be `SKIPPED`, but its checkpoint is `COMPLETE` when the stage correctly executed the configured disabled behavior.

That distinction is important: checkpoint status describes execution integrity, while the scientific result status describes scientific behavior.

## Selective invalidation examples

### Change specificity parameter

```text
input          valid
qc             valid
clustering     valid
alignment      valid
conservation   valid
primer_design  valid
inclusivity    valid
specificity    invalid
ranking        dependent: rerun
```

### Change clustering parameter

```text
input          valid
qc             valid
clustering     invalid
alignment      dependent: rerun
conservation   dependent: rerun
primer_design  dependent: rerun
inclusivity    dependent: rerun
specificity    dependent: rerun
ranking        dependent: rerun
```

### Corrupt one alignment output

Alignment fails output-hash validation. Under `--resume`, alignment runs again and all transitive dependents run again. `input`, `qc` and `clustering` remain reusable.

### Change ranking weights

Only `ranking` becomes invalid and reruns.

## Run diagnostics

`run_summary.json` gains a deterministic stage execution summary so the user can see what happened, for example:

```json
{
  "stage_actions": [
    {"stage": "input", "action": "REUSE"},
    {"stage": "specificity", "action": "RUN"},
    {"stage": "ranking", "action": "FORCED"}
  ]
}
```

The existing scientific summary fields remain preserved. Issue #12 may later add richer environment, duration and reproducibility diagnostics.

Invalid checkpoint errors must identify at minimum:

- stage;
- invalidity category, such as missing manifest, fingerprint mismatch, missing output, output hash mismatch, invalid state or non-complete status;
- what the user can do next when using strict `--from-step`.

## Compatibility and safety

- Existing normal CLI behavior remains full recomputation by default.
- Existing frozen NCBI output-directory guard remains in effect before writes.
- Checkpoint paths are internal to `<outdir>/.checkpoints`.
- Existing stage artifacts and scientific report formats remain authoritative user-facing outputs.
- Checkpoint validation never trusts a manifest path/hash claim without checking the filesystem.
- Malformed checkpoint JSON is treated as invalid state, not as a best-effort partial recovery.

## Testing strategy

Implementation follows TDD. Tests must cover both planning logic and actual pipeline behavior.

### Checkpoint unit tests

- canonical causal/result fingerprint determinism;
- relevant config projection isolation;
- manifest `RUNNING`, `COMPLETE`, `FAILED` lifecycle;
- state/output SHA-256 validation;
- missing/corrupt manifest;
- missing/corrupt state;
- missing/corrupt declared output;
- codec rejection of malformed state;
- tool-version invalidation;
- Geison-version invalidation;
- same online acquisition request with changed materialized data yields a changed `result_fingerprint`.

### Planner tests

- normal execution runs all stages;
- `--resume` reuses all valid stages;
- specificity-only parameter change does not rerun alignment/conservation;
- clustering change invalidates clustering and its complete dependent chain;
- ranking-only change invalidates only ranking;
- `--force-step inclusivity` forces inclusivity + ranking while allowing specificity reuse;
- `--from-step specificity` requires valid inclusivity for ranking;
- invalid prerequisites make `--from-step` fail before any run action;
- invalid CLI flag combinations are rejected.

### Pipeline interruption/resume test

A test stage/fake runner deliberately raises after one or more upstream checkpoints have completed. A subsequent `--resume` run must:

- reuse completed valid upstream stages;
- rerun the interrupted stage;
- execute required downstream stages;
- produce a completed final run.

This directly satisfies the issue criterion requiring simulated interruption and continuation.

### Existing regression suite

All current unit tests remain green. Integration tests on `main` continue to cover real external tools according to the existing CircleCI policy.

## Acceptance mapping

| Issue #11 criterion | Design coverage |
| --- | --- |
| Manifest per stage with inputs, relevant parameters, tool version, outputs and status | Stage registry + manifest model |
| `--resume` reuses valid stages | Execution planner + checkpoint validation |
| `--from-step` continues from chosen stage | Strict `--from-step` contract |
| `--force-step` recalculates selected stage and invalidates dependents | DAG transitive forced subgraph |
| Specificity parameter does not recalculate alignment/conservation | Relevant configuration projections + DAG |
| Clustering parameter invalidates clustering and dependent chain | Dependency result fingerprints + DAG |
| Tests simulate interruption and resume | Pipeline interruption/resume test |

## Implementation boundaries

Expected implementation will introduce focused checkpoint/planning modules and make targeted changes to CLI and pipeline orchestration. Scientific modules should only change where a typed codec or injectable tool identity boundary requires it.

No unrelated scientific refactor belongs in issue #11.
