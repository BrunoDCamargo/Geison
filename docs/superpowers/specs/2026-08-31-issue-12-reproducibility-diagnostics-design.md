# Issue #12 design: reproducibility and execution diagnostics

## Status

Approved design for GitHub issue #12, "Fechar reprodutibilidade e diagnóstico da execução".

## Problem

Geison already persists deterministic stage checkpoints and scientific artifacts, but it does not yet provide one durable account of an entire run. A user cannot inspect one place to learn which effective configuration was used, which environment executed the analysis, which stages ran or were reused, why the result is scientifically complete or partial, or where a failure occurred.

The CLI also lacks an environment diagnostic command and a side-effect-free way to validate and preview an execution before starting it.

Issue #12 adds run-level provenance and diagnostics without duplicating the checkpoint and fingerprint rules implemented by issue #11.

## Goals

1. Persist a run-level manifest that explains what happened, under which effective configuration and environment, and with which final status.
2. Preserve local or NCBI input provenance, resolved accession versions and the selected reference without storing sequence payloads in diagnostic files.
3. Distinguish `COMPLETED`, `PARTIAL` and `FAILED` runs using explicit scientific completeness rules.
4. Guarantee that a partial run never publishes an `IN SILICO PASS` classification.
5. Add `qpcr-pipeline doctor` for environment inspection without requiring a pipeline configuration.
6. Add `--dry-run` for validation and execution planning without scientific or filesystem side effects.
7. Emit structured progress and failure events without exposing secrets or full biological sequences.
8. Preserve all checkpoint, resume, force and from-step behavior from issue #11.

## Non-goals

- Reintroducing BLAST+ as a Geison dependency or specificity backend.
- Replacing stage checkpoint manifests or copying their full fingerprints into the run manifest.
- Building a general workflow engine or observability service.
- Uploading telemetry or logs to an external service.
- Persisting raw FASTA, GenBank, NCBI response bodies or sequence strings in diagnostics.
- Guaranteeing reproducibility across unsupported Geison versions.
- Adding a graphical diagnostics interface; that belongs to later desktop work.
- Adding duration benchmarking, resource profiling or performance optimization.

## Architectural overview

```text
CLI
├── doctor -------------------+
└── run / --dry-run           |
                              v
                     EnvironmentInspector
                     Python / Git / tools
                              |
                              v
ExecutionPolicy ------> ExecutionPlan
                              |
            +-----------------+-----------------+
            |                                   |
            v                                   v
       DryRunReport                    RunRecorder / RunManifest
                                                |
                                     +----------+----------+
                                     |                     |
                                     v                     v
                              run_manifest.json      run.log.jsonl
                                     |
                                     v
                                Pipeline stages
                                     |
                                     +--> existing checkpoints
```

The design introduces three focused boundaries:

- `EnvironmentInspector` describes the runtime and external tools through an injectable command runner.
- `ExecutionPlan` turns `ExecutionPolicy` plus checkpoint state into deterministic stage actions used by both real execution and dry-run.
- `RunRecorder` owns run-level lifecycle, manifest publication, structured events and sanitization.

Scientific stage implementations remain unaware of run-level persistence. The pipeline orchestrator calls these boundaries around existing stages.

## Execution planning

Issue #11 already contains the rules that classify each stage as `RUN`, `REUSE` or `FORCED`. Issue #12 extracts that decision logic into a read-only planner in `execution.py` so dry-run and real execution cannot disagree.

An execution plan contains:

```json
{
  "stages": [
    {"stage": "input", "action": "REUSE", "reason": "CHECKPOINT_VALID"},
    {"stage": "specificity", "action": "RUN", "reason": "FINGERPRINT_MISMATCH"},
    {"stage": "ranking", "action": "FORCED", "reason": "DEPENDENCY_RECOMPUTED"}
  ]
}
```

Reasons are stable diagnostic codes, not free-form parsing of checkpoint errors. Existing detailed invalidity information may accompany a code as a sanitized human message.

Planning remains conservative:

- a normal execution plans every stage as `RUN`;
- `--resume` reuses a valid prefix/branch and forces all transitive dependents after an invalid or forced stage;
- `--from-step` validates the full required reusable boundary before any scientific action;
- `--force-step` requires `--resume` and forces the selected stage and all transitive dependents;
- a disabled scientific stage still participates in the graph and can be planned and checkpointed normally.

The planner may load valid checkpoint state to evaluate downstream checkpoint requests. It must not create directories, write manifests, query NCBI or execute scientific tools.

## Environment inspection

`EnvironmentInspector` returns a typed, serializable snapshot. Subprocess execution is behind an injectable boundary so tests do not require installed scientific binaries.

The snapshot contains:

```json
{
  "python": {
    "status": "USED",
    "required": true,
    "installed": true,
    "version": "3.12.8"
  },
  "geison": {
    "status": "USED",
    "required": true,
    "installed": true,
    "version": "0.1.0.dev0"
  },
  "git": {
    "status": "AVAILABLE",
    "required": false,
    "installed": true,
    "commit": "43b23367...",
    "dirty": false
  },
  "tools": {
    "cd-hit-est": {
      "status": "USED",
      "required": true,
      "installed": true,
      "version": "4.8.1"
    },
    "mafft": {
      "status": "USED",
      "required": true,
      "installed": true,
      "version": "7.x"
    },
    "primer3_core": {
      "status": "USED",
      "required": true,
      "installed": true,
      "version": "2.x"
    },
    "blast+": {
      "status": "NOT_USED",
      "required": false,
      "installed": false,
      "version": null
    }
  }
}
```

Exact version strings are preserved as reported after whitespace normalization. Diagnostic files store command names and reported versions, not full executable paths.

Tool requirement is configuration-aware:

- CD-HIT-EST is required only when clustering is enabled;
- MAFFT is required only when alignment is enabled;
- Primer3 is required only when primer design invokes the external binary;
- BLAST+ is always `NOT_USED` in issue #12, whether installed or absent.

Git metadata is optional. Outside a Git checkout, `commit` is `null`, `dirty` is `null` and status is `UNAVAILABLE`; this does not fail a run. A dirty checkout is recorded but is not itself an execution failure.

The existing checkpoint tool identity provider remains authoritative for checkpoint validity. `EnvironmentInspector` may share its subprocess command runner but must not replace or weaken checkpoint validation.

## Run identity and retry history

An output directory represents one scientific run across interruptions and resumes. The first real invocation creates a `run_id`. Later `--resume`, `--from-step` or `--force-step` invocations in the same output directory retain that `run_id` and append a new attempt.

This prevents a successful resume from erasing evidence of an earlier failure.

Each attempt has:

- `attempt_id`;
- execution policy;
- start and finish timestamps;
- planned and actual stage actions;
- attempt status;
- sanitized failure details when applicable.

A normal non-resume execution targeting an existing run directory starts a new attempt under the existing run identity and records that every stage was recomputed. Issue #12 does not add automatic archival of prior run directories.

## Run manifest model

`<outdir>/run_manifest.json` is the authoritative run-level diagnostic artifact.

Conceptual schema:

```json
{
  "schema_version": 1,
  "run_id": "uuid",
  "status": "COMPLETED",
  "created_at": "2026-08-31T18:00:00Z",
  "updated_at": "2026-08-31T18:05:00Z",
  "completed_at": "2026-08-31T18:05:00Z",
  "target_name": "target",
  "effective_config": {},
  "environment": {},
  "input_provenance": {},
  "reference": {
    "id": "NC_000000.1",
    "mode": "explicit"
  },
  "stage_actions": [
    {"stage": "input", "action": "RUN", "outcome": "COMPLETE"}
  ],
  "scientific_completeness": {
    "complete": true,
    "missing_evidence": []
  },
  "attempts": [],
  "failure": null
}
```

### Effective configuration

The manifest stores the validated effective configuration, including defaults, in a deterministic JSON-compatible representation. Paths are serialized as strings. Values are taken from `PipelineConfig`, not copied from raw YAML, so the manifest describes the configuration the program actually used.

Environment variables and secrets are never part of the effective configuration. NCBI API keys and e-mail values are not copied into the manifest.

### Stage actions

The run-level stage list records only:

- stage name;
- planned/actual action;
- execution outcome;
- checkpoint manifest path relative to `<outdir>` when available.

It does not duplicate checkpoint parameters, dependency fingerprints, tool fingerprints, output hashes or typed stage state. Those remain authoritative in `.checkpoints/<stage>/manifest.json`.

### Timestamps and fingerprints

Run and attempt timestamps are UTC ISO-8601 values used only for audit. They never participate in scientific or checkpoint fingerprints.

The run manifest itself is diagnostic and is not a checkpoint-critical stage output. Updating timestamps or retry history must not invalidate scientific checkpoints.

## Manifest lifecycle and atomicity

Before the recorder starts, the pipeline performs output-directory safety checks, including the frozen NCBI dataset guard. A rejected unsafe destination is not considered a started run because Geison must not write there.

After safety validation:

1. the output directory is created if needed;
2. `RunRecorder` creates or loads the run identity;
3. a new attempt with status `RUNNING` is published atomically;
4. the environment and execution plan are recorded;
5. stage outcomes are published after each stage action;
6. final provenance and completeness are recorded;
7. the attempt and top-level run become `COMPLETED`, `PARTIAL` or `FAILED`.

Manifest writes use temporary files plus atomic replacement in the same directory.

If an exception occurs after the recorder starts, Geison makes a best effort to publish `FAILED`, records the failed stage or `null` when failure preceded the first stage, and re-raises the original exception so the CLI exits non-zero. A diagnostic publication error must not hide the original scientific error.

A killed process may leave a `RUNNING` attempt. On the next invocation, the recorder marks the previous unfinished attempt `FAILED` with code `INTERRUPTED` before appending the new attempt.

## Input provenance

Input provenance is metadata only; no sequence bodies are copied into diagnostics.

### Local input

The manifest records:

- input kind (`fasta` or `genbank`);
- configured source path;
- source SHA-256 already used by checkpoint identity when available;
- accepted/rejected counts after QC.

### Frozen NCBI dataset

The manifest records:

- mode `frozen_dataset`;
- configured dataset path;
- dataset manifest identity;
- original query or requested accessions when present in the frozen manifest;
- resolved accession versions;
- acquisition metadata needed to identify the materialized dataset.

### Online NCBI acquisition

The manifest records:

- mode `query` or `accessions`;
- the exact query or requested accession list;
- resolved accessions and accession versions from the existing acquisition result/manifest;
- materialized dataset identity.

The NCBI API key, request headers, raw payloads and full sequences are forbidden.

## Reference provenance

After alignment, the manifest records the selected `reference_id` and `reference_mode` from the alignment result. If downstream stages are skipped or the run fails before a reference exists, both are `null`.

The run manifest references the selected identity only. Alignment coordinates and annotations remain in existing scientific artifacts and checkpoints.

## Scientific completeness and final status

Final run status is computed by one explicit completeness evaluator after all stages return successfully.

`COMPLETED` requires all of the following:

- Evaluation Set is non-empty;
- primer design produced at least one assay for final evaluation;
- inclusivity completed with the required evidence;
- specificity completed with the required evidence;
- ranking completed;
- no completeness invariant is violated.

If execution returns without exception but one or more requirements are absent, status is `PARTIAL`. Stable missing-evidence codes include:

- `EMPTY_EVALUATION_SET`;
- `NO_ASSAYS`;
- `INCLUSIVITY_NOT_COMPLETE`;
- `SPECIFICITY_NOT_COMPLETE`;
- `RANKING_NOT_COMPLETE`.

Disabled stages produce the corresponding missing-evidence code rather than being treated as errors.

Any exception after run start produces `FAILED` regardless of the amount of valid upstream work. Valid completed stage checkpoints remain reusable.

The current ranking rules already convert missing critical evidence into `REVIEW` with `score_status = INCOMPLETE`. Issue #12 preserves those rules and supplies ranking with a typed pre-ranking completeness context covering Evaluation Set, assay availability, inclusivity and specificity. Ranking must deny `IN SILICO PASS` whenever that context contains a missing-evidence code.

Before ranking artifacts are atomically published, an independent safety invariant checks that no partial-evidence result contains `IN SILICO PASS`. A violation raises `ScientificCompletenessError`, makes the run `FAILED` and prevents the inconsistent ranking artifacts from being promoted. This guard is defensive; normal missing-evidence behavior remains `REVIEW`/`INCOMPLETE`, not an exception.

`RunSummary.status` and `run_summary.json` use the same final status as `run_manifest.json`.

## Structured event log

`<outdir>/run.log.jsonl` is append-only across attempts. Every complete line is one JSON object with an explicit schema version.

Event types are intentionally small:

- `run_started`;
- `environment_inspected`;
- `plan_created`;
- `stage_started`;
- `stage_reused`;
- `stage_completed`;
- `stage_failed`;
- `run_completed`;
- `run_failed`.

Common fields:

```json
{
  "schema_version": 1,
  "timestamp": "2026-08-31T18:00:00Z",
  "level": "INFO",
  "event": "stage_completed",
  "run_id": "uuid",
  "attempt_id": "uuid",
  "stage": "alignment",
  "action": "RUN",
  "status": "COMPLETE",
  "message": "Alignment stage completed"
}
```

Events use per-type field allowlists. The logger never accepts arbitrary scientific result objects, raw configuration dictionaries or external response payloads.

Each line is serialized as compact JSON and flushed after append. The manifest remains authoritative if a process termination leaves an incomplete final log line.

## Sanitization policy

All user-visible failure details pass through one sanitizer before entering the manifest or structured log.

The sanitizer:

- redacts values whose normalized field names match sensitive segments such as `token`, `api_key`, `secret`, `password`, `authorization` or `email`;
- redacts long nucleotide-like strings rather than logging sequence content;
- limits message and field lengths;
- preserves exception type, stable diagnostic code and a concise actionable message;
- never serializes traceback locals or arbitrary object representations.

Normal progress events are constructed from trusted templates and stable identifiers. Sanitization is defense in depth, not permission to pass arbitrary payloads to the logger.

Console errors may remain more detailed for local debugging, but secrets and full sequences are also forbidden there when the error is produced by new issue #12 code.

## `doctor` command

The CLI gains:

```text
qpcr-pipeline doctor
```

`doctor` requires no configuration and performs no network calls. It reports:

- Python and Geison versions;
- Git availability and commit when inside a checkout;
- CD-HIT-EST, MAFFT and Primer3 availability/version;
- BLAST+ as `NOT_USED`, while optionally reporting whether it is installed.

Because no configuration is supplied, external scientific tools are diagnostic/optional in this command. Missing external tools produce clear warnings but exit code `0`. Exit code is non-zero only when Python/Geison inspection itself cannot produce a valid diagnostic report or an internal inspection failure prevents the command from completing.

The command renders a concise human-readable table. Its underlying result is a typed serializable object used directly by tests and future interfaces; issue #12 does not add a separate `--json` CLI flag.

## `--dry-run`

The run command gains:

```text
qpcr-pipeline run CONFIG --dry-run [--outdir OUTDIR]
```

Dry-run:

1. loads and validates the effective configuration;
2. validates execution-policy flag combinations;
3. computes the same execution plan used by a real run;
4. inspects the environment and tools required by that configuration;
5. prints target, stage actions, reasons and missing required tools;
6. exits without starting a scientific run.

If no output directory and no resume control are supplied, every stage is described as `RUN`. Resume, from-step and force-step retain the existing requirement for `--outdir`, because checkpoint inspection needs a concrete run directory.

Dry-run must not:

- create or modify the output directory;
- create a run manifest or log;
- create or replace checkpoint manifests;
- query NCBI;
- run CD-HIT, MAFFT, Primer3 or other scientific binaries;
- publish scientific artifacts.

Version probes used by `EnvironmentInspector` are allowed because they are diagnostics, not scientific execution.

Missing tools required by the effective configuration are reported and produce a non-zero exit code. Missing tools for disabled stages and BLAST+ do not fail dry-run.

## CLI failure behavior

- Configuration and argument errors remain concise and non-zero.
- A real run failure after recorder initialization publishes a `FAILED` manifest and `run_failed` event, then returns non-zero through the existing exception path.
- A dry-run failure never creates diagnostic or scientific artifacts.
- `doctor` never changes repository or run state.

## Compatibility with issue #11

- Stage checkpoint layout and schemas remain unchanged unless a narrowly required planner extraction needs an internal API adjustment.
- Checkpoint fingerprints continue to exclude timestamps, run IDs, attempt IDs and run-level status.
- `--resume`, `--from-step` and `--force-step` retain their existing semantics.
- Completed upstream checkpoints survive a failed run and remain eligible for later resume.
- The run manifest records checkpoint paths and actions but never becomes a dependency of stage fingerprints.
- Existing `run_summary.json` and `qc_report.json` remain user-facing outputs; issue #12 extends status consistency without removing existing fields.

## Testing strategy

Implementation follows TDD. Every behavior is first demonstrated by a failing test.

### Environment inspector tests

- Python and Geison versions are reported;
- Git commit available, unavailable and dirty states;
- present and absent CD-HIT, MAFFT and Primer3;
- BLAST+ always remains `NOT_USED` and optional;
- subprocess output normalization and failure isolation.

### Execution planner and dry-run tests

- normal dry-run plans all stages as `RUN`;
- resume reports valid `REUSE` and invalid/dependent `RUN/FORCED` actions;
- from-step validates reusable boundaries without writes;
- force-step matches real execution semantics;
- dry-run does not create an absent output directory;
- dry-run does not alter an existing output directory;
- dry-run does not call NCBI or scientific runners;
- missing required tools fail dry-run, while unused tools do not;
- invalid flag combinations remain rejected.

### Run recorder tests

- first attempt creates a stable run identity;
- resume retains the run identity and appends a new attempt;
- an unfinished prior attempt becomes `FAILED`/`INTERRUPTED`;
- manifest writes are atomic;
- failure before the first stage records `stage: null`;
- failure during a stage records the correct stage and sanitized error;
- recorder publication failure does not replace the original exception;
- timestamps do not affect checkpoint fingerprints.

### Completeness tests

- complete evidence produces `COMPLETED`;
- disabled ranking produces `PARTIAL`;
- disabled inclusivity produces `PARTIAL`;
- disabled specificity produces `PARTIAL`;
- empty Evaluation Set produces `PARTIAL`;
- each missing-evidence code is deterministic;
- partial execution cannot contain `IN SILICO PASS`;
- a scientific stage exception produces `FAILED` while preserving completed checkpoints.

### Provenance tests

- effective configuration includes validated defaults;
- local input kind/path/hash and QC counts are recorded without sequences;
- NCBI query, requested accessions, resolved accession versions and frozen dataset identity are preserved;
- explicit and automatically selected reference modes are preserved;
- missing Git metadata produces `null` fields without failing the run.

### Logging and sanitization tests

- expected lifecycle events are emitted in order;
- reused and forced stages are distinguishable;
- API-like secrets, e-mail fields and long nucleotide strings are redacted;
- arbitrary scientific result objects cannot be logged;
- messages are bounded in size;
- log append behavior preserves prior attempts.

### Regression and integration tests

- the complete existing unit suite remains green;
- interruption and resume behavior from issue #11 remains green;
- a fixture run produces manifest, event log, run summary and scientific artifacts consistently;
- CircleCI on `main` continues to exercise real CD-HIT, MAFFT and Primer3 integration tests.

## Acceptance mapping

| Issue #12 criterion | Design coverage |
| --- | --- |
| Manifest with effective configuration, Geison/Git/Python/tool versions | Run manifest + environment inspection |
| Query, accessions/versioning and selected reference preserved | Input and reference provenance |
| Explicit COMPLETED, PARTIAL or FAILED | Completeness evaluator + lifecycle |
| Incomplete runs cannot produce IN SILICO PASS | Ranking preservation + final safety invariant |
| `doctor` verifies dependencies/environment | `doctor` command contract |
| `--dry-run` validates and describes without execution | Shared execution plan + side-effect prohibitions |
| Structured logs without secrets/full sequences | Event allowlists + sanitization policy |

## Expected implementation boundaries

Expected new focused modules:

- `qpcr_pipeline/diagnostics.py` for environment inspection and tool reports;
- `qpcr_pipeline/run_recording.py` for manifest models, lifecycle, events, sanitization and completeness evaluation.

Expected targeted modifications:

- `qpcr_pipeline/execution.py` for the shared read-only execution plan;
- `qpcr_pipeline/pipeline.py` for recorder lifecycle and plan execution;
- `qpcr_pipeline/cli.py` for `doctor` and `--dry-run`;
- existing config/ranking/NCBI integration points only where typed provenance or invariant enforcement requires them;
- focused new tests plus minimal updates to existing expectations.

No unrelated scientific algorithm changes, dependency additions or large-scale module refactors belong in issue #12.

