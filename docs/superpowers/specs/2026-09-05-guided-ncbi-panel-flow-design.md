# Guided NCBI Panel Workflow Design

## Goal

Make Geison's default Guided Colab workflow match the intended product experience: the researcher chooses a supported biological target, Geison acquires the target and challenge sequence datasets from NCBI, presents a scientific panel for explicit approval, then runs the existing assay-discovery pipeline without requiring manual FASTA download/upload or path management.

The existing local FASTA/GenBank workflow remains available as an advanced path.

## Why this change

Geison already has a reproducible NCBI acquisition subsystem for target inputs (`NcbiInputConfig`, `BioEntrezClient`, `acquire_ncbi_dataset`, frozen manifests, retries and resume). The guided notebook currently bypasses it and makes `Project` mode depend on local FASTA files. Off-targets/challenges currently accept only local FASTA or an already-frozen NCBI dataset. This design connects the guided experience to the existing target acquisition capability and extends the same acquisition model to challenge datasets.

This work implements the previously deferred “panel intelligence and validation” direction from the hybrid assay architecture, beginning with a small versioned arbovirus knowledge seed rather than pretending to provide a universal clinical knowledge base.

## User experience

The notebook exposes three modes:

1. **Guided (NCBI)** — default. User selects/types a supported target and provides the NCBI contact email required by the E-utilities integration. No FASTA paths or uploads are requested.
2. **Demo (synthetic)** — deterministic synthetic walkthrough retained for regression/demo use.
3. **Advanced (local sequences)** — existing bring-your-own FASTA path workflow.

For Guided (NCBI):

```text
Target selection
      |
      v
Geison builds proposal from versioned panel knowledge
      |
      v
PANEL_APPROVAL_REQUIRED
      |
      v
User reviews target + challenge organisms + criticality
      |
      v
APROVAR
      |
      v
Target NCBI acquisition + challenge NCBI acquisition
      |
      v
QC -> clustering/alignment -> conservation -> contrast
      |
      v
primer/probe design -> inclusivity -> specificity -> ranking -> report
```

The only scientific human gate remains panel approval. Network identity configuration is not treated as scientific input.

## Supported guided knowledge

The architecture is generic, but the first knowledge seed is intentionally small and explicit. It supports **West Nile virus** first because the current acceptance exercise uses WNV and the repository already contains a WNV panel fixture.

Initial WNV proposal:

- Target: West Nile virus, `broad_detection`.
- CRITICAL challenge: Usutu virus — reason `phylogenetic_neighbor`.
- CRITICAL challenge: Japanese encephalitis virus — reason `phylogenetic_neighbor`.
- IMPORTANT challenge: Dengue virus — reason `clinical_differential`.

The preset also records a `knowledge_version` and `proposed_by: geison_guided_knowledge` so the panel proposal is auditable. Unsupported targets fail clearly with the supported target list instead of silently fabricating a panel.

This is not an authoritative clinical panel. The user must approve the proposal and the report continues to state that in-silico evidence requires scientific and experimental validation.

## NCBI queries

Each guided dataset uses the existing `NcbiInputConfig` model. Queries are explicit strings and have explicit bounded `max_records` values so a Colab run cannot accidentally attempt an unbounded public-database download.

Initial smoke-test defaults:

- WNV target: `"West Nile virus"[Organism] AND complete genome[Title]`, `max_records: 50`.
- Each WNV challenge: organism-specific `complete genome[Title]`, `max_records: 20`.
- `batch_size: 20`.
- `retries: 5`.

These limits are operational defaults for the guided smoke workflow, not claims that the first N records are a scientifically representative population. The generated proposal/report must identify this limitation. Future representative-selection work can replace the bounded-prefix policy without changing the acquisition contract.

## NCBI policy and contact identity

The current production adapter requires `NCBI_EMAIL`; this remains required. The Guided notebook collects it once in setup and sets `NCBI_EMAIL` in the runtime before any approved live acquisition. `NCBI_API_KEY` remains optional and may be supplied through the environment.

No credentials are persisted in configuration, panel manifests, checkpoint state, reports, logs, or evidence bundles.

## Configuration model

`OffTargetConfig` gains an optional `ncbi: NcbiInputConfig` source.

Exactly one source is allowed per off-target:

- `fasta`
- `frozen_dataset`
- `ncbi`

For off-target `ncbi`, query and accession modes are allowed. `ncbi.frozen_dataset` is rejected because the existing top-level `frozen_dataset` field is the canonical frozen representation for off-targets.

Example:

```yaml
off_targets:
  - name: Usutu virus
    ncbi:
      query: '"Usutu virus"[Organism] AND complete genome[Title]'
      batch_size: 20
      retries: 5
      max_records: 20
```

## Challenge acquisition stage

Add a real checkpointed pipeline stage named `challenge_acquisition` after `qc` and before `clustering` in stage order. It materializes every configured off-target into a deterministic run-local representation.

Why a stage instead of notebook-side downloads:

- network failures are recorded in normal run lifecycle;
- NCBI retry/freeze behavior remains inside Geison;
- challenge dataset composition participates in checkpoint identity;
- resume can reuse already-materialized challenge datasets;
- contrast and specificity consume the exact same frozen challenge evidence;
- the notebook does not become a second acquisition implementation.

### Materialization rules

For each off-target, in configuration order:

- local FASTA: copy to `<outdir>/challenge_datasets/<index>-<slug>/records.fasta`;
- external frozen NCBI dataset: validate it, then copy its `records.gb` and `dataset_manifest.json` into the run-local challenge directory;
- live NCBI query/accessions: call the existing `acquire_ncbi_dataset` directly in the run-local challenge directory.

The stage result records, per dataset:

- name;
- materialized source type (`FASTA` or `NCBI_FROZEN`);
- run-local source path;
- sequence count;
- records SHA-256;
- optional manifest SHA-256.

All stored paths are inside `outdir` so checkpoint codecs remain safe and portable.

### Downstream consumption

`contrastive_conservation` and `specificity` no longer read raw `config.off_targets` directly. They convert `ChallengeAcquisitionResult` into ordinary frozen/local `OffTargetConfig` values and pass those to the existing scientific functions unchanged.

Both stages depend on `challenge_acquisition`, so a changed challenge dataset fingerprint invalidates both scientific branches and downstream ranking.

## Checkpoint/provenance behavior

`challenge_acquisition` parameters include the declared source mode and NCBI query/accession settings but never credentials.

Input identities:

- local FASTA: hash of original file;
- external frozen NCBI: records + manifest hashes;
- live NCBI query/accession: no pre-existing file identity; the explicit query/accession parameters are checkpoint parameters.

Outputs include every run-local challenge `records.fasta` or `records.gb` plus every NCBI `dataset_manifest.json`.

The result fingerprint therefore freezes the exact challenge composition for downstream stages.

The top-level effective configuration retains the original user intent (`off_targets[].ncbi`) while the stage result and copied manifests retain the resolved dataset evidence.

## Guided configuration builder

Add a pure `qpcr_pipeline.guided` module that owns the versioned supported-target presets and creates two YAML-compatible mappings:

- proposal configuration with `panel.proposal` and `contrastive_conservation.enabled: false`;
- approved-template configuration with `panel.frozen_manifest: __APPROVED_PANEL_PATH__` and `contrastive_conservation.enabled: true`.

The builder has no network calls. It is deterministic and unit-testable.

The Guided notebook uses this builder rather than manually assembling the NCBI/panel dictionaries.

## Notebook behavior

The first Python cell explicitly adds `/content/Geison` to the active kernel `sys.path` before direct helper imports. This closes the currently observed Colab editable-install path-hook issue and also makes the final evidence-bundle import reliable in the same kernel.

In Guided (NCBI) mode:

- no target FASTA field;
- no challenge FASTA fields;
- no challenge-name editing fields in the default path;
- target is fed to `qpcr_pipeline.guided`;
- NCBI email is exported to `os.environ["NCBI_EMAIL"]`;
- the generated proposal is displayed before the existing approval gate;
- after approval, the existing CLI runs the full pipeline.

Advanced local mode preserves the current manual fields.

## Error behavior

- Missing NCBI email: fail before live acquisition with a clear action message.
- Unsupported guided target: fail before panel proposal with supported target names.
- NCBI HTTP 429/5xx/network errors: existing bounded retry/backoff behavior applies.
- Empty NCBI query result: acquisition fails explicitly; no empty dataset is silently accepted.
- Missing challenge data: the challenge acquisition stage fails; contrast/specificity do not run.
- Panel approval remains mandatory before any scientific execution or live sequence acquisition.

## Testing

Offline tests cover:

- off-target NCBI config parsing and mutual exclusion;
- challenge acquisition from deterministic fake NCBI client;
- local FASTA materialization;
- external frozen dataset materialization;
- challenge checkpoint codec round-trip;
- dependency invalidation for contrast/specificity;
- guided WNV proposal/config generation;
- notebook contract: Guided (NCBI) is default, no manual FASTA requirement in that mode, kernel source path is inserted, and NCBI email is exported before execution;
- full pipeline integration with fake NCBI target and challenge acquisition, proving both contrast and specificity consume the same frozen challenge stage output.

Live NCBI testing remains opt-in under `network_tests/` and must never be required by the normal test suite.

## Out of scope for this iteration

- universal automatic clinical-panel generation for arbitrary pathogens;
- literature/AI-based panel suggestions;
- metadata-stratified representative sampling;
- DECIPHER integration;
- report visual-theme/contrast redesign;
- claims of analytical or clinical validation.

These are intentionally separate from restoring the expected guided acquisition workflow.