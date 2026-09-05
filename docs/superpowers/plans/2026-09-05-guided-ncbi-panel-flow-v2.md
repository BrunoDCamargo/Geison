# Guided NCBI Panel Workflow V2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the default Guided Colab path run a real West Nile workflow from target selection to report without manual FASTA downloads/uploads.

**Architecture:** Keep the scientific pipeline unchanged. Use the existing target `input.ncbi` support directly. Add a `qpcr_pipeline.guided` orchestration layer and CLI commands that, after panel approval, acquire challenge datasets using the existing `acquire_ncbi_dataset()` service and write a standard approved config using `off_targets[].frozen_dataset`.

**Tech Stack:** Python 3.10+, PyYAML, existing Biopython/NCBI acquisition service, argparse CLI, pytest, Google Colab notebook JSON.

**Spec:** `docs/superpowers/specs/2026-09-05-guided-ncbi-panel-flow-design.md`

## Global Constraints

- No NCBI HTTP code in the notebook.
- Panel approval occurs before challenge acquisition.
- No changes to contrast, specificity, ranking, or scientific thresholds.
- No credentials persisted in artifacts.
- Keep Demo and Advanced local modes.
- Guided WNV sequence caps are smoke/workbench limits, not representative-sampling claims.

---

### Task 1: Guided knowledge/config module

**Files:**
- Create: `qpcr_pipeline/guided.py`
- Create: `tests/test_guided_config.py`

**Produces:**
- `supported_guided_targets()`
- `build_guided_proposal_config(target_name)`
- `finalize_guided_project(target_name, approved_panel_path, workspace, ncbi_client=None)`

- [ ] Write tests first for WNV proposal: target uses `input.ncbi`; panel proposes Usutu CRITICAL, JEV CRITICAL, Dengue IMPORTANT; no local FASTA fields; unsupported targets fail explicitly.
- [ ] Run focused test and observe RED due to missing module.
- [ ] Implement deterministic WNV preset with knowledge version `2026-09-05`, target max_records 50, challenge max_records 20, batch_size 20, retries 5.
- [ ] Add finalizer tests with a deterministic fake `NcbiClient`. After an approved panel is supplied, each challenge is acquired under `<workspace>/guided_challenges/NNN-slug`, `config-approved.yaml` uses those directories as `frozen_dataset`, target remains `input.ncbi`, contrastive mode is enabled, and `guided_acquisition_manifest.json` contains hashes/counts but no email/API key.
- [ ] Implement finalizer by reusing `acquire_ncbi_dataset()` and `load_approved_panel_manifest()`; fail if the approved target or challenge names do not match the supported preset.
- [ ] Run focused tests GREEN.

---

### Task 2: CLI commands

**Files:**
- Modify: `qpcr_pipeline/cli.py`
- Create: `tests/test_guided_cli.py`

**Produces:**
- `qpcr-pipeline guided prepare --target ... --workspace ...`
- `qpcr-pipeline guided finalize --target ... --workspace ... --approved-panel ...`

- [ ] Write CLI tests first. `prepare` must write `config-proposal.yaml`; `finalize` must call the guided finalizer and report the approved config path. Tests inject/patch the guided boundary, not NCBI transport internals.
- [ ] Verify RED because `guided` parser/action does not exist.
- [ ] Implement nested `guided` subparser without changing existing `run`, `doctor`, or `panel approve` behavior.
- [ ] Run CLI tests GREEN plus existing `tests/test_cli.py`.

---

### Task 3: Guided Colab contract

**Files:**
- Modify: `notebooks/geison_guided_colab.ipynb`
- Modify: `tests/test_guided_colab_notebook.py`
- Modify: `docs/guided-colab.md`

- [ ] Strengthen notebook test first to require `Guided (NCBI)` default, `qpcr-pipeline guided prepare`, `qpcr-pipeline guided finalize`, `NCBI_EMAIL`, and explicit `sys.path.insert(0, str(geison_repo))`.
- [ ] Verify RED against current notebook.
- [ ] Update notebook modes to `Guided (NCBI)`, `Demo (synthetic)`, `Advanced (local sequences)`.
- [ ] Guided mode accepts target, NCBI email, workspace only; it calls `guided prepare`, then existing run for panel gate; after `APROVAR`, it approves the panel, calls `guided finalize`, displays `config-approved.yaml`, then resumes the normal pipeline.
- [ ] Keep manual FASTA/challenge fields only inside Advanced mode.
- [ ] Add kernel `sys.path` insertion before any direct package helper import so the report/evidence bundle section works in the same Colab kernel.
- [ ] Update docs and run notebook tests GREEN.

---

### Task 4: Offline regression and live acceptance handoff

**Files:**
- Modify only for real regressions.

- [ ] Run focused suite: `tests/test_guided_config.py tests/test_guided_cli.py tests/test_guided_colab_notebook.py tests/test_ncbi_acquisition.py tests/test_cli.py`.
- [ ] Run `python -m pytest tests -q`.
- [ ] Run integration suite in CI environment.
- [ ] Restore any temporary CI branch filter before merge.
- [ ] Inspect diff for credentials, unrelated UI/theme changes, or scientific-threshold changes.
- [ ] After automated verification, ask the user only to reset Colab and run the new Guided (NCBI) flow from the top with West Nile virus; no FASTA uploads.