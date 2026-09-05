# Guided NCBI Panel Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the default Guided Colab path fetch target and challenge sequence datasets from NCBI after panel approval, without manual FASTA upload/path management.

**Architecture:** Reuse the existing target `NcbiInputConfig`/`acquire_ncbi_dataset` path. Extend off-target configuration with live NCBI sources and add a checkpointed `challenge_acquisition` stage that materializes every challenge into run-local frozen evidence before contrast/specificity. Add a deterministic guided preset builder (WNV first) and wire the notebook to it, while keeping synthetic and advanced-local modes.

**Tech Stack:** Python 3.10+, dataclasses, pathlib, PyYAML, Biopython Entrez/SeqIO, existing Geison checkpoint framework, pytest, Google Colab notebook JSON.

**Spec:** `docs/superpowers/specs/2026-09-05-guided-ncbi-panel-flow-design.md`

## Global Constraints

- Do not duplicate NCBI HTTP logic in the notebook; all live acquisition goes through `qpcr_pipeline.ncbi`.
- Panel approval remains mandatory before live scientific acquisition.
- NCBI credentials/contact data stay in environment only and never enter durable artifacts.
- Guided WNV limits are operational smoke-test limits, not claims of representative population sampling.
- Existing local FASTA and frozen off-target inputs remain backward compatible.
- Contrastive conservation and final specificity must consume the exact same materialized challenge datasets.
- No report visual-theme redesign in this plan.

---

### Task 1: Extend off-target configuration with live NCBI input

**Files:**
- Modify: `qpcr_pipeline/config/legacy.py`
- Test: `tests/test_specificity_config.py`
- Test: `tests/test_config.py`

**Interfaces:**
- Consumes: existing `NcbiInputConfig`, `_parse_ncbi_input()`, `validate_ncbi_input_config()`.
- Produces: `OffTargetConfig.ncbi: NcbiInputConfig | None`; off-target YAML accepts exactly one of `fasta`, `frozen_dataset`, `ncbi`.

- [ ] **Step 1: Write failing parsing/validation tests**

Add tests equivalent to:

```python
def test_off_target_accepts_live_ncbi_query(tmp_path):
    path = tmp_path / "config.yaml"
    path.write_text("""
target: {name: Target}
input: {fasta: target.fasta}
off_targets:
  - name: Neighbor
    ncbi:
      query: '\"Neighbor virus\"[Organism] AND complete genome[Title]'
      batch_size: 20
      retries: 5
      max_records: 20
""")
    config = load_config(path)
    assert config.off_targets[0].ncbi is not None
    assert config.off_targets[0].ncbi.max_records == 20
```

Also assert that `fasta + ncbi`, `frozen_dataset + ncbi`, and `off_targets[].ncbi.frozen_dataset` are rejected.

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
python -m pytest tests/test_specificity_config.py tests/test_config.py -q
```

Expected: failures because `OffTargetConfig` has no `ncbi` field and parser rejects `ncbi`.

- [ ] **Step 3: Implement config support**

Change the model to:

```python
@dataclass(frozen=True, slots=True)
class OffTargetConfig:
    name: str
    fasta: Path | None = None
    frozen_dataset: Path | None = None
    ncbi: NcbiInputConfig | None = None
```

Update `_parse_off_targets()` allowed fields and parse `ncbi` through `_parse_ncbi_input`. Update `validate_off_target_config()` so exactly one source is configured. When `ncbi` is set, validate it and reject `ncbi.frozen_dataset` with an actionable message directing users to top-level `frozen_dataset`.

- [ ] **Step 4: Run focused tests and verify GREEN**

```bash
python -m pytest tests/test_specificity_config.py tests/test_config.py -q
```

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/config/legacy.py tests/test_specificity_config.py tests/test_config.py
git commit -m "feat: allow NCBI challenge sources"
```

---

### Task 2: Add checkpointed challenge acquisition stage

**Files:**
- Modify: `qpcr_pipeline/off_targets.py`
- Modify: `qpcr_pipeline/execution.py`
- Modify: `qpcr_pipeline/checkpoint_codecs.py`
- Modify: `qpcr_pipeline/checkpoint_stages.py`
- Modify: `qpcr_pipeline/pipeline.py`
- Test: `tests/test_off_target_acquisition.py` (create)
- Test: `tests/test_checkpoint_codecs.py`
- Test: `tests/test_checkpoint_stages.py`
- Test: `tests/test_execution.py`
- Test: `tests/test_execution_plan.py`

**Interfaces:**
- Produces `PreparedChallengeDataset`, `ChallengeAcquisitionResult`, `materialize_challenge_datasets()`, `challenge_configs_from_result()`.
- Produces stage `challenge_acquisition`.
- Downstream `contrastive_conservation` and `specificity` consume `challenge_configs_from_result(results["challenge_acquisition"])`.

- [ ] **Step 1: Write failing materialization tests**

Cover:

```python
def test_live_ncbi_challenge_is_materialized_with_fake_client(...):
    result = materialize_challenge_datasets(
        (OffTargetConfig(name="Neighbor", ncbi=NcbiInputConfig(accessions=("NC_1.1",))),),
        outdir,
        ncbi_client=fake_client,
    )
    dataset = result.datasets[0]
    assert dataset.source_type == "NCBI_FROZEN"
    assert dataset.source_path == outdir / "challenge_datasets/001-neighbor"
    assert dataset.sequence_count == 1
    assert dataset.manifest_sha256
```

Also test local FASTA is copied under `outdir`, external frozen datasets are copied/validated, names are slugged deterministically, and all returned paths are under `outdir`.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_off_target_acquisition.py -q
```

Expected: import failures for new models/functions.

- [ ] **Step 3: Implement materialization models/functions**

Add dataclasses:

```python
@dataclass(frozen=True, slots=True)
class PreparedChallengeDataset:
    name: str
    source_type: Literal["FASTA", "NCBI_FROZEN"]
    source_path: Path
    sequence_count: int
    records_sha256: str
    manifest_sha256: str | None = None

@dataclass(frozen=True, slots=True)
class ChallengeAcquisitionResult:
    status: Literal["COMPLETE"]
    datasets: tuple[PreparedChallengeDataset, ...]
```

Implement deterministic run-local materialization under `challenge_datasets/NNN-slug` and conversion back to ordinary `OffTargetConfig` values.

- [ ] **Step 4: Add stage/checkpoint tests first**

Update execution tests to expect:

```python
STAGE_ORDER == (
    "panel", "input", "qc", "challenge_acquisition", "clustering",
    "alignment", "conservation", "contrastive_conservation",
    "primer_design", "inclusivity", "specificity", "ranking",
)
```

Dependencies:

```python
"challenge_acquisition": ("panel",),
"contrastive_conservation": ("panel", "conservation", "challenge_acquisition"),
"specificity": ("primer_design", "challenge_acquisition"),
```

Add codec round-trip and stage-parameter/input-identity tests.

- [ ] **Step 5: Verify RED for stage tests**

```bash
python -m pytest tests/test_execution.py tests/test_execution_plan.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py -q
```

- [ ] **Step 6: Implement the stage**

Add `CHALLENGE_ACQUISITION_CODEC`, stage definitions/parameters/input identities/outputs, `_run_stage("challenge_acquisition")`, and downstream use of materialized configs.

For `stage_input_identities("challenge_acquisition")`:

- local FASTA -> source hash;
- external frozen -> records + manifest hashes;
- live NCBI -> no file identity; query/accession settings are in stage parameters.

- [ ] **Step 7: Run all Task 2 tests**

```bash
python -m pytest tests/test_off_target_acquisition.py tests/test_execution.py tests/test_execution_plan.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py -q
```

- [ ] **Step 8: Commit**

```bash
git add qpcr_pipeline/off_targets.py qpcr_pipeline/execution.py qpcr_pipeline/checkpoint_codecs.py qpcr_pipeline/checkpoint_stages.py qpcr_pipeline/pipeline.py tests/test_off_target_acquisition.py tests/test_execution.py tests/test_execution_plan.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py
git commit -m "feat: checkpoint challenge dataset acquisition"
```

---

### Task 3: Prove target and challenge NCBI acquisition in one pipeline run

**Files:**
- Create: `tests/test_pipeline_ncbi_challenges.py`
- Modify if required: `tests/pipeline_checkpoint_fixtures.py`

**Interfaces:**
- Consumes existing `NcbiClient` injection on `run_pipeline` for both target and challenge acquisition.
- Proves challenge stage feeds contrast and specificity using one frozen result.

- [ ] **Step 1: Write failing integration test with deterministic fake NCBI client**

Construct an approved panel and config with:

```python
input_ncbi=NcbiInputConfig(accessions=("TARGET.1", "TARGET.2")),
off_targets=(
    OffTargetConfig(name="Neighbor", ncbi=NcbiInputConfig(accessions=("NEIGHBOR.1",))),
),
```

Use the existing fake-runner pattern for MAFFT/Primer3 where needed. Assert:

- target `ncbi_dataset/dataset_manifest.json` exists;
- challenge `challenge_datasets/001-neighbor/dataset_manifest.json` exists;
- run reaches at least the normal scientific stages rather than failing on a missing local path;
- contrastive and specificity checkpoint dependencies include the challenge-acquisition fingerprint;
- only one materialized challenge dataset is used by both stages.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_pipeline_ncbi_challenges.py -q
```

- [ ] **Step 3: Make minimal integration corrections**

Only change production code where the new integration test exposes a real contract mismatch.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_pipeline_ncbi_challenges.py -q
```

- [ ] **Step 5: Commit**

```bash
git add tests/test_pipeline_ncbi_challenges.py qpcr_pipeline
git commit -m "test: cover NCBI target and challenge pipeline"
```

---

### Task 4: Add deterministic guided WNV panel/config builder

**Files:**
- Create: `qpcr_pipeline/guided.py`
- Create: `tests/test_guided_config.py`

**Interfaces:**
- Produces `supported_guided_targets() -> tuple[str, ...]`.
- Produces `build_guided_configs(target_name: str) -> tuple[dict[str, object], dict[str, object]]`.

- [ ] **Step 1: Write failing WNV config tests**

Assert the proposal config contains:

```python
proposal["target"]["name"] == "West Nile virus"
proposal["input"]["ncbi"]["query"] == '"West Nile virus"[Organism] AND complete genome[Title]'
proposal["contrastive_conservation"]["enabled"] is False
```

Assert challenge names/criticalities are Usutu CRITICAL, JEV CRITICAL, Dengue IMPORTANT and every off-target uses `ncbi`, not `fasta`.

Assert approved template switches to frozen panel placeholder and enables contrastive conservation. Assert unsupported target errors name the supported target(s).

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_guided_config.py -q
```

- [ ] **Step 3: Implement `qpcr_pipeline.guided`**

Use immutable dataclasses/constants with `knowledge_version = "2026-09-05"`. Keep all proposal generation deterministic and network-free.

- [ ] **Step 4: Verify GREEN**

```bash
python -m pytest tests/test_guided_config.py -q
```

- [ ] **Step 5: Commit**

```bash
git add qpcr_pipeline/guided.py tests/test_guided_config.py
git commit -m "feat: add guided West Nile NCBI preset"
```

---

### Task 5: Rewire Guided Colab around automatic NCBI acquisition

**Files:**
- Modify: `notebooks/geison_guided_colab.ipynb`
- Modify: `tests/test_guided_colab_notebook.py`
- Modify: `docs/guided-colab.md`

**Interfaces:**
- Default mode becomes `Guided (NCBI)`.
- Direct helper imports work in active Colab kernel through explicit `/content/Geison` insertion.
- Guided mode writes configs from `build_guided_configs()` and exports `NCBI_EMAIL`.
- Advanced local mode preserves existing manual FASTA fields.

- [ ] **Step 1: Strengthen notebook contract tests first**

Require notebook code to contain all of:

```python
'Guided (NCBI)'
'build_guided_configs'
'os.environ["NCBI_EMAIL"]'
'sys.path.insert(0, str(geison_repo))'
```

Require Guided mode not to construct `input: {fasta: project_target_fasta}`. Preserve assertions for panel approval, contrast view, report, evidence bundle and CLI execution.

- [ ] **Step 2: Verify RED**

```bash
python -m pytest tests/test_guided_colab_notebook.py -q
```

- [ ] **Step 3: Modify notebook**

Setup/kernel cell:

```python
from pathlib import Path
import os, sys
geison_repo = Path("/content/Geison")
if str(geison_repo) not in sys.path:
    sys.path.insert(0, str(geison_repo))
```

Guided form:

```python
mode = "Guided (NCBI)"  #@param ["Guided (NCBI)", "Demo (synthetic)", "Advanced (local sequences)"]
target_name = "West Nile virus"  #@param {type:"string"}
ncbi_email = ""  #@param {type:"string"}
```

Before approved live run, require nonblank email and set `os.environ["NCBI_EMAIL"]`.

Generate proposal/approved template using `build_guided_configs(target_name)`. Keep synthetic and advanced paths intact.

- [ ] **Step 4: Update docs**

Document Guided NCBI as the recommended workflow, NCBI email requirement, panel approval, automatic challenge acquisition, and Advanced local mode.

- [ ] **Step 5: Verify GREEN**

```bash
python -m pytest tests/test_guided_colab_notebook.py -q
```

- [ ] **Step 6: Commit**

```bash
git add notebooks/geison_guided_colab.ipynb tests/test_guided_colab_notebook.py docs/guided-colab.md
git commit -m "feat: make guided Colab NCBI-first"
```

---

### Task 6: Regression and completion verification

**Files:**
- Modify only if a regression exposes a real issue.

- [ ] **Step 1: Run focused feature suite**

```bash
python -m pytest tests/test_config.py tests/test_specificity_config.py tests/test_off_target_acquisition.py tests/test_checkpoint_codecs.py tests/test_checkpoint_stages.py tests/test_execution.py tests/test_execution_plan.py tests/test_pipeline_ncbi_challenges.py tests/test_guided_config.py tests/test_guided_colab_notebook.py -q
```

Expected: all pass.

- [ ] **Step 2: Run normal offline suite**

```bash
python -m pytest tests -q
```

Expected: all pass, no live network required.

- [ ] **Step 3: Run existing real-tool integration suite when available**

```bash
python -m pytest integration_tests -q
```

Expected: all configured integration tests pass in the CI environment that has MAFFT/Primer3/CD-HIT.

- [ ] **Step 4: Inspect branch diff**

Confirm no report-theme redesign, no credentials, no unrelated branch-config changes, and no synthetic fixture changes intended to manufacture a PASS.

- [ ] **Step 5: Commit any verification-only corrections**

Use a narrow message describing the actual regression fix; do not make a blanket cleanup commit.