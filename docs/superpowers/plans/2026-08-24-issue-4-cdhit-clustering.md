# Traceable CD-HIT-EST Clustering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add optional, reproducible CD-HIT-EST clustering that reduces the Discovery Set while preserving the complete Evaluation Set and publishing an auditable original-to-cluster-to-representative mapping.

**Architecture:** Add typed clustering configuration and an isolated `qpcr_pipeline.clustering` service behind an injectable runner protocol. The service uses stable internal FASTA IDs, strictly parses CD-HIT `.clstr` output, restores original IDs, and atomically publishes a Discovery Set FASTA and JSON report; the pipeline calls it after QC without changing Evaluation Set semantics.

**Tech Stack:** Python 3.10+, standard-library `subprocess`/`tempfile`/`shutil`/`json`/`pathlib`, Biopython `SeqIO`, PyYAML, `unittest`, CD-HIT-EST external binary for the optional integration test.

**Spec:** `docs/superpowers/specs/2026-08-24-issue-4-cdhit-clustering-design.md`

## Global Constraints

- Clustering defaults to disabled; existing runs must not require CD-HIT-EST.
- Enabled clustering uses only QC-approved records and never mutates, filters, or reorders the Evaluation Set.
- `identity` is configurable from 0.75 through 1.0; compatible CD-HIT word length is derived, never user-supplied.
- Production execution uses an argument list with `shell=False`, `-d 0`, `-g 1`, and `-r 0`.
- Standard tests remain offline and do not require CD-HIT-EST.
- The real-tool test lives under `integration_tests/` and skips when `cd-hit-est` is absent.
- Final artifacts use original sequence IDs and atomic sibling replacement.
- Duplicate approved sequence IDs and ambiguous/malformed CD-HIT compositions fail explicitly.
- No alignment, reverse-complement normalization, checkpoint/resume, or `doctor` behavior is added.

---

### Task 1: Parse and validate clustering configuration

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Produces: `ClusteringConfig(enabled: bool = False, identity: float = 0.95, threads: int = 1, memory_mb: int = 800)`
- Produces: `validate_clustering_config(config: ClusteringConfig) -> None`
- Extends: `PipelineConfig.clustering: ClusteringConfig`
- Consumes later: Task 2 clustering service and Task 4 pipeline integration

- [ ] **Step 1: Add RED configuration tests**

Add tests that prove omitted configuration uses defaults and YAML values are parsed:

```python
config = self._load_yaml(
    "target:\n  name: target\n"
    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
    "clustering:\n"
    "  enabled: true\n"
    "  identity: 0.9\n"
    "  threads: 4\n"
    "  memory_mb: 2048\n"
)
self.assertEqual(
    config.clustering,
    ClusteringConfig(enabled=True, identity=0.9, threads=4, memory_mb=2048),
)
```

Add table-driven invalid YAML cases for unknown keys, non-boolean `enabled`, boolean/non-number/out-of-range identity, boolean/non-integer/out-of-range threads, and non-positive memory. Add direct-construction cases and assert `config.selected_input` raises before any pipeline work:

```python
for clustering in (
    ClusteringConfig(enabled=1),
    ClusteringConfig(identity=True),
    ClusteringConfig(identity=0.74),
    ClusteringConfig(identity=1.01),
    ClusteringConfig(threads=0),
    ClusteringConfig(memory_mb=0),
):
    with self.subTest(clustering=clustering):
        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            clustering=clustering,
        )
        with self.assertRaises(ValueError):
            _ = config.selected_input
```

- [ ] **Step 2: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_config.PipelineConfigTests -v
```

Expected: import/attribute failures because `ClusteringConfig` and `PipelineConfig.clustering` do not exist.

Commit:

```powershell
git add -- tests/test_config.py
git commit -m "test: define clustering configuration contract"
```

- [ ] **Step 3: Implement minimal typed parsing and shared validation**

Add:

```python
@dataclass(frozen=True, slots=True)
class ClusteringConfig:
    enabled: bool = False
    identity: float = 0.95
    threads: int = 1
    memory_mb: int = 800
```

Add `clustering: ClusteringConfig = field(default_factory=ClusteringConfig)` to `PipelineConfig`. Parse only `enabled`, `identity`, `threads`, and `memory_mb` from an optional mapping. Reject unknown fields before reading values.

Use a shared validator:

```python
def validate_clustering_config(config: ClusteringConfig) -> None:
    if not isinstance(config.enabled, bool):
        raise ValueError("Clustering enabled must be a boolean.")
    if (
        isinstance(config.identity, bool)
        or not isinstance(config.identity, (int, float))
        or not 0.75 <= config.identity <= 1.0
    ):
        raise ValueError("Clustering identity must be a number between 0.75 and 1.0.")
    if (
        isinstance(config.threads, bool)
        or not isinstance(config.threads, int)
        or not 1 <= config.threads <= 256
    ):
        raise ValueError("Clustering threads must be an integer between 1 and 256.")
    if (
        isinstance(config.memory_mb, bool)
        or not isinstance(config.memory_mb, int)
        or config.memory_mb < 1
    ):
        raise ValueError("Clustering memory_mb must be a positive integer.")
```

Call it from `validate_pipeline_config` so YAML and direct configs have the same contract.

- [ ] **Step 4: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_config -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
```

Expected: all pass with no warnings.

Commit:

```powershell
git add -- qpcr_pipeline/config.py
git commit -m "feat: configure optional CD-HIT clustering"
```

---

### Task 2: Build the deterministic clustering service and parser

**Files:**
- Create: `qpcr_pipeline/clustering.py`
- Create: `tests/test_clustering.py`

**Interfaces:**
- Consumes: `ClusteringConfig`, `LocalSequenceRecord`, `DiscoverySet`, `EvaluationSet`
- Produces: `CdHitError`, `CdHitRunner`, `ClusterMember`, `SequenceCluster`, `ClusteringResult`
- Produces: `derive_word_length(identity: float) -> int`
- Produces: `cluster_sequences(records, evaluation_set, config, output_dir, *, runner=None) -> ClusteringResult`
- Task 3 supplies: production `SubprocessCdHitRunner`
- Task 4 consumes: `cluster_sequences` and `ClusteringResult.discovery_set`

- [ ] **Step 1: Add RED pure-function and parser tests**

Test every word-length boundary:

```python
self.assertEqual(derive_word_length(1.0), 10)
self.assertEqual(derive_word_length(0.95), 10)
self.assertEqual(derive_word_length(0.949), 8)
self.assertEqual(derive_word_length(0.90), 8)
self.assertEqual(derive_word_length(0.88), 7)
self.assertEqual(derive_word_length(0.85), 6)
self.assertEqual(derive_word_length(0.80), 5)
self.assertEqual(derive_word_length(0.75), 4)
```

Create a fake runner with the same contract as production:

```python
class FakeCdHitRunner:
    def __init__(self, representative_internal_ids, cluster_text):
        self.representative_internal_ids = representative_internal_ids
        self.cluster_text = cluster_text
        self.calls = []

    def run(self, input_path, output_path, config):
        self.calls.append((input_path, output_path, config))
        records = {record.id: record for record in SeqIO.parse(input_path, "fasta")}
        SeqIO.write(
            [records[record_id] for record_id in self.representative_internal_ids],
            output_path,
            "fasta",
        )
        Path(str(output_path) + ".clstr").write_text(
            self.cluster_text, encoding="utf-8"
        )
```

Use realistic `.clstr` text with output order different from Evaluation Set order:

```text
>Cluster 9
0 8nt, >geison-00000002... *
1 8nt, >geison-00000000... at +/99.00%
>Cluster 2
0 8nt, >geison-00000001... *
```

Assert final stable clusters, restored original IDs, representative membership,
identity/strand parsing, and Discovery Set ordered by original Evaluation Set
positions.

Add malformed cases for duplicate members, missing members, unknown IDs, zero/multiple representatives, representative absent from the cluster, malformed headers/member lines, and duplicate approved original IDs.

- [ ] **Step 2: Add RED disabled, empty, artifact, and atomicity tests**

Prove disabled and empty-enabled modes never call the runner:

```python
result = cluster_sequences(records, evaluation_set, ClusteringConfig(), outdir)
self.assertEqual(result.discovery_set, DiscoverySet(evaluation_set.sequence_ids))
self.assertEqual(json.loads(result.report_path.read_text())["clusters"], [])
```

For enabled clustering, assert:

- `discovery_set.fasta` contains original IDs/sequences;
- `clustering_report.json` has schema version 1, effective parameters, exact ordered Evaluation/Discovery IDs, before/after counts, clusters, and relative artifact names;
- `clustering/cd-hit-est.clstr` contains the raw fake output;
- internal IDs do not appear in final FASTA/JSON;
- a fake runner failure or invalid composition publishes no `COMPLETE` report;
- replacing final artifacts does not follow a pre-existing hardlink to an input artifact.

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_clustering -v
```

Expected: import failure because `qpcr_pipeline.clustering` does not exist.

Commit:

```powershell
git add -- tests/test_clustering.py
git commit -m "test: define traceable clustering service"
```

- [ ] **Step 4: Implement immutable models, strict parsing, and artifact publication**

Define:

```python
class CdHitError(RuntimeError):
    pass

class CdHitRunner(Protocol):
    def run(
        self,
        input_path: Path,
        output_path: Path,
        config: ClusteringConfig,
    ) -> None: ...

@dataclass(frozen=True, slots=True)
class ClusterMember:
    sequence_id: str
    representative: bool
    identity: float | None
    strand: Literal["+", "-"] | None

@dataclass(frozen=True, slots=True)
class SequenceCluster:
    cluster_id: str
    representative_id: str
    members: tuple[ClusterMember, ...]

@dataclass(frozen=True, slots=True)
class ClusteringResult:
    discovery_set: DiscoverySet
    clusters: tuple[SequenceCluster, ...]
    discovery_fasta_path: Path
    report_path: Path
    raw_cluster_path: Path | None
```

Implement strict parsing with full-line regular expressions. Convert raw internal
cluster labels into stable cluster IDs after sorting by earliest Evaluation Set
position. Require `set(parsed_ids) == set(expected_internal_ids)` and exactly one
occurrence per expected ID.

Write temporary input/output under `TemporaryDirectory`. Publish final FASTA, raw
cluster text, and JSON through unique temporary siblings and `Path.replace` only
after all validation succeeds. Write the report last.

When disabled or empty, build the Discovery Set directly without constructing the
default production runner.

- [ ] **Step 5: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_clustering -v
python -W default::ResourceWarning -m unittest tests.test_config tests.test_clustering -v
```

Expected: all pass, no handle warnings or real subprocess calls.

Commit:

```powershell
git add -- qpcr_pipeline/clustering.py
git commit -m "feat: build traceable Discovery Set clusters"
```

---

### Task 3: Add the production CD-HIT-EST runner and optional real-tool test

**Files:**
- Modify: `qpcr_pipeline/clustering.py`
- Modify: `tests/test_clustering.py`
- Create: `integration_tests/test_cdhit_est.py`

**Interfaces:**
- Produces: `SubprocessCdHitRunner(executable: str = "cd-hit-est")`
- `cluster_sequences(..., runner=None)` lazily constructs `SubprocessCdHitRunner`
- Consumes: Task 2 runner protocol and service

- [ ] **Step 1: Add RED production-runner tests**

Patch `shutil.which` and `subprocess.run`. Assert the exact structured argument list:

```python
[
    resolved_executable,
    "-i", str(input_path),
    "-o", str(output_path),
    "-c", "0.95",
    "-n", "10",
    "-d", "0",
    "-g", "1",
    "-r", "0",
    "-T", "4",
    "-M", "2048",
]
```

Assert `check=False`, `capture_output=True`, `text=True`, and no `shell` value other
than absent/false. Cover:

- missing executable raises an actionable `CdHitError` before `subprocess.run`;
- nonzero exit includes exit code and only bounded stderr;
- success without representative FASTA or `.clstr` raises `CdHitError`;
- paths containing spaces remain single arguments;
- default runner is created only for enabled, non-empty clustering.

- [ ] **Step 2: Run RED and commit unit tests**

Run:

```powershell
python -m unittest tests.test_clustering.SubprocessCdHitRunnerTests -v
```

Expected: import/attribute failure for `SubprocessCdHitRunner`.

Commit:

```powershell
git add -- tests/test_clustering.py
git commit -m "test: define CD-HIT subprocess boundary"
```

- [ ] **Step 3: Implement the production runner**

Resolve once per `run` call:

```python
executable = shutil.which(self.executable)
if executable is None:
    raise CdHitError(
        "cd-hit-est was not found on PATH; install CD-HIT or disable clustering."
    )
```

Invoke `subprocess.run(args, capture_output=True, text=True, check=False)`. On
nonzero exit, normalize whitespace and retain at most 2,000 stderr characters in the
error. Verify both output paths are regular files before returning.

- [ ] **Step 4: Add the optional real-tool integration test**

Create `integration_tests/test_cdhit_est.py`. Skip unless
`shutil.which("cd-hit-est")` returns a path. Build three records where two are
identical/highly similar and one differs, run enabled clustering at identity 0.95 in
a temporary directory, and assert:

```python
self.assertEqual(result.discovery_set.sequence_ids, ("seq-1", "seq-3"))
self.assertEqual(
    {member.sequence_id for member in result.clusters[0].members},
    {"seq-1", "seq-2"},
)
self.assertEqual(
    json.loads(result.report_path.read_text())["counts"],
    {"evaluation": 3, "discovery": 2},
)
```

Use sufficiently long deterministic nucleotide sequences for CD-HIT word matching.
Do not add the directory to pytest `testpaths`.

- [ ] **Step 5: Verify GREEN, separation, and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_clustering -v
python -m pytest --collect-only -q
python -W default::ResourceWarning -m unittest discover -s integration_tests -p "test_cdhit_est.py" -v
```

Expected on the current Windows host: unit tests pass; root pytest collects only
`tests/`; the real-tool test skips because `cd-hit-est` is absent.

Commit:

```powershell
git add -- qpcr_pipeline/clustering.py integration_tests/test_cdhit_est.py
git commit -m "feat: execute CD-HIT-EST through a safe runner"
```

---

### Task 4: Integrate Discovery Set clustering into pipeline outputs

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: `cluster_sequences`, `CdHitRunner`, `ClusteringResult`
- Produces: `run_pipeline(config, outdir, *, ncbi_client=None, cdhit_runner=None)`
- Produces: top-level `discovery_set.fasta`, `clustering_report.json`, and enabled raw cluster artifact
- Extends: `qc_report.json` with `discovery_set.sequence_ids`
- Preserves: `RunSummary` and `evaluation_set.sequence_ids` contracts

- [ ] **Step 1: Add RED pipeline integration tests**

Add a fake runner that records the internal input FASTA and emits two clusters. Run
the pipeline with one rejected input plus three approved inputs, two of which cluster
together. Assert:

```python
self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
self.assertEqual(
    json.loads((outdir / "clustering_report.json").read_text())["counts"],
    {"evaluation": 3, "discovery": 2},
)
```

Assert the fake runner input excludes the rejected sequence. Assert default-disabled
local and NCBI pipeline tests do not construct/call a runner, still create a
Discovery Set equal to Evaluation Set, and retain their existing summary IDs.

Add a direct invalid clustering config case and assert failure before output directory
creation.

- [ ] **Step 2: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_minimal_run -v
```

Expected: failures because `run_pipeline` lacks `cdhit_runner`, clustering is not
called, and Discovery artifacts/report fields do not exist.

Commit:

```powershell
git add -- tests/test_minimal_run.py
git commit -m "test: define Discovery Set pipeline integration"
```

- [ ] **Step 3: Implement pipeline routing**

After QC, select approved records without a dictionary that could hide duplicate
IDs:

```python
approved_ids = result.evaluation_set.sequence_ids
approved_id_set = set(approved_ids)
approved_records = tuple(
    record for record in records if record.sequence_id in approved_id_set
)
clustering = cluster_sequences(
    approved_records,
    result.evaluation_set,
    config.clustering,
    output_dir,
    runner=cdhit_runner,
)
```

The clustering service performs the final duplicate/membership validation. Add
`discovery_set` to the QC report from `clustering.discovery_set`. Do not change
`RunSummary.sequence_count`, `RunSummary.sequence_ids`, or Evaluation Set fields.

- [ ] **Step 4: Document configuration and outputs**

Extend `README.md` with:

- clustering is optional and disabled by default;
- YAML example with identity, threads, and memory;
- enabled runs require `cd-hit-est` on `PATH`;
- output meanings for `discovery_set.fasta`, `clustering_report.json`, and raw
  `.clstr`;
- explicit statement that the Evaluation Set remains the full QC-approved population.

- [ ] **Step 5: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_minimal_run tests.test_config tests.test_clustering -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
```

Expected: all standard tests pass with the one existing Windows symlink skip and no
real CD-HIT execution.

Commit:

```powershell
git add -- qpcr_pipeline/pipeline.py README.md
git commit -m "feat: publish clustered Discovery Set artifacts"
```

---

### Task 5: Review, verify, and publish issue #4

**Files:**
- Review: all changes from issue #4 base through HEAD
- No production changes unless review identifies a demonstrated defect

**Interfaces:**
- Produces: reviewed issue #4 implementation and fast-forwarded `origin/develop`

- [ ] **Step 1: Run fresh standard verification**

Run:

```powershell
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
python -W default::ResourceWarning -m unittest discover -s integration_tests -p "test_cdhit_est.py" -v
git diff --check 9f22cf3..HEAD
```

Expected: standard suites pass; real-tool test passes when `cd-hit-est` exists and
otherwise skips; no whitespace errors.

- [ ] **Step 2: Review the complete issue diff**

Review against GitHub issue #4 and the approved spec. Critical and Important findings
return through focused RED/GREEN fixes and scoped re-review. Minor observations are
recorded but do not block publication unless they expose an acceptance failure.

- [ ] **Step 3: Confirm remote fast-forward and push**

Run:

```powershell
git fetch origin develop
git merge-base --is-ancestor origin/develop HEAD
git rev-list --left-right --count origin/develop...HEAD
git push origin develop
git ls-remote origin refs/heads/develop
```

Expected: zero commits behind, a fast-forward push, and remote `develop` equal to
local HEAD. Do not open a PR, modify `main`, or close/mutate the issue.
