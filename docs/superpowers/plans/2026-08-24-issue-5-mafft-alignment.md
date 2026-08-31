# Traceable MAFFT Alignment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Align the Discovery Set with MAFFT, record reference selection and reverse-complement transformations, and publish a 1-based alignment-to-reference coordinate map while preserving Discovery and Evaluation Set contracts.

**Architecture:** Add typed alignment configuration and an isolated `qpcr_pipeline.alignment` service behind an injectable runner protocol. The service selects a deterministic Discovery reference, uses stable internal IDs with MAFFT direction adjustment, strictly restores original IDs and orientations, constructs coordinates from the aligned reference, and atomically publishes FASTA, TSV, and JSON artifacts; the pipeline calls it after clustering.

**Tech Stack:** Python 3.10+, standard-library `subprocess`/`tempfile`/`shutil`/`json`/`pathlib`/`uuid`, Biopython `SeqIO` and `Seq.reverse_complement`, PyYAML, `unittest`, MAFFT external binary for the optional integration test.

**Spec:** `docs/superpowers/specs/2026-08-24-issue-5-mafft-alignment-design.md`

## Global Constraints

- Alignment defaults to disabled; existing runs must not require MAFFT.
- Enabled alignment consumes exactly the unchanged Discovery Set and never mutates, filters, or reorders Discovery or Evaluation membership.
- An explicit reference is an existing Discovery Set sequence ID; no external reference acquisition is added.
- Automatic selection minimizes ambiguous-base fraction, then maximizes sequence length, then preserves Discovery order.
- Multi-sequence execution uses `--auto --nuc --inputorder --adjustdirectionaccurately --thread <threads> --threadit 0 --quiet` through an argument list with no shell.
- Final coordinates are 1-based; reference gaps have empty/null reference positions and bases.
- Final FASTA, TSV, and JSON artifacts use original IDs and atomic sibling replacement; the report is written last.
- Standard tests remain offline; the real MAFFT test lives under `integration_tests/` and skips when `mafft` is absent.
- Empty and singleton enabled runs do not invoke MAFFT; disabled runs publish `SKIPPED` and remove only exact stale alignment data artifacts.
- No external canonical reference, consensus, primer discovery, variant analysis, visualization, checkpointing, or dependency `doctor` behavior is added.

---

### Task 1: Parse and validate alignment configuration

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Produces: `AlignmentConfig(enabled: bool = False, threads: int = 1, reference_id: str | None = None)`
- Produces: `validate_alignment_config(config: AlignmentConfig) -> None`
- Extends: `PipelineConfig.alignment: AlignmentConfig`
- Consumes later: Task 2 alignment service and Task 4 pipeline integration

- [ ] **Step 1: Add RED defaults and YAML tests**

Import `AlignmentConfig`. Prove omitted configuration uses defaults and explicit YAML parses exactly:

```python
config = self._load_yaml(
    "target:\n  name: target\n"
    f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
    "alignment:\n"
    "  enabled: true\n"
    "  threads: 4\n"
    "  reference_id: seq-1\n"
)
self.assertEqual(
    config.alignment,
    AlignmentConfig(enabled=True, threads=4, reference_id="seq-1"),
)
```

The production mutation these tests catch is ignoring the section or changing its defaults.

- [ ] **Step 2: Add RED invalid YAML and direct-construction tables**

Add invalid YAML cases for an unknown key, integer `enabled`, boolean/fractional/zero/257 threads, numeric reference ID, empty reference ID, and whitespace-only reference ID. Add direct cases and prove `selected_input` validates before pipeline work:

```python
for alignment in (
    AlignmentConfig(enabled=1),
    AlignmentConfig(threads=True),
    AlignmentConfig(threads=0),
    AlignmentConfig(reference_id=7),
    AlignmentConfig(reference_id=" "),
):
    with self.subTest(alignment=alignment):
        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            alignment=alignment,
        )
        with self.assertRaises(ValueError):
            _ = config.selected_input
```

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_config.PipelineConfigTests -v
```

Expected: import failure because `AlignmentConfig` does not exist.

Commit:

```powershell
git add -- tests/test_config.py
git commit -m "test: define alignment configuration contract"
```

- [ ] **Step 4: Implement typed parsing and shared validation**

Add the frozen, slotted dataclass and `alignment: AlignmentConfig = field(default_factory=AlignmentConfig)` to `PipelineConfig`. Parse only `enabled`, `threads`, and `reference_id` from an optional mapping, rejecting unknown fields before values are read.

Use shared validation:

```python
def validate_alignment_config(config: AlignmentConfig) -> None:
    if not isinstance(config, AlignmentConfig):
        raise ValueError("Alignment configuration must be an AlignmentConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Alignment enabled must be a boolean.")
    if (
        isinstance(config.threads, bool)
        or not isinstance(config.threads, int)
        or not 1 <= config.threads <= 256
    ):
        raise ValueError("Alignment threads must be an integer between 1 and 256.")
    if config.reference_id is not None and (
        not isinstance(config.reference_id, str) or not config.reference_id.strip()
    ):
        raise ValueError("Alignment reference_id must be a non-blank string when configured.")
```

Call it from `validate_pipeline_config`.

- [ ] **Step 5: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_config -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
```

Expected: all pass with only the existing Windows symlink skip and no warnings.

Commit:

```powershell
git add -- qpcr_pipeline/config.py
git commit -m "feat: configure optional MAFFT alignment"
```

---

### Task 2: Build the deterministic alignment service and artifacts

**Files:**
- Create: `qpcr_pipeline/alignment.py`
- Create: `tests/test_alignment.py`

**Interfaces:**
- Consumes: `AlignmentConfig`, `LocalSequenceRecord`, `DiscoverySet`
- Produces: `MafftError`, `MafftRunner`, `AlignedSequence`, `AlignmentCoordinate`, `AlignmentResult`
- Produces: `align_discovery(records, discovery_set, config, output_dir, *, runner=None) -> AlignmentResult`
- Task 3 supplies: production `SubprocessMafftRunner`
- Task 4 consumes: `align_discovery` and the result's status/reference fields

- [ ] **Step 1: Add RED fake-runner happy-path tests**

Create a fake that writes a literal aligned FASTA and records the temporary input:

```python
class FakeMafftRunner:
    def __init__(self, output_records):
        self.output_records = output_records
        self.calls = []
        self.input_records = []

    def run(self, input_path, output_path, config):
        self.calls.append((input_path, output_path, config))
        self.input_records = [
            (record.id, str(record.seq))
            for record in SeqIO.parse(input_path, "fasta")
        ]
        SeqIO.write(
            [
                SeqRecord(Seq(sequence), id=sequence_id, description="")
                for sequence_id, sequence in self.output_records
            ],
            output_path,
            "fasta",
        )
```

Use Discovery records in order `ref`, `reverse`, `other`, an explicit `reference_id="ref"`, output records in a different order, and `_R_geison-00000001` for a sequence whose ungapped aligned value is the reverse complement of its original. Assert:

- MAFFT temporary input puts the reference first and uses stable IDs derived from Discovery positions;
- `AlignmentResult.discovery_set` is exactly the input Discovery Set;
- final sequences restore original IDs in Discovery order;
- orientations are `forward`, `reverse_complemented`, `forward`;
- internal and `_R_` IDs do not appear in the final FASTA or JSON.

- [ ] **Step 2: Add RED automatic-reference and coordinate tests**

Prove the automatic key in separate table-driven cases: lower ambiguity fraction wins, then greater length, then earlier Discovery position. Assert report mode `automatic` and the exact rule name `lowest_ambiguity_fraction_then_longest_then_discovery_order`.

Use an aligned reference `A-CGT` and assert literal coordinates:

```python
self.assertEqual(
    result.coordinates,
    (
        AlignmentCoordinate(1, 1, "A"),
        AlignmentCoordinate(2, None, None),
        AlignmentCoordinate(3, 2, "C"),
        AlignmentCoordinate(4, 3, "G"),
        AlignmentCoordinate(5, 4, "T"),
    ),
)
```

Assert `coordinate_map.tsv` has the exact header and blank reference fields on row 2.

- [ ] **Step 3: Add RED strict-validation cases**

Add cases for:

- duplicate Discovery IDs and duplicate record IDs;
- record/Discovery membership mismatch;
- explicit reference missing from Discovery;
- empty input record;
- missing, unknown, or duplicate MAFFT output IDs;
- double `_R_` prefix and a reversed reference;
- unequal alignment lengths, all-gap record, and invalid output symbol;
- forward output not equal to the original after gap removal;
- reversed output not equal to the original's IUPAC reverse complement;
- runner failure and missing/empty output file.

Each case must assert no new `COMPLETE` report in a fresh output directory.

- [ ] **Step 4: Add RED disabled, empty, singleton, report, and atomicity tests**

Prove disabled, enabled-empty, and enabled-singleton paths never call a failing runner. Assert:

- disabled status `SKIPPED`, only `alignment/alignment_report.json`, null artifacts, and exact stale FASTA/TSV removal while preserving unrelated siblings;
- enabled empty status `COMPLETE`, empty FASTA, header-only TSV, null reference, and no orientations;
- enabled singleton identity alignment, automatic or matching explicit reference, forward orientation, and one coordinate per base;
- schema version 1, effective parameters, ordered IDs, reference decision, orientations, counts, and relative artifact names;
- replacing final artifacts does not follow a pre-existing hardlink to a source file;
- report publication occurs after data artifacts, so a replace failure does not publish a new `COMPLETE` report.

- [ ] **Step 5: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_alignment -v
```

Expected: import failure because `qpcr_pipeline.alignment` does not exist.

Commit:

```powershell
git add -- tests/test_alignment.py
git commit -m "test: define traceable MAFFT alignment service"
```

- [ ] **Step 6: Implement immutable models, selection, parsing, coordinates, and publication**

Define the public interfaces exactly as the spec. Validate records before branching. For automatic reference selection, use:

```python
def selection_key(item):
    position, record = item
    ambiguous = sum(base not in "ACGT" for base in record.sequence)
    return (ambiguous / len(record.sequence), -len(record.sequence), position)
```

For two or more enabled records, require an injected runner in this task and raise an actionable `MafftError` when it is absent; Task 3 replaces this temporary boundary with lazy production construction.

Write temporary input/output under `TemporaryDirectory`. Put the chosen reference first while keeping internal IDs tied to original Discovery positions. Parse output with `SeqIO`, split exactly one optional `_R_` prefix, validate complete membership and sequence identity, restore Discovery order, and build the coordinate tuple from the aligned forward reference.

Write `alignment/discovery_alignment.fasta`, `alignment/coordinate_map.tsv`, and `alignment/alignment_report.json` through unique temporary siblings and `Path.replace`. Write the report last. Disabled publication unlinks only the exact stale FASTA/TSV paths and publishes `SKIPPED` last.

- [ ] **Step 7: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_alignment -v
python -W default::ResourceWarning -m unittest tests.test_config tests.test_alignment -v
```

Expected: all pass with no real subprocess call or warnings.

Commit:

```powershell
git add -- qpcr_pipeline/alignment.py
git commit -m "feat: build traceable Discovery Set alignment"
```

---

### Task 3: Add the production MAFFT runner and optional real-tool test

**Files:**
- Modify: `qpcr_pipeline/alignment.py`
- Modify: `tests/test_alignment.py`
- Create: `integration_tests/test_mafft.py`

**Interfaces:**
- Produces: `SubprocessMafftRunner(executable: str = "mafft")`
- `align_discovery(..., runner=None)` lazily constructs `SubprocessMafftRunner` only for enabled inputs with two or more records
- Consumes: Task 2 runner protocol and service

- [ ] **Step 1: Add RED production-runner tests**

Patch `shutil.which` and `subprocess.run`. Assert the exact list:

```python
[
    resolved_executable,
    "--auto",
    "--nuc",
    "--inputorder",
    "--adjustdirectionaccurately",
    "--thread", "4",
    "--threadit", "0",
    "--quiet",
    str(input_path),
]
```

Assert `capture_output=True`, `text=True`, `check=False`, and no true `shell` option. Cover paths with spaces, a missing configured executable with install-or-disable guidance, a nonzero exit with at most 2,000 normalized stderr characters, and success with empty stdout. On success the runner writes stdout exactly to `output_path`.

- [ ] **Step 2: Add RED lazy-default tests and commit**

Patch `SubprocessMafftRunner`. Call disabled, empty enabled, singleton enabled, and two-record enabled services. Assert the factory and runner are called only for the last case.

Run:

```powershell
python -m unittest tests.test_alignment.SubprocessMafftRunnerTests tests.test_alignment.DefaultRunnerTests -v
```

Expected: import/attribute failure for `SubprocessMafftRunner`.

Commit:

```powershell
git add -- tests/test_alignment.py
git commit -m "test: define MAFFT subprocess boundary"
```

- [ ] **Step 3: Implement the safe production runner**

Resolve once per call:

```python
executable = shutil.which(self.executable)
if executable is None:
    raise MafftError(
        f"{self.executable!r} was not found on PATH; install MAFFT or disable alignment."
    )
```

Invoke `subprocess.run(args, capture_output=True, text=True, check=False)`. Normalize and bound stderr on failure. Reject empty/whitespace stdout. Write stdout to the temporary service output path using UTF-8 with newline preservation.

- [ ] **Step 4: Add the optional real MAFFT integration test**

Create `integration_tests/test_mafft.py` with `@unittest.skipUnless(shutil.which("mafft"), "mafft is not installed")`. Build a deterministic reference of at least 80 nt, its exact reverse complement, and one forward variant. Enable alignment with `reference_id="seq-1"` and two threads in a temporary directory. Assert:

```python
self.assertEqual(result.reference_id, "seq-1")
self.assertEqual(
    [sequence.orientation for sequence in result.sequences],
    ["forward", "reverse_complemented", "forward"],
)
self.assertEqual(result.coordinates[0].alignment_position, 1)
self.assertEqual(
    result.coordinates[-1].reference_position,
    len(reference),
)
```

Also assert final IDs, artifact existence, and report counts. Do not add `integration_tests` to pytest `testpaths`.

- [ ] **Step 5: Verify GREEN, test separation, and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_alignment -v
python -m pytest --collect-only -q
python -W default::ResourceWarning -m unittest discover -s integration_tests -p "test_mafft.py" -v
```

Expected on the current Windows host: unit tests pass; root pytest collects only `tests/`; the real test skips because MAFFT is absent.

Commit:

```powershell
git add -- qpcr_pipeline/alignment.py integration_tests/test_mafft.py
git commit -m "feat: execute MAFFT through a safe runner"
```

---

### Task 4: Integrate alignment into pipeline outputs and documentation

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `README.md`

**Interfaces:**
- Consumes: `align_discovery`, `MafftRunner`, `AlignmentResult`
- Produces: `run_pipeline(config, outdir, *, ncbi_client=None, cdhit_runner=None, mafft_runner=None)`
- Produces: `alignment/alignment_report.json` on every run and enabled alignment FASTA/TSV
- Extends: `qc_report.json.alignment` with status, reference ID, and reference mode
- Preserves: `RunSummary`, Evaluation Set, Discovery Set, and clustering artifact contracts

- [ ] **Step 1: Add RED enabled pipeline integration test**

Use one rejected input and three approved inputs. Inject a fake CD-HIT runner that reduces Discovery to `s1` and `s3`, then inject a fake MAFFT runner with explicit reference `s3`. Assert MAFFT's internal input contains only Discovery representatives and puts the stable internal ID for `s3` first:

```python
self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["s1", "s2", "s3"])
self.assertEqual(qc_report["discovery_set"]["sequence_ids"], ["s1", "s3"])
self.assertEqual(qc_report["alignment"], {
    "status": "COMPLETE",
    "reference_id": "s3",
    "reference_mode": "explicit",
})
self.assertEqual(summary.sequence_ids, ["s1", "s2", "s3"])
```

Assert final aligned FASTA IDs are `s1`, `s3` in Discovery order and coordinate TSV columns are present.

- [ ] **Step 2: Add RED disabled and validation pipeline tests**

Extend local and NCBI default tests to patch `SubprocessMafftRunner`, assert it is never constructed, and assert a `SKIPPED` alignment report plus top-level QC traceability. Add a directly invalid `AlignmentConfig` case and assert failure before output directory creation. Add enabled explicit-reference-not-in-Discovery coverage proving no MAFFT call and no new `COMPLETE` report.

- [ ] **Step 3: Run RED and commit tests**

Run:

```powershell
python -m unittest tests.test_minimal_run -v
```

Expected: import/signature/artifact failures because alignment is not routed through `run_pipeline`.

Commit:

```powershell
git add -- tests/test_minimal_run.py
git commit -m "test: define alignment pipeline integration"
```

- [ ] **Step 4: Implement pipeline routing**

After clustering, select the Discovery records without a dictionary:

```python
discovery_ids = clustering.discovery_set.sequence_ids
discovery_id_set = set(discovery_ids)
discovery_records = tuple(
    record for record in approved_records if record.sequence_id in discovery_id_set
)
alignment = align_discovery(
    discovery_records,
    clustering.discovery_set,
    config.alignment,
    output_dir,
    runner=mafft_runner,
)
```

Add the keyword-only injection and the QC `alignment` object. Do not change summary or set semantics.

- [ ] **Step 5: Document configuration, artifacts, and coordinate meaning**

Extend `README.md` with the default-disabled alignment section, YAML example, explicit/automatic reference behavior, MAFFT-on-PATH requirement, accurate direction normalization, artifact descriptions, and the exact statement that both coordinate columns are 1-based while reference gaps have blank reference fields.

- [ ] **Step 6: Verify GREEN and commit**

Run:

```powershell
python -W default::ResourceWarning -m unittest tests.test_minimal_run tests.test_config tests.test_alignment tests.test_clustering -v
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
```

Expected: all standard tests pass with the existing Windows symlink skip and no real MAFFT call.

Commit:

```powershell
git add -- qpcr_pipeline/pipeline.py README.md
git commit -m "feat: publish reference-coordinate alignment artifacts"
```

---

### Task 5: Review, verify, and publish issue #5

**Files:**
- Review: all changes from issue #5 base `cfe0a95` through HEAD
- No production changes unless review demonstrates a defect

**Interfaces:**
- Produces: reviewed issue #5 implementation and fast-forwarded `origin/develop`

- [ ] **Step 1: Run fresh final verification**

Run:

```powershell
python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v
python -m pytest -q
python -W default::ResourceWarning -m unittest discover -s integration_tests -p "test_mafft.py" -v
git diff --check cfe0a95..HEAD
```

Expected: standard suites pass; the real MAFFT test passes when installed and otherwise skips; no whitespace errors.

- [ ] **Step 2: Review the complete issue diff**

Review against GitHub issue #5 and the approved design. Critical and Important findings return through focused RED/GREEN fixes and scoped re-review. Minor observations are recorded and triaged by the whole-issue review.

- [ ] **Step 3: Confirm remote fast-forward and publish**

Run:

```powershell
git fetch origin develop
git merge-base --is-ancestor origin/develop HEAD
git rev-list --left-right --count origin/develop...HEAD
git push origin develop
git ls-remote origin refs/heads/develop
```

Expected: zero commits behind, fast-forward push, and remote `develop` equal to local HEAD. The standing authorization is direct `origin/develop` publication only. Do not open a PR, modify `main`, or close/mutate GitHub issues.
