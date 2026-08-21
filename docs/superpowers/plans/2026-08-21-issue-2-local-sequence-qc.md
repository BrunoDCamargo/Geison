# Issue 2 Local Sequence QC Completion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete issue #2 by making local FASTA and GenBank inputs flow through traceable QC into the target and evaluation sequence sets.

**Architecture:** Normalize both Biopython input formats into `LocalSequenceRecord`, including available source metadata, and keep QC format-independent. Extend the existing configuration and pipeline orchestration to select the loader, pass optional QC thresholds, and serialize a dedicated QC report while preserving the current run summary behavior for approved sequences.

**Tech Stack:** Python 3.10+, Biopython, PyYAML, `unittest`

**Spec:** `C:/Users/brunodc.local/Desktop/handoff-geison-codex.md` and GitHub issue `BrunoDCamargo/Geison#2`

## Global Constraints

- Work only on the authorized `develop` branch.
- Keep all sequence processing offline and deterministic.
- Preserve available metadata, but never reject a sequence only because metadata is absent.
- Keep QC reason codes unchanged and traceable.
- `EvaluationSet` must contain every and only QC-approved sequence.
- Use strict RED/GREEN TDD for each behavior; commit tests only after observing the expected RED, then commit the minimal GREEN.
- Do not open a PR, merge `main`, or close or mutate issues.

---

### Task 1: Close FASTA parser resources deterministically

**Files:**
- Modify: `tests/test_fasta_qc.py`
- Modify: `qpcr_pipeline/local_input.py`

**Interfaces:**
- Consumes: `load_fasta(path: str | Path) -> tuple[LocalSequenceRecord, ...]`
- Produces: the same loader contract with its file handle closed before return

- [ ] **Step 1: Add a regression test that catches the leak**

Add a test that creates a real temporary FASTA, records `ResourceWarning` instances with `warnings.catch_warnings(record=True)`, calls `load_fasta`, forces collection with `gc.collect()`, and asserts both the parsed IDs and an empty list of captured `ResourceWarning` messages. The production regression that must make this test fail is passing a path directly to Biopython 1.88 without closing its parser-owned stream.

```python
def test_load_fasta_closes_its_input_file(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        fasta_path = Path(tmpdir) / "target.fasta"
        fasta_path.write_text(">seq-1\nACGT\n", encoding="utf-8")

        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter("always", ResourceWarning)
            records = load_fasta(fasta_path)
            gc.collect()

    self.assertEqual(tuple(record.sequence_id for record in records), ("seq-1",))
    self.assertEqual(
        [str(item.message) for item in captured if issubclass(item.category, ResourceWarning)],
        [],
    )
```

- [ ] **Step 2: Run the focused test and verify RED**

Run: `python -W default::ResourceWarning -m unittest tests.test_fasta_qc.FastaQcTests.test_load_fasta_closes_its_input_file -v`

Expected: FAIL because the current `SeqIO.parse(Path, "fasta")` leaves the parser-owned stream unclosed.

- [ ] **Step 3: Commit the RED test**

Run: `git add -- tests/test_fasta_qc.py && git commit -m "test: reproduce FASTA input resource leak"`

- [ ] **Step 4: Implement the minimal resource-lifecycle fix**

Open the FASTA path explicitly as a UTF-8 text stream with a `with` block and pass that stream to `SeqIO.parse`. Materialize the normalized tuple inside the block so the stream is deterministically closed before `load_fasta` returns.

```python
def load_fasta(path: str | Path) -> tuple[LocalSequenceRecord, ...]:
    fasta_path = Path(path)
    with fasta_path.open(encoding="utf-8") as handle:
        return tuple(
            LocalSequenceRecord(sequence_id=record.id, sequence=str(record.seq).upper())
            for record in SeqIO.parse(handle, "fasta")
        )
```

- [ ] **Step 5: Verify GREEN and the focused QC suite**

Run: `python -W default::ResourceWarning -m unittest tests.test_fasta_qc -v`

Expected: PASS with no `ResourceWarning` output.

- [ ] **Step 6: Commit the GREEN implementation**

Run: `git add -- qpcr_pipeline/local_input.py && git commit -m "fix: close FASTA input stream"`

---

### Task 2: Normalize GenBank records and preserve metadata

**Files:**
- Create: `tests/test_local_input.py`
- Modify: `qpcr_pipeline/local_input.py`

**Interfaces:**
- Consumes: Biopython `SeqRecord` objects parsed from local files
- Produces: `LocalSequenceRecord(sequence_id: str, sequence: str, metadata: Mapping[str, object])`, `load_fasta(path)`, `load_genbank(path)`, and `load_local_sequences(path, file_format)`

- [ ] **Step 1: Add RED tests for normalized GenBank input and metadata**

Create deterministic temporary GenBank files with real Biopython `SeqRecord` values. Assert that `load_genbank` uppercases sequence text and preserves literal `name`, `description`, `annotations`, `dbxrefs`, and `features` values under `LocalSequenceRecord.metadata`. Add a FASTA assertion for its available name and description. The production changes that must make these tests fail are the missing `metadata` field and missing GenBank loader.

```python
source = SeqRecord(
    Seq("acgtn"),
    id="GB001",
    name="GENE1",
    description="Synthetic target",
)
source.annotations.update({"molecule_type": "DNA", "organism": "Synthetic construct"})
source.dbxrefs = ["BioProject:PRJ1"]
source.features = [
    SeqFeature(FeatureLocation(0, 5), type="gene", qualifiers={"gene": ["abc"]})
]
SeqIO.write((source,), genbank_path, "genbank")

loaded = load_genbank(genbank_path)

self.assertEqual(loaded[0].sequence_id, "GB001")
self.assertEqual(loaded[0].sequence, "ACGTN")
self.assertEqual(loaded[0].metadata["name"], "GENE1")
self.assertEqual(loaded[0].metadata["description"], "Synthetic target")
self.assertEqual(loaded[0].metadata["annotations"]["organism"], "Synthetic construct")
self.assertEqual(loaded[0].metadata["dbxrefs"], ("BioProject:PRJ1",))
self.assertEqual(loaded[0].metadata["features"][0].type, "gene")
```

- [ ] **Step 2: Add a RED test proving missing metadata does not reject**

Create a minimal GenBank record whose only required annotation is `molecule_type`, load it, pass it to `evaluate_sequences`, and assert `QCStatus.ACCEPTED` with empty reason codes and its ID present in the evaluation set. The expected value must be literal and must not be computed by loader helpers.

```python
source = SeqRecord(Seq("ACGT"), id="metadata-light")
source.annotations["molecule_type"] = "DNA"
SeqIO.write((source,), genbank_path, "genbank")

result = evaluate_sequences(load_genbank(genbank_path))

self.assertEqual(result.records[0].status, QCStatus.ACCEPTED)
self.assertEqual(result.records[0].reason_codes, ())
self.assertEqual(result.evaluation_set.sequence_ids, ("metadata-light",))
```

- [ ] **Step 3: Run the focused tests and verify RED**

Run: `python -m unittest tests.test_local_input -v`

Expected: FAIL because `load_genbank`, `load_local_sequences`, and record metadata do not exist.

- [ ] **Step 4: Commit the RED tests**

Run: `git add -- tests/test_local_input.py && git commit -m "test: define GenBank metadata input contract"`

- [ ] **Step 5: Implement the minimal normalized loaders**

Add a default-empty metadata mapping to `LocalSequenceRecord`. Factor one private normalizer that copies the Biopython record's `name`, `description`, `annotations`, `dbxrefs`, and `features`; use it from both explicit UTF-8 context-managed loaders. Add `load_local_sequences` with exactly the supported format strings `"fasta"` and `"genbank"`, raising `ValueError` for anything else.

```python
@dataclass(frozen=True, slots=True)
class LocalSequenceRecord:
    sequence_id: str
    sequence: str
    metadata: Mapping[str, object] = field(default_factory=dict)


def _normalize(record: SeqRecord) -> LocalSequenceRecord:
    return LocalSequenceRecord(
        sequence_id=record.id,
        sequence=str(record.seq).upper(),
        metadata={
            "name": record.name,
            "description": record.description,
            "annotations": dict(record.annotations),
            "dbxrefs": tuple(record.dbxrefs),
            "features": tuple(record.features),
        },
    )


def load_genbank(path: str | Path) -> tuple[LocalSequenceRecord, ...]:
    with Path(path).open(encoding="utf-8") as handle:
        return tuple(_normalize(record) for record in SeqIO.parse(handle, "genbank"))


def load_local_sequences(
    path: str | Path,
    file_format: Literal["fasta", "genbank"],
) -> tuple[LocalSequenceRecord, ...]:
    if file_format == "fasta":
        return load_fasta(path)
    if file_format == "genbank":
        return load_genbank(path)
    raise ValueError(f"Unsupported local sequence format: {file_format}")
```

- [ ] **Step 6: Verify GREEN**

Run: `python -W default::ResourceWarning -m unittest tests.test_local_input tests.test_fasta_qc -v`

Expected: PASS with no warnings.

- [ ] **Step 7: Commit the GREEN implementation**

Run: `git add -- qpcr_pipeline/local_input.py && git commit -m "feat: load local GenBank records with metadata"`

---

### Task 3: Integrate local input and traceable QC into the pipeline

**Files:**
- Modify: `tests/test_config.py`
- Modify: `tests/test_minimal_run.py`
- Modify: `qpcr_pipeline/config.py`
- Modify: `qpcr_pipeline/pipeline.py`

**Interfaces:**
- Consumes: `load_local_sequences(path, file_format)` and `evaluate_sequences(records, *, min_length, max_ambiguous_fraction, expected_length, length_tolerance_fraction)`
- Produces: `PipelineConfig` with exactly one local FASTA or GenBank path plus optional QC thresholds; `run_pipeline` writes `run_summary.json` and `qc_report.json`

- [ ] **Step 1: Add RED configuration tests**

Add literal YAML cases showing `input.genbank` is accepted, exactly one of `input.fasta` and `input.genbank` is required, and optional `qc.min_length`, `qc.max_ambiguous_fraction`, `qc.expected_length`, and `qc.length_tolerance_fraction` values are loaded. Preserve the current FASTA configuration contract.

```python
config_path.write_text(
    "target:\n"
    "  name: synthetic-target\n"
    "input:\n"
    "  genbank: tests/fixtures/target.gb\n"
    "qc:\n"
    "  min_length: 100\n"
    "  max_ambiguous_fraction: 0.05\n"
    "  expected_length: 150\n"
    "  length_tolerance_fraction: 0.10\n",
    encoding="utf-8",
)

config = load_config(config_path)

self.assertIsNone(config.input_fasta)
self.assertEqual(config.input_genbank, Path("tests/fixtures/target.gb"))
self.assertEqual(config.selected_input, (Path("tests/fixtures/target.gb"), "genbank"))
self.assertEqual(config.qc.min_length, 100)
self.assertEqual(config.qc.max_ambiguous_fraction, 0.05)
self.assertEqual(config.qc.expected_length, 150)
self.assertEqual(config.qc.length_tolerance_fraction, 0.10)
```

- [ ] **Step 2: Run configuration tests and verify RED**

Run: `python -m unittest tests.test_config -v`

Expected: FAIL because GenBank selection and QC configuration are not represented.

- [ ] **Step 3: Implement minimal configuration support and verify GREEN**

Represent QC thresholds in a frozen `QCConfig` with `None` defaults. Keep `PipelineConfig.input_fasta` compatible, add `input_genbank`, require exactly one configured input in `load_config`, and expose the selected path and format to the pipeline without filesystem or network access.

```python
@dataclass(frozen=True, slots=True)
class QCConfig:
    min_length: int | None = None
    max_ambiguous_fraction: float | None = None
    expected_length: int | None = None
    length_tolerance_fraction: float | None = None


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path | None = None
    input_genbank: Path | None = None
    qc: QCConfig = field(default_factory=QCConfig)

    @property
    def selected_input(self) -> tuple[Path, Literal["fasta", "genbank"]]:
        if self.input_fasta is not None:
            return self.input_fasta, "fasta"
        if self.input_genbank is not None:
            return self.input_genbank, "genbank"
        raise ValueError("Exactly one local sequence input must be configured.")
```

Run: `python -m unittest tests.test_config -v`

Expected: PASS.

- [ ] **Step 4: Commit the configuration RED/GREEN slice**

Run: `git add -- tests/test_config.py qpcr_pipeline/config.py && git commit -m "feat: configure local sequence input and QC"`

- [ ] **Step 5: Add a RED end-to-end QC report test**

Extend the real CLI/pipeline test with a deterministic local file containing approved and rejected sequences. Assert that `run_summary.json` lists only all approved IDs and that `qc_report.json` contains, for every input ID in order, literal `status` strings and literal `reason_codes` arrays. Assert that both serialized target and evaluation set IDs equal all approved IDs. The production regression that must fail this test is bypassing `evaluate_sequences` or dropping an approved sequence from `EvaluationSet`.

```python
self.assertEqual(summary["sequence_ids"], ["accepted-1", "accepted-2"])
self.assertEqual(
    qc_report["records"],
    [
        {"sequence_id": "accepted-1", "status": "ACCEPTED", "reason_codes": []},
        {"sequence_id": "invalid", "status": "REJECTED", "reason_codes": ["INVALID_NUCLEOTIDE"]},
        {"sequence_id": "accepted-2", "status": "ACCEPTED", "reason_codes": []},
    ],
)
self.assertEqual(qc_report["target_sequence_set"]["sequence_ids"], ["accepted-1", "accepted-2"])
self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["accepted-1", "accepted-2"])
```

- [ ] **Step 6: Run the focused pipeline test and verify RED**

Run: `python -m unittest tests.test_minimal_run -v`

Expected: FAIL because the pipeline still reads FASTA headers directly and does not write `qc_report.json`.

- [ ] **Step 7: Implement minimal pipeline integration**

Replace `_read_fasta_ids` with the selected normalized loader and `evaluate_sequences`, passing every configured QC threshold. Build the summary from `result.evaluation_set.sequence_ids`. Serialize `qc_report.json` with records in source order and with `QCStatus.value`, reason-code lists, target IDs, and evaluation IDs. Keep output deterministic and UTF-8 encoded.

```python
input_path, input_format = config.selected_input
records = load_local_sequences(input_path, input_format)
result = evaluate_sequences(
    records,
    min_length=config.qc.min_length,
    max_ambiguous_fraction=config.qc.max_ambiguous_fraction,
    expected_length=config.qc.expected_length,
    length_tolerance_fraction=config.qc.length_tolerance_fraction,
)
approved_ids = result.evaluation_set.sequence_ids

qc_report = {
    "records": [
        {
            "sequence_id": record.sequence_id,
            "status": record.status.value,
            "reason_codes": list(record.reason_codes),
        }
        for record in result.records
    ],
    "target_sequence_set": {"sequence_ids": list(result.target_sequence_set.sequence_ids)},
    "evaluation_set": {"sequence_ids": list(approved_ids)},
}
```

- [ ] **Step 8: Verify GREEN and commit**

Run: `python -W default::ResourceWarning -m unittest tests.test_minimal_run tests.test_cli -v`

Expected: PASS with no warnings.

Run: `git add -- tests/test_minimal_run.py qpcr_pipeline/pipeline.py && git commit -m "feat: run traceable QC for local inputs"`

---

### Task 4: Completion gate

**Files:**
- Verify only; modify files solely if review findings require a new RED/GREEN cycle

**Interfaces:**
- Consumes: all issue #2 implementation commits
- Produces: complete green suite, clean review, and updated remote `develop`

- [ ] **Step 1: Run the complete suite with warnings visible**

Run: `python -W default::ResourceWarning -m unittest discover -s tests -v`

Expected: all tests PASS with no warnings or errors.

- [ ] **Step 2: Review the complete issue #2 diff**

Review against GitHub issue #2, the global constraints above, and the full branch diff. Any Critical or Important finding returns to a focused RED/GREEN fix and scoped re-review.

- [ ] **Step 3: Publish only the authorized branch**

After verification and clean review, push the current `develop` commits to `origin/develop`. Do not create a PR, merge, or mutate the issue.
