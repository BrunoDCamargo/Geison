# Issue 3 Reproducible NCBI Acquisition Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a resumable, reproducible NCBI nucleotide acquisition path that accepts a query or accession list, freezes the exact dataset, and feeds the existing QC pipeline.

**Architecture:** Extend configuration with one typed NCBI source, implement a client boundary and checksummed manifest/batch service in `qpcr_pipeline.ncbi`, then route the resulting GenBank file through the existing normalizer and QC. Keep production Entrez calls behind dependency injection so the standard suite remains offline and deterministic; keep the optional live test outside `tests/`.

**Tech Stack:** Python 3.10+, Biopython `Bio.Entrez`/`SeqIO`, PyYAML, standard-library JSON/hashlib/pathlib/urllib, `unittest`

**Spec:** `docs/superpowers/specs/2026-08-21-issue-3-ncbi-acquisition-design.md`

## Global Constraints

- Work only on the authorized `develop` branch.
- Exactly one input source is configured: local `fasta`, local `genbank`, or `ncbi`.
- NCBI database is fixed to `nuccore` for this issue.
- `NCBI_EMAIL` is required for live acquisition; `NCBI_API_KEY` is optional; neither may appear in manifests, logs, exceptions, or reports.
- Standard tests are fully offline and deterministic; optional live tests live under `network_tests/` and require explicit environment opt-in.
- Preserve exact ordered requested/resolved composition, accession versions, useful GenBank metadata, batch checksums, and consolidated checksum.
- Retry only transient failures with a bounded 1, 2, 4, 8, 16 second schedule capped at 16 seconds.
- Frozen mode performs no network requests and rejects incomplete, corrupt, reordered, or composition-mismatched artifacts.
- Continue using the existing GenBank normalizer and QC; do not duplicate scientific logic.
- Do not implement desktop UI, general step checkpoints, or broad run-state handling in this issue.
- Use strict RED/GREEN TDD and retain command/output evidence for every behavior.
- Do not open a PR, merge `main`, or close or mutate issues.

---

### Task 1: Parse and validate NCBI input configuration

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Modify: `tests/test_config.py`

**Interfaces:**
- Produces: `NcbiInputConfig(query, accessions, frozen_dataset, batch_size, retries, max_records)`
- Produces: `PipelineConfig.input_ncbi: NcbiInputConfig | None`
- Produces: `PipelineConfig.selected_input -> tuple[Path, Literal["fasta", "genbank"]] | NcbiInputConfig`

- [ ] **Step 1: Add RED tests for all three NCBI modes**

Add literal YAML cases and assertions equivalent to:

```python
def test_loads_ncbi_query_configuration(self):
    config = self._load_yaml(
        "target:\n"
        "  name: query-target\n"
        "input:\n"
        "  ncbi:\n"
        "    query: example[Organism]\n"
        "    batch_size: 25\n"
        "    retries: 2\n"
        "    max_records: 50\n"
    )
    self.assertEqual(config.input_ncbi.query, "example[Organism]")
    self.assertEqual(config.input_ncbi.accessions, ())
    self.assertIsNone(config.input_ncbi.frozen_dataset)
    self.assertEqual(config.input_ncbi.batch_size, 25)
    self.assertEqual(config.input_ncbi.retries, 2)
    self.assertEqual(config.input_ncbi.max_records, 50)


def test_loads_ncbi_accessions_and_frozen_modes(self):
    accessions = self._load_yaml(
        "target:\n  name: accessions\ninput:\n  ncbi:\n"
        "    accessions: [NC_000001.11, AB123456.2]\n"
    )
    frozen = self._load_yaml(
        "target:\n  name: frozen\ninput:\n  ncbi:\n"
        "    frozen_dataset: datasets/frozen\n"
    )
    self.assertEqual(accessions.input_ncbi.accessions, ("NC_000001.11", "AB123456.2"))
    self.assertEqual(accessions.input_ncbi.batch_size, 100)
    self.assertEqual(accessions.input_ncbi.retries, 3)
    self.assertEqual(frozen.input_ncbi.frozen_dataset, Path("datasets/frozen"))
```

Use a test-only `_load_yaml(text)` helper that writes a temporary config and returns `load_config`; it must remain in the test class, not production.

- [ ] **Step 2: Add RED validation tests**

Table-drive literal invalid cases for multiple top-level sources; zero or multiple NCBI modes; empty/duplicate/non-string accessions; query-only `max_records`; frozen mode with `batch_size`, `retries`, or `max_records`; `batch_size` outside 1..500; `retries` outside 0..10; and non-positive `max_records`. Each case must assert a stable `ValueError` fragment naming the invalid field or exclusivity rule.

- [ ] **Step 3: Run RED and commit tests**

Run: `python -m unittest tests.test_config -v`

Expected: FAIL because `NcbiInputConfig` and `input.ncbi` parsing do not exist.

Run:

```powershell
git add -- tests/test_config.py
git commit -m "test: define NCBI input configuration contract"
```

- [ ] **Step 4: Implement minimal typed parsing**

Add these types and selection behavior:

```python
@dataclass(frozen=True, slots=True)
class NcbiInputConfig:
    query: str | None = None
    accessions: tuple[str, ...] = ()
    frozen_dataset: Path | None = None
    batch_size: int = 100
    retries: int = 3
    max_records: int | None = None


@dataclass(frozen=True, slots=True)
class PipelineConfig:
    target_name: str
    input_fasta: Path | None = None
    input_genbank: Path | None = None
    input_ncbi: NcbiInputConfig | None = None
    qc: QCConfig = field(default_factory=QCConfig)

    @property
    def selected_input(
        self,
    ) -> tuple[Path, Literal["fasta", "genbank"]] | NcbiInputConfig:
        if self.input_fasta is not None:
            return self.input_fasta, "fasta"
        if self.input_genbank is not None:
            return self.input_genbank, "genbank"
        if self.input_ncbi is not None:
            return self.input_ncbi
        raise ValueError("Exactly one sequence input must be configured.")
```

Implement `_parse_ncbi_input(raw)` with the exact defaults and bounds from the spec. Count configured top-level sources with boolean predicates and require a total of one. In frozen mode, detect whether network-only keys were present in raw YAML rather than accepting their default values silently.

- [ ] **Step 5: Verify GREEN and commit**

Run: `python -m unittest tests.test_config tests.test_minimal_run -v`

Expected: PASS, including unchanged FASTA and GenBank behavior.

Run:

```powershell
git add -- qpcr_pipeline/config.py
git commit -m "feat: configure reproducible NCBI inputs"
```

---

### Task 2: Validate complete frozen datasets

**Files:**
- Create: `qpcr_pipeline/ncbi.py`
- Create: `tests/test_ncbi_acquisition.py`

**Interfaces:**
- Produces: `AcquiredNcbiDataset(records_path: Path, manifest_path: Path)`
- Produces: `validate_frozen_dataset(dataset_dir: str | Path) -> AcquiredNcbiDataset`
- Produces private atomic JSON and SHA-256 helpers reused by acquisition tasks

- [ ] **Step 1: Write a real frozen-dataset fixture helper and RED tests**

In the test file, create real GenBank with `SeqIO.write` and a literal schema-version-1 manifest. Compute the expected checksum only in the fixture builder; test assertions use hand-checked IDs and status. Cover a valid dataset plus each rejected invariant: missing manifest, `PARTIAL`, unsupported schema, missing `records.gb`, wrong byte size, wrong SHA-256, wrong record count, and accession-version order mismatch.

The valid behavior must assert:

```python
dataset = validate_frozen_dataset(dataset_dir)
self.assertEqual(dataset.records_path, dataset_dir / "records.gb")
self.assertEqual(dataset.manifest_path, dataset_dir / "dataset_manifest.json")
with dataset.records_path.open(encoding="utf-8") as handle:
    self.assertEqual(tuple(record.id for record in SeqIO.parse(handle, "genbank")), ("NC_1.2", "NC_2.3"))
```

- [ ] **Step 2: Run RED and commit tests**

Run: `python -m unittest tests.test_ncbi_acquisition.FrozenDatasetTests -v`

Expected: ERROR because `qpcr_pipeline.ncbi` does not exist.

Run:

```powershell
git add -- tests/test_ncbi_acquisition.py
git commit -m "test: define frozen NCBI dataset validation"
```

- [ ] **Step 3: Implement frozen validation**

Create:

```python
MANIFEST_NAME = "dataset_manifest.json"
RECORDS_NAME = "records.gb"
SCHEMA_VERSION = 1


@dataclass(frozen=True, slots=True)
class AcquiredNcbiDataset:
    records_path: Path
    manifest_path: Path


def _sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _read_manifest(path: Path) -> dict[str, object]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(raw, dict):
        raise ValueError("NCBI dataset manifest root must be an object.")
    return raw


def _write_json_atomic(path: Path, value: dict[str, object]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    temporary.replace(path)
```

`validate_frozen_dataset` must validate schema/status, the consolidated metadata fields, bytes/checksum/count, parse through a context-managed GenBank handle, and require parsed record IDs to equal manifest `resolved_entries[*].accession_version` in order. Raise field-specific `ValueError` messages.

- [ ] **Step 4: Verify GREEN and commit**

Run: `python -W default::ResourceWarning -m unittest tests.test_ncbi_acquisition.FrozenDatasetTests -v`

Expected: PASS with no warnings.

Run:

```powershell
git add -- qpcr_pipeline/ncbi.py
git commit -m "feat: validate frozen NCBI datasets"
```

---

### Task 3: Acquire accessions with batches, retries, and resume

**Files:**
- Modify: `qpcr_pipeline/ncbi.py`
- Modify: `tests/test_ncbi_acquisition.py`

**Interfaces:**
- Produces: `NcbiTransientError`, `NcbiFetchedRecord(request_id, record)`, and `NcbiClient` protocol
- Produces: `acquire_ncbi_dataset(config, dataset_dir, *, client=None, sleep=time.sleep, clock=utc_now)`
- Consumes: `NcbiInputConfig` and Task 2 frozen/manifest helpers

- [ ] **Step 1: Add a deterministic fake client and RED accession tests**

The fake implements:

```python
class FakeNcbiClient:
    def __init__(self, records_by_request, failures=()):
        self.records_by_request = records_by_request
        self.failures = list(failures)
        self.fetch_calls = []
        self.resolve_calls = []

    def resolve_query(self, query, max_records):
        self.resolve_calls.append((query, max_records))
        raise AssertionError("query resolution was not expected")

    def fetch_records(self, identifiers, *, identifier_kind):
        requested = tuple(identifiers)
        self.fetch_calls.append((requested, identifier_kind))
        if self.failures:
            failure = self.failures.pop(0)
            if failure is not None:
                raise failure
        return tuple(
            NcbiFetchedRecord(request_id=identifier, record=self.records_by_request[identifier])
            for identifier in reversed(requested)
        )
```

Add grouped RED tests proving:

- five accessions with `batch_size=2` create calls of 2, 2, and 1, but consolidated records and `resolved_entries` follow requested order;
- manifest starts `PARTIAL`, records exact requested accessions, excludes `NCBI_EMAIL`/`NCBI_API_KEY`, and finishes `COMPLETE` with batch and consolidated checksums;
- two `NcbiTransientError` failures retry with injected sleeps `[1, 2]` and then succeed;
- non-transient `ValueError` is attempted once and propagated;
- an interrupted second batch leaves the first batch reusable; a resumed run calls only remaining identifiers;
- a corrupted completed batch is selectively refetched while other valid batches are reused.

- [ ] **Step 2: Run RED and commit tests**

Run: `python -m unittest tests.test_ncbi_acquisition.AccessionAcquisitionTests -v`

Expected: FAIL because acquisition/client interfaces do not exist.

Run:

```powershell
git add -- tests/test_ncbi_acquisition.py
git commit -m "test: define resumable accession acquisition"
```

- [ ] **Step 3: Implement client types and bounded retry**

Add:

```python
class NcbiTransientError(RuntimeError):
    pass


@dataclass(frozen=True, slots=True)
class NcbiFetchedRecord:
    request_id: str
    record: SeqRecord


@dataclass(frozen=True, slots=True)
class ResolvedNcbiQuery:
    uids: tuple[str, ...]
    reported_count: int
    query_translation: str


class NcbiClient(Protocol):
    def resolve_query(self, query: str, max_records: int | None) -> ResolvedNcbiQuery: ...
    def fetch_records(
        self,
        identifiers: tuple[str, ...],
        *,
        identifier_kind: Literal["uid", "accession"],
    ) -> tuple[NcbiFetchedRecord, ...]: ...
```

The protocol method bodies use `...` only as Python protocol declarations, not implementation placeholders.

Implement `_with_retries(operation, retries, sleep)` so attempts equal `retries + 1`, only `NcbiTransientError` is caught, and delays are `min(2 ** retry_index, 16)` for retry indexes starting at zero.

- [ ] **Step 4: Implement accession persistence and resume**

For accession mode, initialize a manifest with this stable shape before fetching:

```python
{
    "schema_version": 1,
    "status": "PARTIAL",
    "source": {"mode": "accessions", "database": "nuccore", "requested_accessions": list(config.accessions)},
    "batch_size": config.batch_size,
    "retries": config.retries,
    "resolved_entries": [],
    "completed_batches": [],
    "consolidated": None,
    "created_at": clock(),
    "updated_at": clock(),
    "tool": "geison-qpcr",
}
```

Implement stable batch planning; validate cached batch checksum/size/count/request IDs before reuse; require the client response to contain exactly one unique `request_id` for every requested identifier; reorder client results to the request order; write GenBank through a temporary file; and update the manifest atomically after each batch.

For every record create a resolved entry with `requested_accession`, `uid: None`, accession without a terminal numeric version suffix, and full `record.id` as `accession_version`. Consolidate batches in global request order, record bytes/SHA/count, mark `COMPLETE`, then call `validate_frozen_dataset` before returning.

- [ ] **Step 5: Verify GREEN and commit**

Run: `python -W default::ResourceWarning -m unittest tests.test_ncbi_acquisition.AccessionAcquisitionTests tests.test_ncbi_acquisition.FrozenDatasetTests -v`

Expected: PASS with no waits and no warnings.

Run:

```powershell
git add -- qpcr_pipeline/ncbi.py
git commit -m "feat: acquire and resume accession datasets"
```

---

### Task 4: Resolve queries and implement the Bio.Entrez adapter

**Files:**
- Modify: `qpcr_pipeline/ncbi.py`
- Modify: `tests/test_ncbi_acquisition.py`

**Interfaces:**
- Produces: query mode in `acquire_ncbi_dataset`
- Produces: `BioEntrezClient.from_environment(environ=os.environ)` implementing `NcbiClient`
- Produces: sanitized `NcbiRequestError` for non-retryable HTTP failures

- [ ] **Step 1: Add RED query-resume tests**

Extend the fake client so configured `ResolvedNcbiQuery` values are returned. Test that query mode:

- persists ordered UIDs and reported count before fetch;
- applies an explicit `max_records` prefix and records selected count/query translation;
- uses UID batches and preserves UID-to-accession-version mapping;
- resumes a partial dataset without resolving the query again;
- rejects reuse when the query, maximum, batch size, or retry limit differs from the partial manifest.

The primary assertion shape is:

```python
self.assertEqual(client.resolve_calls, [("example[Organism]", 2)])
self.assertEqual(client.fetch_calls, [(('101', '102'), "uid")])
self.assertEqual(manifest["source"]["resolved_uids"], ["101", "102"])
self.assertEqual(
    [(entry["uid"], entry["accession_version"]) for entry in manifest["resolved_entries"]],
    [("101", "NC_1.2"), ("102", "NC_2.3")],
)
```

- [ ] **Step 2: Add RED adapter/security tests**

Inject a test Entrez module into `BioEntrezClient`. Its `esearch` returns deterministic dictionaries through `read`; its `efetch` returns a real in-memory GenBank handle. Assert observable adapter results: paged ordered UIDs, reported count/translation, request-to-record mapping, and closed handles. Assert missing/blank `NCBI_EMAIL` fails before the fake module receives a call. Assert no returned value or raised message contains the API-key literal.

- [ ] **Step 3: Run RED and commit tests**

Run: `python -m unittest tests.test_ncbi_acquisition.QueryAcquisitionTests tests.test_ncbi_acquisition.BioEntrezClientTests -v`

Expected: FAIL because query mode and `BioEntrezClient` do not exist.

Run:

```powershell
git add -- tests/test_ncbi_acquisition.py
git commit -m "test: define query and Entrez acquisition contract"
```

- [ ] **Step 4: Implement query mode**

Use the same acquisition loop as accessions. On a new query dataset, call `resolve_query` once, save `resolved_uids`, `reported_count`, `selected_count`, `max_records`, and `query_translation` before fetching, and use identifier kind `"uid"`. On matching partial resume, reuse persisted UIDs without resolving again. A different effective source/config raises `ValueError` rather than mixing artifacts.

- [ ] **Step 5: Implement BioEntrezClient**

Construction follows:

```python
@classmethod
def from_environment(
    cls,
    environ: Mapping[str, str] = os.environ,
) -> "BioEntrezClient":
    email = environ.get("NCBI_EMAIL", "").strip()
    if not email:
        raise ValueError("NCBI_EMAIL must be set for live NCBI acquisition.")
    return cls(email=email, api_key=environ.get("NCBI_API_KEY") or None)
```

Before each operation assign `Entrez.email`, `Entrez.tool = "geison-qpcr"`, and `Entrez.api_key`. Resolve UIDs in stable pages of at most 10,000. Fetch with `Entrez.efetch(db="nuccore", id=comma_joined_ids, rettype="gb", retmode="text")`, parse inside an explicitly closed handle, validate returned count, and map accessions by version/base or UIDs by returned order.

Translate HTTP 429, HTTP 5xx, `URLError`, `TimeoutError`, and connection-level `OSError` into `NcbiTransientError` without embedding credentials or retaining the raw URL-bearing exception as visible context. Translate other HTTP 4xx into `NcbiRequestError("NCBI request failed with HTTP status <code>.") from None`. Propagate parse/composition errors unchanged because they contain local data invariants rather than credential-bearing request URLs.

- [ ] **Step 6: Verify GREEN and commit**

Run: `python -W default::ResourceWarning -m unittest tests.test_ncbi_acquisition -v`

Expected: PASS offline with no real network and no warnings.

Run:

```powershell
git add -- qpcr_pipeline/ncbi.py
git commit -m "feat: resolve NCBI queries through Entrez"
```

---

### Task 5: Route acquired datasets through pipeline QC

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `tests/test_minimal_run.py`
- Create: `network_tests/test_ncbi_live.py`

**Interfaces:**
- Consumes: `PipelineConfig.selected_input`, `NcbiInputConfig`, `acquire_ncbi_dataset`, and `validate_frozen_dataset`
- Produces: `run_pipeline(config, outdir, *, ncbi_client=None) -> RunSummary`
- Produces: `<outdir>/ncbi_dataset_manifest.json` for every NCBI-backed run

- [ ] **Step 1: Add RED offline pipeline tests**

Call `run_pipeline` directly with an accession-mode config and deterministic fake client returning one valid and one invalid nucleotide record. Assert:

```python
self.assertEqual(summary.sequence_ids, ["NC_VALID.1"])
self.assertEqual(qc_report["records"], [
    {"sequence_id": "NC_VALID.1", "status": "ACCEPTED", "reason_codes": []},
    {"sequence_id": "NC_INVALID.1", "status": "REJECTED", "reason_codes": ["INVALID_NUCLEOTIDE"]},
])
self.assertEqual(qc_report["evaluation_set"]["sequence_ids"], ["NC_VALID.1"])
self.assertEqual(effective_manifest["status"], "COMPLETE")
```

Add a frozen-mode test that supplies a client whose every method raises `AssertionError`; assert the run succeeds, the supplied frozen directory remains byte-for-byte unchanged, and the effective manifest is copied into the run output.

- [ ] **Step 2: Run RED and commit tests**

Run: `python -m unittest tests.test_minimal_run -v`

Expected: FAIL because `run_pipeline` still assumes a local `(path, format)` tuple.

Run:

```powershell
git add -- tests/test_minimal_run.py
git commit -m "test: define NCBI pipeline integration"
```

- [ ] **Step 3: Implement pipeline routing**

Create the output directory before source resolution. For local tuple input, preserve existing behavior. For `NcbiInputConfig`:

```python
if selected.frozen_dataset is not None:
    acquired = validate_frozen_dataset(selected.frozen_dataset)
else:
    acquired = acquire_ncbi_dataset(
        selected,
        output_dir / "ncbi_dataset",
        client=ncbi_client,
    )
records = load_genbank(acquired.records_path)
shutil.copyfile(acquired.manifest_path, output_dir / "ncbi_dataset_manifest.json")
```

Then call the existing `evaluate_sequences` once and serialize the existing summary/QC contracts unchanged. Ensure a frozen source manifest is only read and copied, never rewritten.

- [ ] **Step 4: Add the separate opt-in live test**

Create `network_tests/test_ncbi_live.py`. Skip unless `GEISON_RUN_NETWORK_TESTS == "1"` and `NCBI_EMAIL` is non-empty. Acquire the single stable accession `NC_001416.1` into a temporary directory with `batch_size=1`, assert a complete manifest, exact accession version, a parseable non-empty sequence, and successful frozen validation. Do not require or print an API key.

- [ ] **Step 5: Verify GREEN and commit**

Run: `python -W default::ResourceWarning -m unittest tests.test_minimal_run tests.test_config tests.test_ncbi_acquisition -v`

Expected: PASS offline with no warnings.

Run:

```powershell
git add -- qpcr_pipeline/pipeline.py network_tests/test_ncbi_live.py
git commit -m "feat: run QC on reproducible NCBI datasets"
```

---

### Task 6: Completion gate

**Files:**
- Verify only; modify files only through a new focused RED/GREEN cycle if review finds a defect

**Interfaces:**
- Produces: complete offline evidence, reviewed issue diff, and fast-forwarded `origin/develop`

- [ ] **Step 1: Run the standard offline suite**

Run: `python -W default::ResourceWarning -m unittest discover -s tests -p "test_*.py" -v`

Expected: all tests PASS with no network calls and no warnings.

- [ ] **Step 2: Verify optional network separation**

Run without opt-in: `python -m unittest discover -s network_tests -p "test_*.py" -v`

Expected: the live test is skipped and no request is made. Run live only when `GEISON_RUN_NETWORK_TESTS=1` and `NCBI_EMAIL` are supplied by the environment; report that evidence separately.

- [ ] **Step 3: Review the entire issue #3 range**

Review against the written spec and GitHub issue #3. Any Critical or Important finding returns through one focused RED/GREEN fix and scoped re-review.

- [ ] **Step 4: Publish only the authorized branch**

Fetch `origin/develop`, prove it is an ancestor of the local HEAD, then push a normal fast-forward to `origin/develop`. Do not force-push, create a PR, merge `main`, or mutate issues.
