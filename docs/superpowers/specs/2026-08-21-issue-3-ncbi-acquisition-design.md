# Issue 3: Reproducible NCBI Acquisition Design

## Objective

Allow Geison to build the same validated Target Sequence Set from either an NCBI nucleotide query or an explicit accession list. Downloads must be batched, bounded in retry behavior, resumable after interruption, frozen for offline reuse, and traceable down to accession versions and source metadata.

This design implements GitHub issue `BrunoDCamargo/Geison#3` on `develop`. It does not add a desktop UI, general pipeline checkpoints, or broad run-state management; those belong to later issues.

## Design choice

Use Biopython `Bio.Entrez`, which is already a project dependency, rather than introducing NCBI Datasets CLI or maintaining direct E-utilities HTTP and XML parsing. Network access remains behind a small client interface so the standard suite can exercise acquisition deterministically without contacting NCBI.

## Configuration

Existing local inputs remain valid. Exactly one top-level input source is allowed: `fasta`, `genbank`, or `ncbi`.

Query mode:

```yaml
target:
  name: example-target
input:
  ncbi:
    query: '"Example virus"[Organism] AND complete genome[Title]'
    batch_size: 100
    retries: 3
    max_records: 500
```

Explicit accession mode:

```yaml
target:
  name: example-target
input:
  ncbi:
    accessions:
      - NC_000001.11
      - AB123456.2
    batch_size: 100
    retries: 3
```

Frozen reuse mode:

```yaml
target:
  name: example-target
input:
  ncbi:
    frozen_dataset: datasets/example-target
```

Rules:

- Exactly one of `query`, `accessions`, or `frozen_dataset` is required inside `input.ncbi`.
- `accessions` must be a non-empty ordered list of non-empty unique strings. Versioned and unversioned accessions are accepted; the resolved accession version is always recorded.
- `batch_size` is an integer from 1 through 500 and defaults to 100.
- `retries` is an integer from 0 through 10 and defaults to 3. It is the number of retries after the initial attempt.
- `max_records` is an optional positive integer valid only for query mode. When omitted, the complete query result is selected. When set, truncation is explicit and recorded in the manifest.
- `frozen_dataset` cannot specify network options.
- The NCBI database is fixed to `nuccore` for this issue.
- `NCBI_EMAIL` is required for live acquisition. `NCBI_API_KEY` is optional. Both are read from environment variables only and never persisted or logged.

## Components

### Configuration model

`qpcr_pipeline.config` gains a frozen `NcbiInputConfig` and extends `PipelineConfig` with `input_ncbi`. Local input compatibility is preserved. The pipeline exposes a selected input descriptor rather than assuming every source is a local path.

Configuration validation rejects ambiguous input sources, empty queries, duplicate or empty accessions, invalid numeric bounds, and frozen/network option mixtures before any network request.

### NCBI client boundary

`qpcr_pipeline.ncbi` contains a small `NcbiClient` interface and a `BioEntrezClient` implementation. The interface provides two operations:

1. Resolve a query into an ordered list of NCBI nucleotide UIDs, preserving the reported total count and query translation.
2. Fetch a requested batch of UIDs or accessions as GenBank records.

The production adapter sets `Entrez.email`, `Entrez.tool`, and the optional API key immediately before each request. It does not expose credentials through exceptions, manifests, logs, or return values.

Tests use a deterministic fake client at this boundary. They do not mock or assert Biopython internals.

### Acquisition service

`acquire_ncbi_dataset(config, dataset_dir, client=None)` orchestrates acquisition and returns the path to a complete frozen GenBank dataset plus its manifest.

For query mode:

1. Resolve the complete ordered UID composition, or the explicit `max_records` prefix.
2. Persist that composition before record downloads begin.
3. Fetch unresolved UIDs in stable batches.

For accession mode:

1. Persist the exact ordered requested list before downloads begin.
2. Fetch unresolved accessions in stable batches.
3. Match returned records to requests and preserve the requested ordering in the consolidated dataset.

Every fetched record preserves its NCBI UID when available, requested accession when applicable, accession, accession version, and normalized GenBank metadata.

### Retry behavior

Only transient transport failures are retried: HTTP 429, HTTP 5xx, URL/network errors, timeouts, and connection-level `OSError` failures. Parse errors, invalid responses, missing requested records, configuration errors, and HTTP 4xx other than 429 fail immediately.

Retry delay is bounded exponential backoff: 1, 2, 4, 8, then 16 seconds, capped at 16 seconds. The client receives injectable sleep and clock functions so offline tests do not wait and manifests remain deterministic.

### Resumable persistence

The dataset directory contains:

```text
dataset_manifest.json
batches/
  batch-00000.gb
  batch-00001.gb
records.gb
```

Each batch is written to a temporary sibling and atomically renamed only after it parses successfully and contains the expected records. The manifest is atomically rewritten after each completed batch.

On resume, the service validates every completed batch against its stored SHA-256, byte size, record count, and requested identifiers. Valid batches are reused. A missing, corrupt, or mismatched batch is downloaded again without discarding other valid batches.

After all batches pass validation, the service concatenates records in the frozen composition order into `records.gb`, validates the consolidated file, records its SHA-256, and marks the manifest `COMPLETE`. An interrupted or failed acquisition leaves a `PARTIAL` manifest and completed batches for the next attempt.

### Manifest

`dataset_manifest.json` uses `schema_version: 1` and contains:

- `status`: `PARTIAL` or `COMPLETE`;
- source mode and database;
- original query or ordered requested accessions;
- NCBI reported count, explicit maximum, selected count, and query translation when applicable;
- resolved ordered entries with UID when available, requested accession when applicable, accession, and accession version;
- batch size and retry limit;
- completed batch filenames, requested identifiers, record count, byte size, and SHA-256;
- consolidated GenBank filename, record count, byte size, and SHA-256 when complete;
- Geison tool identifier and retrieval timestamps supplied by an injectable UTC clock.

Secrets are structurally absent. The manifest records effective dataset composition, not environment credentials.

### Frozen mode

Frozen mode performs no network operation. It requires a `COMPLETE` schema-version-1 manifest and validates the consolidated file checksum, size, record count, accession versions, and ordering before loading it.

Any validation mismatch raises an actionable error instead of silently accepting a changed dataset.

### Pipeline integration

For NCBI query or accession input, `run_pipeline` acquires or resumes the dataset under `<outdir>/ncbi_dataset`, then loads its `records.gb` through the existing GenBank normalizer and runs the existing QC unchanged.

For frozen input, the pipeline validates and loads the supplied dataset directory. It copies the effective acquisition manifest to `<outdir>/ncbi_dataset_manifest.json` so the run remains auditable without modifying the frozen source.

`run_summary.json` and `qc_report.json` retain their current contracts. Acquisition failures raise a clear exception after persisting partial state; general `PARTIAL` and `FAILED` run states are deferred to issue #12.

## Error handling

- Configuration errors fail before network or output mutation.
- Missing `NCBI_EMAIL` fails before the first live request with an actionable message.
- A missing accession or unexpected record composition fails the affected acquisition rather than silently dropping records.
- Exhausted retries preserve the underlying error type as the cause and report the operation and attempt count without credentials or full sequences.
- Frozen validation reports the exact artifact or invariant that failed.

## Testing strategy

The standard offline suite adds:

- configuration tests for query, accessions, frozen mode, mutual exclusion, bounds, duplicate accessions, and environment-independent parsing;
- acquisition tests with a real deterministic fake client for stable batching, bounded retry, non-retryable failures, partial persistence, resume reuse, selective corrupt-batch refetch, ordering, accession versions, manifest contents, and absence of secrets;
- frozen-mode tests for zero client calls and checksum/composition validation;
- pipeline integration tests proving NCBI records flow through the existing QC into all approved Evaluation Set IDs and traceable reason codes.

Optional live network tests live outside the standard `tests/` discovery tree under `network_tests/`. They run only when `GEISON_RUN_NETWORK_TESTS=1` and `NCBI_EMAIL` are present, use a small stable accession fixture, and never require an API key.

The completion gate runs the entire offline suite with `ResourceWarning` visible. The live test is reported separately and is never required for ordinary offline development.

## Acceptance-criteria mapping

- Query and accessions: configuration, client boundary, and acquisition service.
- Batches and bounded retries: stable batch planner and retry policy.
- API key secrecy: environment-only adapter configuration and secret-free artifacts/tests.
- Accession and accession version: ordered resolved manifest entries and GenBank validation.
- Metadata: existing GenBank normalizer receives fetched records unchanged.
- Frozen reuse: complete manifest plus checksummed consolidated GenBank, validated without network.
- Network-test separation: opt-in `network_tests/`, excluded from the standard offline suite.
