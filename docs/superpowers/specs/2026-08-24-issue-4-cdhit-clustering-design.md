# Issue #4: CD-HIT-EST clustering design

## Goal

Reduce redundant sequences before expensive discovery stages while preserving the
complete QC-approved population for later assay evaluation. The pipeline will expose
an explicit `DiscoverySet` of CD-HIT-EST representatives and leave the
`EvaluationSet` unchanged.

## Scope and non-goals

This issue adds configuration, a CD-HIT-EST process boundary, cluster parsing,
traceable artifacts, and pipeline integration. It does not add alignment, reverse
orientation normalization, checkpoint/resume behavior, or dependency diagnostics;
those belong to issues #5, #11, and #12.

Clustering is disabled by default so the current local and NCBI workflows remain
usable without external scientific binaries. When disabled, the `DiscoverySet`
equals the `EvaluationSet`. When enabled, absence of `cd-hit-est` is an actionable
error before clustering artifacts are published.

## Configuration

`PipelineConfig` gains a frozen `ClusteringConfig`:

```python
@dataclass(frozen=True, slots=True)
class ClusteringConfig:
    enabled: bool = False
    identity: float = 0.95
    threads: int = 1
    memory_mb: int = 800
```

YAML uses an optional top-level section:

```yaml
clustering:
  enabled: true
  identity: 0.95
  threads: 1
  memory_mb: 800
```

Validation is shared by parsed and directly constructed configurations:

- `enabled` is a boolean;
- `identity` is a real number from 0.75 through 1.0, excluding booleans;
- `threads` is an integer from 1 through 256;
- `memory_mb` is a positive integer;
- unknown `clustering` fields are rejected.

The user does not configure CD-HIT word length directly. The pipeline derives a
compatible value from the official CD-HIT-EST ranges:

| Identity | `-n` |
| --- | ---: |
| `0.95 <= c <= 1.0` | 10 |
| `0.90 <= c < 0.95` | 8 |
| `0.88 <= c < 0.90` | 7 |
| `0.85 <= c < 0.88` | 6 |
| `0.80 <= c < 0.85` | 5 |
| `0.75 <= c < 0.80` | 4 |

## Clustering boundary

A new `qpcr_pipeline.clustering` module owns the external-tool integration. It
defines a small runner protocol and a production subprocess implementation so unit
tests can inject a deterministic fake without installing CD-HIT.

The service interface is conceptually:

```python
cluster_sequences(
    records: tuple[LocalSequenceRecord, ...],
    evaluation_set: EvaluationSet,
    config: ClusteringConfig,
    output_dir: Path,
    *,
    runner: CdHitRunner | None = None,
) -> ClusteringResult
```

`ClusteringResult` contains the `DiscoverySet`, ordered clusters, and paths to the
published artifacts. Public result and cluster models are immutable.

The production runner resolves `cd-hit-est` from `PATH` and invokes it with an
argument list, never through a shell. The initial command is equivalent to:

```text
cd-hit-est -i input.fasta -o representatives.fasta
  -c <identity> -n <derived> -d 0 -g 1 -r 0
  -T <threads> -M <memory_mb>
```

`-g 1` requests the more accurate cluster assignment mode. `-r 0` avoids silently
normalizing reverse-complement orientation before issue #5 can record that
transformation. The runner captures output, requires exit code zero, and verifies
that both representative FASTA and `.clstr` files exist.

## Data flow and identity mapping

1. Input acquisition and QC remain unchanged.
2. The pipeline selects the approved records in `EvaluationSet` order.
3. Duplicate approved sequence IDs are rejected because an unambiguous
   original-to-cluster mapping would be impossible.
4. A temporary CD-HIT input FASTA uses stable internal identifiers such as
   `geison-00000000`; an in-memory map preserves original IDs and records.
5. CD-HIT-EST runs only when clustering is enabled and the Evaluation Set is
   non-empty.
6. The `.clstr` parser requires every internal ID exactly once, one representative
   per cluster, no unknown IDs, and representatives present in their own clusters.
7. Final clusters are ordered by the earliest Evaluation Set member. Cluster IDs are
   stable (`cluster-00000`, `cluster-00001`, ...).
8. `DiscoverySet` representatives are ordered by their original Evaluation Set
   positions, independent of raw CD-HIT output order.
9. `EvaluationSet` is never replaced, filtered, or reordered.

The empty Evaluation Set produces an empty Discovery Set without invoking the
external binary.

## Artifacts

Every run publishes atomically:

```text
discovery_set.fasta
clustering_report.json
```

Enabled clustering additionally publishes:

```text
clustering/
  cd-hit-est.clstr
```

`discovery_set.fasta` uses original sequence IDs and sequences. The report uses
schema version 1 and contains:

- status `COMPLETE` and whether clustering was enabled;
- tool name and effective parameters, including derived word length;
- ordered Evaluation and Discovery Set IDs;
- counts before and after clustering;
- ordered clusters with cluster ID, representative ID, and ordered original member
  IDs;
- each non-representative member's reported identity and strand when available;
- relative artifact filenames.

When clustering is disabled, the report records the default parameters,
`DiscoverySet == EvaluationSet`, and an empty cluster list. The QC report gains a
`discovery_set.sequence_ids` field. Existing summary fields and Evaluation Set
semantics remain unchanged.

All final files use temporary siblings plus atomic replacement. Raw temporary input,
internal-ID representative FASTA, stdout, and stderr are not published. The raw
`.clstr` file is retained for audit, while the JSON report is the human-readable
original-ID mapping.

## Error handling

- Missing `cd-hit-est` raises a `CdHitError` naming the missing executable and
  explaining that clustering can be disabled or the dependency installed.
- Nonzero exit raises a sanitized `CdHitError` with the exit code and bounded stderr;
  no shell command string is evaluated.
- Missing output, malformed `.clstr`, duplicate membership, missing members, unknown
  internal IDs, or inconsistent representatives fail before final artifacts are
  published.
- Invalid configuration fails before output mutation or subprocess execution.
- A clustering failure does not modify the Evaluation Set or produce a `COMPLETE`
  clustering report.

## Testing

The standard suite remains offline and deterministic:

- configuration tests cover defaults, valid YAML, unknown fields, direct configs,
  bounds, booleans, and word-length boundaries;
- parser tests use literal realistic `.clstr` text and cover malformed composition;
- service tests use a fake runner that writes real FASTA and `.clstr` files, proving
  representative ordering, original-ID restoration, atomic artifacts, missing-tool
  diagnostics, and empty/disabled behavior;
- pipeline tests prove clustering consumes only QC-approved records, the
  `EvaluationSet` is unchanged, the `DiscoverySet` is reduced, mappings are
  persisted, and existing local/NCBI runs remain compatible.

A real integration test lives under `integration_tests/`, outside default unittest
and pytest discovery. It runs only when `cd-hit-est` is present, uses a small
deterministic FASTA fixture, and verifies the representative and mapping artifacts.
No fallback algorithm is implemented because a Python substitute would not preserve
CD-HIT-EST scientific behavior.

## Acceptance mapping

- CD-HIT-EST integration: isolated production subprocess runner.
- Configurable identity: validated `identity` with automatically compatible `-n`.
- Discovery representatives: explicit `DiscoverySet` and `discovery_set.fasta`.
- Original to cluster to representative mapping: `clustering_report.json`.
- Evaluation Set preserved: immutable pipeline contract and integration tests.
- Legible persisted artifacts: original-ID FASTA/JSON plus raw `.clstr` audit file.
- Separation tests: deterministic fake-runner suite and optional real-tool test.
