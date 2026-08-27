# Issue #9 Specificity Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add deterministic, offline off-target specificity evaluation for Primer3 assays, distinguishing isolated oligo hits, plausible primer amplicons, and probe-detectable off-targets across FASTA and frozen NCBI datasets.

**Architecture:** Add typed off-target/specificity configuration, isolate dataset loading and provenance in `off_targets.py`, keep public models and orchestration in `specificity.py`, and put exhaustive IUPAC matching plus F/R geometry in the focused internal module `specificity_matching.py`. Integrate after inclusivity in `run_pipeline()` without requiring inclusivity to be enabled, publish three auditable artifacts, and keep all standard tests offline.

**Tech Stack:** Python >=3.10, existing Biopython >=1.85,<2 and PyYAML >=6,<7, frozen/slotted dataclasses, standard library (`hashlib`, `json`, `pathlib`, `typing`, `uuid`), unittest-style tests collected by pytest. No new dependency, BLAST executable, or network call.

**Spec:** `docs/superpowers/specs/2026-08-27-issue-9-specificity-design.md`

## Global Constraints

- Work on `feature/issue-9-specificity`; do not change `main` while implementing this issue.
- `specificity.enabled` defaults to `false` and existing configurations must retain current behavior.
- Enabled specificity requires enabled primer design and at least one `off_targets` dataset.
- Each off-target has one unique nonblank `name` and exactly one source: local `fasta` or `frozen_dataset`.
- Specificity never performs live NCBI acquisition; frozen NCBI data must pass existing `validate_frozen_dataset()`.
- Search every off-target sequence and its reverse complement using the existing conservative IUPAC target-support-subset semantics.
- Mismatch positions are 1-based in oligo synthesis orientation; final-primer-3-prime rejection is independently configurable for specificity.
- Find all compatible hits needed for geometry before applying `max_hits_per_oligo_per_dataset`; truncation may reduce only the published individual-hit artifact and must never reduce scientific risk.
- `primer_amplicon_plausible` requires compatible F/R sites facing inward on the same oriented sequence with amplicon size <= `max_amplicon_size`.
- `detectable_off_target` additionally requires a compatible probe site within that F/R interval.
- Preserve multiple plausible amplicons rather than silently choosing one.
- Search and tests are deterministic and offline; BLAST, `blastn-short`, BLAST databases, and large-bank indexing are out of scope.
- Do not call private helpers from `inclusivity.py`; a semantic-equivalence test is allowed to compare public observable behavior.
- Follow RED -> GREEN -> REFACTOR TDD for production behavior changes.
- CircleCI runs automatically only after merge to `develop`; feature-branch work should not consume executor credits.

## File map

- Modify `qpcr_pipeline/config.py`: define `OffTargetConfig` and `SpecificityConfig`, parse `off_targets` and `specificity`, and enforce dependencies/validation.
- Create `qpcr_pipeline/off_targets.py`: validate/load configured FASTA or frozen NCBI data, compute provenance, and expose immutable loaded datasets.
- Create `qpcr_pipeline/specificity_matching.py`: exhaustive IUPAC-compatible hit enumeration, compatibility classification, deterministic ordering, and complete F/R/probe geometry.
- Create `qpcr_pipeline/specificity.py`: public result models, stage orchestration, retention accounting, TSV/JSON rendering, stale cleanup, and atomic publication.
- Modify `qpcr_pipeline/pipeline.py`: call specificity after inclusivity and add the top-level QC summary.
- Create `tests/test_specificity_config.py`: defaults, YAML parsing, invalid types/sources, duplicate names, and dependencies.
- Create `tests/test_specificity_matching.py`: exact/mismatch/IUPAC/orientation ordering and full-hit retention semantics.
- Create `tests/test_specificity_geometry.py`: isolated hits, invalid orientation/distance, probe placement, reverse-complement cases, and multiple amplicons.
- Create `tests/test_specificity_artifacts.py`: FASTA/frozen loading, provenance, empty cases, truncation accounting, artifact content, stale cleanup, and atomic failure behavior.
- Create `tests/test_pipeline_specificity.py`: disabled and enabled pipeline routing, no-network guarantee, QC summary, and IUPAC semantic-equivalence regression.
- Modify `README.md`: document effective configuration, artifacts, risk semantics, offline behavior, and scaling limit.

---

### Task 1: Configuration contract for off-target datasets and specificity

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Create: `tests/test_specificity_config.py`

**Interfaces:**
- Consumes: existing `PipelineConfig`, `PrimerDesignConfig`, `load_config()`, `validate_pipeline_config()`.
- Produces:
  - `OffTargetConfig(name: str, fasta: Path | None = None, frozen_dataset: Path | None = None)`
  - `SpecificityConfig(enabled: bool = False, max_hits_per_oligo_per_dataset: int = 20, max_primer_mismatches: int = 2, max_probe_mismatches: int = 1, reject_primer_3_prime_mismatch: bool = True, primer_3_prime_bases: int = 5, max_amplicon_size: int = 1000)`
  - `PipelineConfig.off_targets: tuple[OffTargetConfig, ...] = ()`
  - `PipelineConfig.specificity: SpecificityConfig`
  - `validate_off_target_config()` and `validate_specificity_config()`.

- [ ] **Step 1: Write failing default and YAML parsing tests**

Create `tests/test_specificity_config.py` using the existing temporary-YAML pattern:

```python
import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import (
    OffTargetConfig,
    PipelineConfig,
    SpecificityConfig,
    load_config,
)

FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class SpecificityConfigTests(unittest.TestCase):
    def _load_yaml(self, text: str) -> PipelineConfig:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "config.yaml"
            path.write_text(text, encoding="utf-8")
            return load_config(path)

    def test_defaults_are_disabled_and_have_no_off_targets(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.off_targets, ())
        self.assertEqual(config.specificity, SpecificityConfig())

    def test_loads_fasta_and_frozen_off_targets_with_specificity(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
            "primer_design:\n  enabled: true\n"
            "off_targets:\n"
            "  - name: human\n    fasta: data/human.fasta\n"
            "  - name: neighbors\n    frozen_dataset: data/neighbors\n"
            "specificity:\n"
            "  enabled: true\n"
            "  max_hits_per_oligo_per_dataset: 7\n"
            "  max_primer_mismatches: 1\n"
            "  max_probe_mismatches: 0\n"
            "  reject_primer_3_prime_mismatch: false\n"
            "  primer_3_prime_bases: 4\n"
            "  max_amplicon_size: 600\n"
        )
        self.assertEqual(
            config.off_targets,
            (
                OffTargetConfig(name="human", fasta=Path("data/human.fasta")),
                OffTargetConfig(name="neighbors", frozen_dataset=Path("data/neighbors")),
            ),
        )
        self.assertEqual(
            config.specificity,
            SpecificityConfig(
                enabled=True,
                max_hits_per_oligo_per_dataset=7,
                max_primer_mismatches=1,
                max_probe_mismatches=0,
                reject_primer_3_prime_mismatch=False,
                primer_3_prime_bases=4,
                max_amplicon_size=600,
            ),
        )
```

- [ ] **Step 2: Run focused tests and verify RED**

Run:

```bash
python -m unittest tests.test_specificity_config.SpecificityConfigTests.test_defaults_are_disabled_and_have_no_off_targets -v
```

Expected: import failure for `OffTargetConfig` or `SpecificityConfig`.

- [ ] **Step 3: Add immutable config models and parser wiring**

Add exact models near the other stage config dataclasses:

```python
@dataclass(frozen=True, slots=True)
class OffTargetConfig:
    name: str
    fasta: Path | None = None
    frozen_dataset: Path | None = None


@dataclass(frozen=True, slots=True)
class SpecificityConfig:
    enabled: bool = False
    max_hits_per_oligo_per_dataset: int = 20
    max_primer_mismatches: int = 2
    max_probe_mismatches: int = 1
    reject_primer_3_prime_mismatch: bool = True
    primer_3_prime_bases: int = 5
    max_amplicon_size: int = 1000
```

Add `off_targets` and `specificity` fields to `PipelineConfig`. Parse `off_targets` only as a list of mappings; reject unknown keys; parse paths as nonblank strings. Parse `specificity` only as a mapping and reject unknown keys.

- [ ] **Step 4: Write failing validation/dependency tests**

Cover: `off_targets` not a list, entry not mapping, blank name, duplicate name, both sources, no source, non-path source values, specificity unknown field, boolean passed as integer, negative mismatch limits, zero positive fields, enabled without primer design, and enabled with zero off-targets.

Use direct-construction dependency assertions:

```python
with self.assertRaisesRegex(ValueError, "specificity.*requires enabled primer design"):
    PipelineConfig(
        target_name="target",
        input_fasta=FIXTURE_FASTA,
        off_targets=(OffTargetConfig(name="human", fasta=Path("human.fasta")),),
        specificity=SpecificityConfig(enabled=True),
    ).selected_input

with self.assertRaisesRegex(ValueError, "specificity.*at least one off-target"):
    PipelineConfig(
        target_name="target",
        input_fasta=FIXTURE_FASTA,
        primer_design=__import__("qpcr_pipeline.config", fromlist=["PrimerDesignConfig"]).PrimerDesignConfig(enabled=True),
        specificity=SpecificityConfig(enabled=True),
    ).selected_input
```

- [ ] **Step 5: Implement exact validation**

Implement:

```python
def validate_off_target_config(config: OffTargetConfig) -> None:
    if not isinstance(config, OffTargetConfig):
        raise ValueError("Off-target configuration must be an OffTargetConfig.")
    if not isinstance(config.name, str) or not config.name.strip():
        raise ValueError("Off-target name must be a non-blank string.")
    for field_name, value in (("fasta", config.fasta), ("frozen_dataset", config.frozen_dataset)):
        if value is not None and not isinstance(value, Path):
            raise ValueError(f"Off-target {field_name} must be a Path when configured.")
    if sum(value is not None for value in (config.fasta, config.frozen_dataset)) != 1:
        raise ValueError("Off-target must configure exactly one of fasta or frozen_dataset.")


def validate_specificity_config(config: SpecificityConfig) -> None:
    if not isinstance(config, SpecificityConfig):
        raise ValueError("Specificity configuration must be a SpecificityConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Specificity enabled must be a boolean.")
    if not isinstance(config.reject_primer_3_prime_mismatch, bool):
        raise ValueError("Specificity reject_primer_3_prime_mismatch must be a boolean.")
    for name in ("max_primer_mismatches", "max_probe_mismatches"):
        value = getattr(config, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 0:
            raise ValueError(f"Specificity {name} must be a non-negative integer.")
    for name in ("max_hits_per_oligo_per_dataset", "primer_3_prime_bases", "max_amplicon_size"):
        value = getattr(config, name)
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"Specificity {name} must be a positive integer.")
```

In `validate_pipeline_config()`, validate all off-targets, reject duplicate names, validate specificity, then enforce enabled dependencies.

- [ ] **Step 6: Run config tests and commit**

Run:

```bash
python -m unittest tests.test_specificity_config tests.test_config -q
```

Expected: all configuration tests pass.

Commit:

```bash
git add qpcr_pipeline/config.py tests/test_specificity_config.py
git commit -m "feat: configure off-target specificity"
```

---

### Task 2: Deterministic off-target dataset loading and provenance

**Files:**
- Create: `qpcr_pipeline/off_targets.py`
- Create: `tests/test_specificity_artifacts.py`

**Interfaces:**
- Consumes: `OffTargetConfig`, `LocalSequenceRecord`, `load_fasta()`, `load_genbank()`, `validate_frozen_dataset()`.
- Produces:
  - `OffTargetDataset(name, source_type, source_path, sha256, sequence_ids, records, frozen_manifest_path, frozen_manifest)`
  - `load_off_target_dataset(config: OffTargetConfig) -> OffTargetDataset`
  - `load_off_target_datasets(configs: tuple[OffTargetConfig, ...]) -> tuple[OffTargetDataset, ...]`.

- [ ] **Step 1: Write failing FASTA provenance tests**

```python
class OffTargetDatasetTests(unittest.TestCase):
    def test_loads_fasta_with_stable_sha_and_sequence_ids(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "off.fasta"
            path.write_text(">off-1\nACGT\n>off-2\nTTTT\n", encoding="utf-8")
            dataset = load_off_target_dataset(
                OffTargetConfig(name="human", fasta=path)
            )
        self.assertEqual(dataset.name, "human")
        self.assertEqual(dataset.source_type, "FASTA")
        self.assertEqual(dataset.sequence_ids, ("off-1", "off-2"))
        self.assertEqual(len(dataset.sha256), 64)
        self.assertIsNone(dataset.frozen_manifest_path)
```

- [ ] **Step 2: Verify RED**

Run:

```bash
python -m unittest tests.test_specificity_artifacts.OffTargetDatasetTests.test_loads_fasta_with_stable_sha_and_sequence_ids -v
```

Expected: import failure for `qpcr_pipeline.off_targets`.

- [ ] **Step 3: Implement FASTA loading, file checksum, and duplicate-ID validation**

Define:

```python
@dataclass(frozen=True, slots=True)
class OffTargetDataset:
    name: str
    source_type: Literal["FASTA", "NCBI_FROZEN"]
    source_path: Path
    sha256: str
    sequence_ids: tuple[str, ...]
    records: tuple[LocalSequenceRecord, ...]
    frozen_manifest_path: Path | None
    frozen_manifest: dict[str, object] | None
```

Use `hashlib.sha256(path.read_bytes()).hexdigest()`. Require every record ID to be nonblank and unique. Empty FASTA is valid and yields empty records/IDs. Wrap file/parse errors as `ValueError` identifying dataset name/path without dumping sequence contents.

- [ ] **Step 4: Write failing frozen-dataset validation/provenance tests**

Build frozen data using the existing deterministic `acquire_ncbi_dataset()` test pattern, then assert `load_off_target_dataset(OffTargetConfig(name="neighbors", frozen_dataset=dir))` calls the existing validator, loads `records.gb`, reads `dataset_manifest.json`, exposes `source_type == "NCBI_FROZEN"`, and preserves the manifest `source` plus `resolved_entries` in `frozen_manifest`.

Add a corrupt-manifest case and assert a contextual `ValueError` mentioning off-target `neighbors`.

- [ ] **Step 5: Implement frozen loading without network access**

Use only:

```python
acquired = validate_frozen_dataset(config.frozen_dataset)
records = load_genbank(acquired.records_path)
manifest = json.loads(acquired.manifest_path.read_text(encoding="utf-8"))
```

Do not import or call `acquire_ncbi_dataset`, `BioEntrezClient`, `Entrez`, or any client boundary from production specificity loading.

- [ ] **Step 6: Add order and invalid-data tests, run, commit**

Test configured dataset order is preserved, duplicate record IDs are rejected, missing FASTA fails, empty FASTA succeeds, and malformed FASTA content containing invalid sequence symbols is left for matching to reject contextually rather than silently normalized away.

Run:

```bash
python -m unittest tests.test_specificity_artifacts.OffTargetDatasetTests -q
```

Commit:

```bash
git add qpcr_pipeline/off_targets.py tests/test_specificity_artifacts.py
git commit -m "feat: load reproducible off-target datasets"
```

---

### Task 3: Exhaustive IUPAC-compatible oligo matching

**Files:**
- Create: `qpcr_pipeline/specificity_matching.py`
- Create: `tests/test_specificity_matching.py`

**Interfaces:**
- Consumes: `AssayCandidate`, `DesignedOligo`, `SpecificityConfig`, `LocalSequenceRecord`, and public helpers from `qpcr_pipeline.iupac`.
- Produces:
  - aliases `OligoRole = Literal["FORWARD", "PROBE", "REVERSE"]`, `Orientation = Literal["FORWARD", "REVERSE_COMPLEMENT"]`
  - `MatchHit` internal immutable model
  - `enumerate_compatible_hits(dataset_name, record, assay, role, config) -> tuple[MatchHit, ...]`
  - `all_assay_hits(dataset_name, records, assays, config) -> tuple[MatchHit, ...]`.

- [ ] **Step 1: Write failing exact, mismatch, and 3-prime tests**

Use fixed `DesignedOligo`/`AssayCandidate` builders and assert exact sites plus mismatch semantics:

```python
hits = enumerate_compatible_hits(
    "human",
    LocalSequenceRecord("off-1", "TTACGTACTT"),
    self._assay(forward="ACGT"),
    "FORWARD",
    SpecificityConfig(max_primer_mismatches=1, primer_3_prime_bases=2),
)
self.assertIn(
    (3, 6, (), True),
    [(h.source_start, h.source_end, h.mismatch_positions, h.compatible) for h in hits],
)
```

Add one internal mismatch allowed and one final-3-prime mismatch rejected.

- [ ] **Step 2: Verify RED**

Run:

```bash
python -m unittest tests.test_specificity_matching.SpecificityMatchingTests -v
```

Expected: import failure for `specificity_matching`.

- [ ] **Step 3: Implement internal match model and compatibility rule**

Use:

```python
@dataclass(frozen=True, slots=True)
class MatchHit:
    dataset_name: str
    assay_id: str
    sequence_id: str
    role: OligoRole
    orientation: Orientation
    source_start: int
    source_end: int
    oriented_start: int
    oriented_end: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    compatible: bool
```

Normalize every sequence/oligo through `normalize_iupac()`. For every position of both the original record and `reverse_complement_iupac(record.sequence)`, extract a segment of oligo length. For reverse-primer comparison, reverse-complement the segment into primer synthesis orientation before `mismatch_positions()`; forward/probe compare directly. Convert reverse-complement oriented coordinates back to source coordinates with:

```python
source_start = record_length - oriented_end + 1
source_end = record_length - oriented_start + 1
```

Retain only compatible hits from this API because specificity geometry is defined from compatible sites; tests should verify incompatible mismatch sites are not returned.

- [ ] **Step 4: Write failing IUPAC, reverse-complement, and stable-order tests**

Cover:
- oligo `R` matching target `A`/`G` but conservative target `N` not matching `R`;
- degenerate oligo evaluated directly without expansion;
- whole record reverse-complement hit;
- forward/probe/reverse role behavior;
- two equal hits sorted by orientation rank then coordinate;
- multiple records and assays sorted deterministically.

Require global order key:

```python
(
    hit.dataset_name,
    hit.assay_id,
    hit.sequence_id,
    {"FORWARD": 0, "PROBE": 1, "REVERSE": 2}[hit.role],
    0 if hit.orientation == "FORWARD" else 1,
    hit.source_start,
    hit.source_end,
    hit.mismatch_count,
    hit.mismatch_positions,
)
```

- [ ] **Step 5: Implement `all_assay_hits()` and contextual IUPAC errors**

Loop dataset record order, assay order, role order, but sort the final tuple by the documented key. Raise `SpecificityError` from `specificity.py` only through a dependency-neutral local `SpecificityMatchingError` here; the public orchestrator will translate it. Error text must identify dataset, assay, sequence and role while not including the whole sequence.

- [ ] **Step 6: Run matching tests and commit**

Run:

```bash
python -m unittest tests.test_specificity_matching -q
```

Commit:

```bash
git add qpcr_pipeline/specificity_matching.py tests/test_specificity_matching.py
git commit -m "feat: find deterministic off-target oligo hits"
```

---

### Task 4: F/R amplicon geometry and probe detectability

**Files:**
- Modify: `qpcr_pipeline/specificity_matching.py`
- Create: `tests/test_specificity_geometry.py`

**Interfaces:**
- Consumes: complete, untruncated `MatchHit` tuples for one dataset.
- Produces:
  - `GeometryAmplicon` internal immutable model
  - `find_plausible_amplicons(hits, config) -> tuple[GeometryAmplicon, ...]`.

- [ ] **Step 1: Write failing isolated-hit and invalid-geometry tests**

Create hand-built `MatchHit` fixtures rather than invoking search. Assert a single forward hit, F/R on different sequence IDs, opposite orientation values, reverse before forward, and size > `max_amplicon_size` all produce no amplicons.

- [ ] **Step 2: Verify RED**

Run:

```bash
python -m unittest tests.test_specificity_geometry -v
```

Expected: import failure for `GeometryAmplicon` or `find_plausible_amplicons`.

- [ ] **Step 3: Implement inward-facing geometry in oriented coordinates**

Use:

```python
@dataclass(frozen=True, slots=True)
class GeometryAmplicon:
    dataset_name: str
    assay_id: str
    sequence_id: str
    orientation: Orientation
    forward: MatchHit
    reverse: MatchHit
    probes: tuple[MatchHit, ...]
    source_start: int
    source_end: int
    amplicon_size: int
    primer_amplicon_plausible: bool
    detectable_off_target: bool
```

Group by `(dataset_name, assay_id, sequence_id, orientation)`. For each forward/reverse pair require `forward.oriented_end < reverse.oriented_start`, then calculate `amplicon_size = reverse.oriented_end - forward.oriented_start + 1` and require `<= max_amplicon_size`. A compatible probe counts only when `forward.oriented_end < probe.oriented_start <= probe.oriented_end < reverse.oriented_start`.

Convert the overall oriented interval back to source coordinates for reverse-complement orientation; source start must always be <= source end.

- [ ] **Step 4: Write failing positive/probe/multi-amplicon tests**

Cover:
- plausible F/R without probe: plausible true, detectable false;
- probe outside interval: detectable false;
- probe inside interval: detectable true;
- same geometry on reverse-complement orientation;
- two forward and two reverse hits yielding multiple valid combinations; preserve every combination;
- exact `max_amplicon_size` boundary is accepted.

Assert deterministic amplicon ordering by dataset, assay, sequence, orientation rank, source_start, source_end, then primer hit coordinates.

- [ ] **Step 5: Run geometry tests and commit**

Run:

```bash
python -m unittest tests.test_specificity_geometry -q
```

Commit:

```bash
git add qpcr_pipeline/specificity_matching.py tests/test_specificity_geometry.py
git commit -m "feat: detect plausible off-target amplicons"
```

---

### Task 5: Public specificity service, truncation accounting, and atomic artifacts

**Files:**
- Create: `qpcr_pipeline/specificity.py`
- Modify: `tests/test_specificity_artifacts.py`

**Interfaces:**
- Consumes: `PrimerDesignResult`, `SpecificityConfig`, configured `OffTargetConfig` tuple, `load_off_target_datasets()`, complete `MatchHit` tuples, and geometry results.
- Produces:
  - `SpecificityError`
  - public `OffTargetHit`, `PlausibleAmplicon`, `SpecificityResult`, `HitRetentionSummary`
  - `evaluate_specificity(primer_design, off_target_configs, config, output_dir) -> SpecificityResult`
  - `specificity/off_target_hits.tsv`, `specificity/plausible_amplicons.tsv`, `specificity/specificity_report.json`.

- [ ] **Step 1: Write failing public-model, disabled, and empty-assay tests**

Define expected public fields with `dataclasses.fields()`. Disabled mode must not inspect deliberately invalid off-target paths and must publish only `specificity_report.json` with status `SKIPPED`. Enabled mode with a COMPLETE PrimerDesignResult containing zero assays must still validate/load configured datasets and publish both header-only TSVs plus COMPLETE JSON.

Use this disabled assertion:

```python
result = evaluate_specificity(
    self._primer_result(status="SKIPPED"),
    (OffTargetConfig(name="bad", fasta=Path("does-not-exist.fasta")),),
    SpecificityConfig(enabled=False),
    output_dir,
)
self.assertEqual(result.status, "SKIPPED")
self.assertEqual(
    {p.name for p in (output_dir / "specificity").iterdir()},
    {"specificity_report.json"},
)
```

- [ ] **Step 2: Verify RED**

Run:

```bash
python -m unittest tests.test_specificity_artifacts.SpecificityArtifactTests.test_disabled_does_not_read_datasets -v
```

Expected: import failure for `evaluate_specificity`.

- [ ] **Step 3: Implement public immutable models and conversion from internal matches**

Use public models that contain no raw sequence strings:

```python
@dataclass(frozen=True, slots=True)
class OffTargetHit:
    dataset_name: str
    assay_id: str
    sequence_id: str
    role: str
    orientation: str
    hit_rank: int
    source_start: int
    source_end: int
    mismatch_positions: tuple[int, ...]
    mismatch_count: int
    exact_match: bool
    three_prime_mismatch: bool
    compatible: bool


@dataclass(frozen=True, slots=True)
class HitRetentionSummary:
    dataset_name: str
    assay_id: str
    role: str
    total_hit_count: int
    retained_hit_count: int
    truncated: bool


@dataclass(frozen=True, slots=True)
class PlausibleAmplicon:
    dataset_name: str
    assay_id: str
    sequence_id: str
    orientation: str
    source_start: int
    source_end: int
    amplicon_size: int
    forward_hit_rank: int | None
    reverse_hit_rank: int | None
    probe_hit_ranks: tuple[int, ...]
    primer_amplicon_plausible: bool
    detectable_off_target: bool
```

`SpecificityResult` contains status, datasets metadata, retained hits, all plausible amplicons, retention summaries, paths to the two TSVs or `None`, and report path.

- [ ] **Step 4: Implement retention after geometry, not before**

For each dataset, compute **all** hits and all geometry first. Then group all hits by `(dataset_name, assay_id, role)`, sort by the documented global key, retain at most `max_hits_per_oligo_per_dataset`, and assign public `hit_rank` within the retained group.

Because geometry can reference unretained internal hits, public `PlausibleAmplicon.forward_hit_rank`, `.reverse_hit_rank`, and `.probe_hit_ranks` may be `None`/omit individual ranks when the supporting hit is not retained. The amplicon row must still exist and its scientific booleans must remain true. This explicitly prevents truncation-induced false safety.

- [ ] **Step 5: Write failing truncation invariance test**

Build a sequence with an early decoy hit plus a later valid F/R/probe geometry and set `max_hits_per_oligo_per_dataset=1`. Assert:

```python
self.assertTrue(any(a.detectable_off_target for a in result.amplicons))
summary = next(
    item for item in result.retention
    if item.dataset_name == "human" and item.assay_id == "a1" and item.role == "FORWARD"
)
self.assertGreater(summary.total_hit_count, summary.retained_hit_count)
self.assertTrue(summary.truncated)
```

- [ ] **Step 6: Implement TSV/JSON schema and dataset provenance**

`off_target_hits.tsv` columns:

```python
OFF_TARGET_HIT_COLUMNS = (
    "dataset_name", "assay_id", "sequence_id", "role", "orientation", "hit_rank",
    "source_start", "source_end", "mismatch_positions", "mismatch_count",
    "exact_match", "three_prime_mismatch", "compatible",
)
```

`plausible_amplicons.tsv` columns:

```python
PLAUSIBLE_AMPLICON_COLUMNS = (
    "dataset_name", "assay_id", "sequence_id", "orientation",
    "source_start", "source_end", "amplicon_size",
    "forward_hit_rank", "reverse_hit_rank", "probe_hit_ranks",
    "primer_amplicon_plausible", "detectable_off_target",
)
```

The JSON top level must contain `schema_version`, `status`, `enabled`, `configuration`, `datasets`, `counts`, `retention`, `artifacts`. Dataset objects contain only provenance: name, source_type, source_path, sha256, sequence_ids, and for frozen data a normalized manifest subset containing `source` and `resolved_entries`; do not embed raw records.

Counts include at least `datasets`, `sequences`, `assays`, `total_compatible_hits`, `retained_hits`, `plausible_amplicons`, `detectable_off_targets`, and `assays_with_detectable_off_target`.

- [ ] **Step 7: Write and implement stale cleanup / atomic JSON-last publication**

Tests must prove:
- disabled run removes only exact old `off_target_hits.tsv` and `plausible_amplicons.tsv`, preserving unrelated siblings;
- enabled empty data writes header-only TSVs;
- old report is invalidated before any data replacement;
- report is replaced last;
- a simulated data replacement failure leaves no new report claiming success;
- LF newlines are stable.

Use the established unique-temporary sibling pattern with `uuid.uuid4().hex` and `finally: temporary.unlink(missing_ok=True)`.

- [ ] **Step 8: Run artifact/service tests and commit**

Run:

```bash
python -m unittest tests.test_specificity_artifacts -q
```

Commit:

```bash
git add qpcr_pipeline/specificity.py tests/test_specificity_artifacts.py
git commit -m "feat: publish off-target specificity results"
```

---

### Task 6: Pipeline integration and semantic regression protection

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Create: `tests/test_pipeline_specificity.py`

**Interfaces:**
- Consumes: `config.off_targets`, `config.specificity`, `PrimerDesignResult`, and `evaluate_specificity()`.
- Produces: pipeline specificity stage after inclusivity and `qc_report.json["specificity"]`.

- [ ] **Step 1: Write failing disabled pipeline summary test**

Patch or fake upstream stages using the same patterns already used by `tests/test_minimal_run.py`. Assert a normal configuration with specificity omitted produces:

```python
{
    "status": "SKIPPED",
    "dataset_count": 0,
    "sequence_count": 0,
    "assay_count": 0,
    "retained_hit_count": 0,
    "plausible_amplicon_count": 0,
    "detectable_off_target_count": 0,
}
```

and only `specificity/specificity_report.json` for the stage.

- [ ] **Step 2: Verify RED**

Run:

```bash
python -m unittest tests.test_pipeline_specificity -v
```

Expected: missing specificity QC key or missing stage call.

- [ ] **Step 3: Integrate service after inclusivity**

Add import and immediately after `evaluate_inclusivity(...)` call:

```python
specificity = evaluate_specificity(
    primer_design,
    config.off_targets,
    config.specificity,
    output_dir,
)
```

Add QC summary from the result, counting upstream assays only when status COMPLETE and counting retained hits/amplicons/detectable rows from the public result.

- [ ] **Step 4: Write enabled FASTA pipeline test**

Create a tiny off-target FASTA in a temporary directory and an upstream fake Primer3 assay whose F/R/probe sequences are present in one off-target record. Enable specificity while allowing inclusivity to remain disabled. Assert:

```python
self.assertEqual(qc["specificity"]["status"], "COMPLETE")
self.assertEqual(qc["specificity"]["dataset_count"], 1)
self.assertGreaterEqual(qc["specificity"]["plausible_amplicon_count"], 1)
self.assertGreaterEqual(qc["specificity"]["detectable_off_target_count"], 1)
```

This proves specificity depends on primer design, not on enabled inclusivity.

- [ ] **Step 5: Write no-network frozen-dataset test**

Materialize a frozen NCBI fixture with existing test helpers, then patch `qpcr_pipeline.ncbi.acquire_ncbi_dataset` and `BioEntrezClient.from_environment` to raise `AssertionError` if touched. Run specificity with the frozen off-target and assert COMPLETE. The test passes only if no network acquisition path is invoked.

- [ ] **Step 6: Add semantic IUPAC equivalence regression**

For small cases where inclusivity's public `evaluate_inclusivity()` can observe the same oligo/target relationship, assert specificity and inclusivity agree on exact/compatible outcomes for:
- canonical exact match;
- degenerate oligo `R` against target `A`;
- conservative failure `R` against ambiguous target `N`;
- one permitted internal mismatch;
- one rejected final-3-prime primer mismatch.

Do not import `_compatibility`, `_enumerate_hits`, or any other private inclusivity helper.

- [ ] **Step 7: Run integration-focused tests and commit**

Run:

```bash
python -m unittest tests.test_pipeline_specificity tests.test_minimal_run -q
```

Commit:

```bash
git add qpcr_pipeline/pipeline.py tests/test_pipeline_specificity.py
git commit -m "feat: integrate off-target specificity stage"
```

---

### Task 7: Documentation and complete acceptance regression

**Files:**
- Modify: `README.md`
- Modify: tests from Tasks 1-6 only if a documented acceptance gap requires a regression assertion.

**Interfaces:**
- Consumes: the final public configuration/artifact contracts.
- Produces: user-facing issue-#9 documentation and evidence for all updated issue criteria.

- [ ] **Step 1: Document configuration and scientific interpretation**

Add a `## Especificidade contra off-targets` section after inclusivity containing the exact effective defaults:

```yaml
off_targets:
  - name: human
    fasta: data/human_subset.fasta
  - name: near_neighbors
    frozen_dataset: runs/ncbi_near_neighbors

specificity:
  enabled: false
  max_hits_per_oligo_per_dataset: 20
  max_primer_mismatches: 2
  max_probe_mismatches: 1
  reject_primer_3_prime_mismatch: true
  primer_3_prime_bases: 5
  max_amplicon_size: 1000
```

Explain:
- FASTA or already-frozen NCBI only;
- no network and no BLAST in the MVP;
- conservative IUPAC matching;
- distinction among isolated hit, plausible F/R amplicon, and probe-detectable off-target;
- truncation affects only `off_target_hits.tsv`, never geometry/risk;
- the two TSVs plus JSON;
- pure-Python exhaustive search is intended for small/moderate curated off-target sets, not whole large genomic banks;
- a future indexed backend must preserve the same scientific output contract.

- [ ] **Step 2: Run focused issue-#9 tests**

Run:

```bash
python -m unittest \
  tests.test_specificity_config \
  tests.test_specificity_matching \
  tests.test_specificity_geometry \
  tests.test_specificity_artifacts \
  tests.test_pipeline_specificity -q
```

Expected: zero failures/errors.

- [ ] **Step 3: Run full offline regression**

Run:

```bash
python -m unittest discover -s tests -q
python -m pytest -q
git diff --check
```

Expected: zero failures/errors; `git diff --check` prints nothing.

- [ ] **Step 4: Commit documentation**

```bash
git add README.md
git commit -m "docs: document off-target specificity"
```

---

## Final review, integration, and CI acceptance

- [ ] Invoke `superpowers:requesting-code-review` against the issue #9 spec, this plan, and the complete branch diff. If no subagent is available, perform an explicit diff review and state that limitation.
- [ ] Resolve every confirmed Critical or Important review finding before integration. Behavior fixes require a focused failing regression test first and `superpowers:systematic-debugging` when the symptom is not immediately understood.
- [ ] Invoke `superpowers:verification-before-completion` and obtain fresh evidence from the available execution environment. Because CircleCI intentionally does not run this feature branch, do not claim the merged state is green before CI executes on `develop`.
- [ ] Open a PR from `feature/issue-9-specificity` to `develop`, review the exact changed-file set, and merge only after the branch implementation satisfies the plan review gate.
- [ ] After merge, query the merged `develop` commit's `ci/circleci: test` status. If it is `pending`, wait for a later status check rather than claiming success. If `failure`, inspect the CI result, invoke `superpowers:systematic-debugging`, fix on a new branch/PR, and rerun.
- [ ] Close GitHub issue #9 as `completed` only after CircleCI reports `success` on the merged `develop` commit and the README/acceptance criteria are present.
- [ ] Do not require the user's manual day-13 validation to close the technical issue; any later user-observed defect becomes a tracked bug/regression.
