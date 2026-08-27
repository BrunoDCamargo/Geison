# Issue #10 Assay Ranking and Classification Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an offline, auditable final assay classification and ranking stage that classifies before scoring, preserves every Primer3 assay, publishes structured reasons plus decomposable scores, and replaces the root report only when ranking is enabled.

**Architecture:** Add a `RankingConfig` contract to configuration, a pure ranking core in `qpcr_pipeline/ranking.py`, and a separate static renderer in `qpcr_pipeline/assay_report_html.py`. The ranking stage consumes `PrimerDesignResult`, `InclusivityResult`, and `SpecificityResult` directly, validates complete evidence matrices, classifies from structured reasons, calculates absolute score components, sorts class-first with deterministic ties, then atomically publishes TSV, JSON, and final HTML.

**Tech Stack:** Python 3.12, stdlib dataclasses/json/html/pathlib/math/uuid, PyYAML already used by the project, unittest-style tests collected by pytest, existing Geison result dataclasses.

**Spec:** `docs/superpowers/specs/2026-08-27-issue-10-ranking-design.md`

## Global Constraints

- Ranking is disabled by default.
- Default inclusivity thresholds are `min_inclusivity_for_pass: 1.0` and `min_inclusivity_before_high_risk: 0.90`.
- Default score weights are inclusivity `0.35`, specificity `0.25`, conservation `0.20`, Primer3 quality `0.10`, robustness `0.10`; weights must sum to `1.0` within absolute tolerance `1e-9`.
- Classification is computed before quantitative ranking; score never overrides class.
- Class order is `IN SILICO PASS`, then `REVIEW`, then `HIGH_RISK`.
- A detectable off-target is always `HIGH_RISK`; this cannot be weakened through configuration.
- Only original Primer3 assays are ranked. Accepted IUPAC proposals are advisory robustness evidence and never become independent ranked assays.
- Missing evidence is never replaced by zero. Any unavailable score component yields `score_status=INCOMPLETE` and `final_score=None`, even when that component has weight zero.
- Ranking is offline: no network access, no BLAST, no remote report resource, and no new third-party dependency.
- When ranking is disabled, it must leave root `report.html` untouched so the conservation report remains backward compatible.
- When ranking is enabled, it owns root `report.html`, removes stale ranking artifacts before upstream evidence validation, and publishes the final assay report only after all output content is computed successfully.
- All text artifacts are written atomically through a temporary sibling file and `Path.replace()`.
- Normal validation runs through `python -m pytest -q`; do not change CircleCI configuration for this issue.

---

## File Structure

**Create**

- `qpcr_pipeline/ranking.py` — ranking dataclasses, evidence integrity checks, reason generation, classification, score components, deterministic ordering, stage artifact publication.
- `qpcr_pipeline/assay_report_html.py` — self-contained static final assay HTML renderer with safe escaping and no remote resources.
- `tests/ranking_fixtures.py` — focused constructors for valid Primer3, inclusivity, proposal, and specificity result objects used across ranking tests.
- `tests/test_ranking_config.py` — configuration parsing and validation.
- `tests/test_ranking_classification.py` — evidence integrity, reason codes, class precedence, inclusivity and specificity rules.
- `tests/test_ranking_scoring.py` — five component formulas and incomplete-score semantics.
- `tests/test_ranking_ordering.py` — class-first and deterministic tie ordering.
- `tests/test_ranking_artifacts.py` — skipped/complete/failure cleanup, deterministic TSV/JSON, zero-assay behavior.
- `tests/test_assay_report_html.py` — final report content, escaping, empty state, and offline guarantees.
- `tests/test_pipeline_ranking.py` — orchestration and `qc_report.json` summary.

**Modify**

- `qpcr_pipeline/config.py` — `RankingWeights`, `RankingConfig`, YAML parsing, validation, and `PipelineConfig.ranking`.
- `qpcr_pipeline/pipeline.py` — call ranking after specificity and publish the ranking QC summary.
- `README.md` — document ranking configuration, classes, score semantics, artifacts, limitations, and root report ownership.

**Do not modify**

- `qpcr_pipeline/report_html.py` — remains the conservation-specific renderer.
- `.circleci/config.yml` — issue #10 adds no CI dependency or job.

---

### Task 1: Add the ranking configuration contract

**Files:**
- Modify: `qpcr_pipeline/config.py`
- Create: `tests/test_ranking_config.py`

**Interfaces:**
- Produces: `RankingWeights`, `RankingConfig`, `validate_ranking_config(config: RankingConfig) -> None`.
- Extends: `PipelineConfig.ranking: RankingConfig`.
- YAML contract:

```yaml
ranking:
  enabled: false
  min_inclusivity_for_pass: 1.0
  min_inclusivity_before_high_risk: 0.90
  weights:
    inclusivity: 0.35
    specificity: 0.25
    conservation: 0.20
    primer3_quality: 0.10
    robustness: 0.10
```

- [ ] **Step 1: Write failing configuration tests**

Create `tests/test_ranking_config.py` with concrete default, parse, invalid-value, and dependency tests:

```python
import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.config import (
    AlignmentConfig,
    ConservationConfig,
    PipelineConfig,
    PrimerDesignConfig,
    RankingConfig,
    RankingWeights,
    load_config,
    validate_ranking_config,
)


FIXTURE_FASTA = Path("tests/fixtures/target_small.fasta")


class RankingConfigTests(unittest.TestCase):
    def _load_yaml(self, text: str) -> PipelineConfig:
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "config.yaml"
            path.write_text(text, encoding="utf-8")
            return load_config(path)

    def _enabled_base(self) -> str:
        return (
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
            "alignment:\n  enabled: true\n"
            "conservation:\n  enabled: true\n"
            "primer_design:\n  enabled: true\n"
        )

    def test_ranking_defaults_are_disabled_and_conservative(self):
        config = self._load_yaml(
            "target:\n  name: target\n"
            f"input:\n  fasta: {FIXTURE_FASTA.as_posix()}\n"
        )
        self.assertEqual(config.ranking, RankingConfig())
        self.assertEqual(config.ranking.weights, RankingWeights())
        self.assertEqual(config.ranking.min_inclusivity_for_pass, 1.0)
        self.assertEqual(config.ranking.min_inclusivity_before_high_risk, 0.90)

    def test_loads_ranking_thresholds_and_weights(self):
        config = self._load_yaml(
            self._enabled_base()
            + "ranking:\n"
            "  enabled: true\n"
            "  min_inclusivity_for_pass: 0.98\n"
            "  min_inclusivity_before_high_risk: 0.85\n"
            "  weights:\n"
            "    inclusivity: 0.40\n"
            "    specificity: 0.30\n"
            "    conservation: 0.15\n"
            "    primer3_quality: 0.10\n"
            "    robustness: 0.05\n"
        )
        self.assertEqual(
            config.ranking,
            RankingConfig(
                enabled=True,
                min_inclusivity_for_pass=0.98,
                min_inclusivity_before_high_risk=0.85,
                weights=RankingWeights(0.40, 0.30, 0.15, 0.10, 0.05),
            ),
        )

    def test_rejects_invalid_ranking_values(self):
        invalid = (
            ("ranking:\n  enabled: yes\n", "enabled"),
            ("ranking:\n  min_inclusivity_for_pass: 1.1\n", "min_inclusivity_for_pass"),
            ("ranking:\n  min_inclusivity_before_high_risk: -0.1\n", "min_inclusivity_before_high_risk"),
            (
                "ranking:\n  min_inclusivity_for_pass: 0.8\n"
                "  min_inclusivity_before_high_risk: 0.9\n",
                "cannot exceed",
            ),
            (
                "ranking:\n  weights:\n"
                "    inclusivity: 0.5\n"
                "    specificity: 0.5\n"
                "    conservation: 0.5\n"
                "    primer3_quality: 0.0\n"
                "    robustness: 0.0\n",
                "sum to 1.0",
            ),
            ("ranking:\n  weights:\n    inclusivity: -0.1\n", "non-negative"),
            ("ranking:\n  surprise: 1\n", "unrecognized"),
            ("ranking:\n  weights:\n    surprise: 1\n", "unrecognized"),
        )
        for suffix, message in invalid:
            with self.subTest(suffix=suffix):
                with self.assertRaisesRegex(ValueError, message):
                    self._load_yaml(self._enabled_base() + suffix)

    def test_enabled_ranking_requires_enabled_primer_design(self):
        config = PipelineConfig(
            target_name="target",
            input_fasta=FIXTURE_FASTA,
            alignment=AlignmentConfig(enabled=True),
            conservation=ConservationConfig(enabled=True),
            primer_design=PrimerDesignConfig(enabled=False),
            ranking=RankingConfig(enabled=True),
        )
        with self.assertRaisesRegex(ValueError, "ranking.*requires enabled primer design"):
            config.selected_input

    def test_direct_validator_rejects_wrong_type(self):
        with self.assertRaisesRegex(ValueError, "Ranking configuration"):
            validate_ranking_config(object())


if __name__ == "__main__":
    unittest.main()
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```bash
python -m pytest tests/test_ranking_config.py -q
```

Expected: collection/import fails because `RankingConfig`, `RankingWeights`, and `validate_ranking_config` do not exist.

- [ ] **Step 3: Add immutable ranking configuration models and validation**

Add beside the existing stage configs in `qpcr_pipeline/config.py`:

```python
@dataclass(frozen=True, slots=True)
class RankingWeights:
    inclusivity: float = 0.35
    specificity: float = 0.25
    conservation: float = 0.20
    primer3_quality: float = 0.10
    robustness: float = 0.10


@dataclass(frozen=True, slots=True)
class RankingConfig:
    enabled: bool = False
    min_inclusivity_for_pass: float = 1.0
    min_inclusivity_before_high_risk: float = 0.90
    weights: RankingWeights = field(default_factory=RankingWeights)
```

Extend `PipelineConfig`:

```python
ranking: RankingConfig = field(default_factory=RankingConfig)
```

Implement exact validation:

```python
def validate_ranking_config(config: RankingConfig) -> None:
    if not isinstance(config, RankingConfig):
        raise ValueError("Ranking configuration must be a RankingConfig.")
    if not isinstance(config.enabled, bool):
        raise ValueError("Ranking enabled must be a boolean.")
    for name in (
        "min_inclusivity_for_pass",
        "min_inclusivity_before_high_risk",
    ):
        value = getattr(config, name)
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(value)
            or not 0.0 <= value <= 1.0
        ):
            raise ValueError(f"Ranking {name} must be a finite number between 0 and 1.")
    if config.min_inclusivity_before_high_risk > config.min_inclusivity_for_pass:
        raise ValueError(
            "Ranking min_inclusivity_before_high_risk cannot exceed "
            "min_inclusivity_for_pass."
        )
    if not isinstance(config.weights, RankingWeights):
        raise ValueError("Ranking weights must be a RankingWeights.")
    values = tuple(getattr(config.weights, field_name) for field_name in (
        "inclusivity", "specificity", "conservation", "primer3_quality", "robustness"
    ))
    if any(
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
        or value < 0.0
        for value in values
    ):
        raise ValueError("Ranking weights must be finite non-negative numbers.")
    if not math.isclose(sum(values), 1.0, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("Ranking weights must sum to 1.0.")
```

- [ ] **Step 4: Parse ranking YAML and enforce the pipeline dependency**

Add `_parse_ranking_config(raw: Any) -> RankingConfig` with allowed top-level fields `enabled`, `min_inclusivity_for_pass`, `min_inclusivity_before_high_risk`, `weights`; allow only the five named weight keys; use `RankingConfig()` and `RankingWeights()` for defaults; call `validate_ranking_config` before returning.

In `load_config()` parse `raw.get("ranking", {})` and pass it into `PipelineConfig`. In `validate_pipeline_config()` type-check `config.ranking`, call `validate_ranking_config(config.ranking)`, and add:

```python
if config.ranking.enabled and not config.primer_design.enabled:
    raise ValueError("Enabled ranking requires enabled primer design.")
```

- [ ] **Step 5: Run the focused config suite and normal config regressions**

Run:

```bash
python -m pytest tests/test_ranking_config.py tests/test_specificity_config.py -q
```

Expected: PASS.

- [ ] **Step 6: Commit the configuration slice**

```bash
git add qpcr_pipeline/config.py tests/test_ranking_config.py
git commit -m "feat: add ranking configuration"
```

---

### Task 2: Build evidence integrity, structured reasons, and classification

**Files:**
- Create: `qpcr_pipeline/ranking.py`
- Create: `tests/ranking_fixtures.py`
- Create: `tests/test_ranking_classification.py`

**Interfaces:**
- Produces public types:

```python
ReasonSeverity = Literal["HIGH_RISK", "REVIEW", "ADVISORY"]
AssayClassification = Literal["IN SILICO PASS", "REVIEW", "HIGH_RISK"]
ScoreStatus = Literal["COMPLETE", "INCOMPLETE"]

@dataclass(frozen=True, slots=True)
class RankingReason:
    code: str
    severity: ReasonSeverity
    source: str
    message: str
    evidence: tuple[tuple[str, object], ...]

@dataclass(frozen=True, slots=True)
class ScoreComponents:
    inclusivity: float | None
    specificity: float | None
    conservation: float | None
    primer3_quality: float | None
    robustness: float | None

@dataclass(frozen=True, slots=True)
class RankedAssay:
    rank: int
    assay_id: str
    region_id: str
    primer3_index: int
    classification: AssayClassification
    final_score: float | None
    score_status: ScoreStatus
    components: ScoreComponents
    reasons: tuple[RankingReason, ...]
    original_compatible_count: int | None
    evaluation_sequence_count: int | None
    compatible_off_target_hit_count: int | None
    plausible_off_target_count: int | None
    detectable_off_target_count: int | None
    pair_penalty: float | None
```

- Produces pure core: `rank_assays(primer_design: PrimerDesignResult, inclusivity: InclusivityResult, specificity: SpecificityResult, config: RankingConfig) -> tuple[RankedAssay, ...]`.
- Task 2 initially fills score fields using Task 3 helpers added next; classification tests should assert class/reasons rather than exact final scores.

- [ ] **Step 1: Add shared valid result constructors**

Create `tests/ranking_fixtures.py` with deterministic builders. Use one candidate region with perfect conservation by default and full evidence matrices:

```python
from pathlib import Path

from qpcr_pipeline.inclusivity import (
    AssayInclusivity,
    DegeneracyProposal,
    InclusivityResult,
)
from qpcr_pipeline.primer_design import (
    AssayCandidate,
    CandidateRegion,
    DesignedOligo,
    PrimerDesignResult,
)
from qpcr_pipeline.specificity import (
    HitRetentionSummary,
    PlausibleAmplicon,
    SpecificityResult,
)


def make_region(region_id: str = "r1") -> CandidateRegion:
    return CandidateRegion(
        region_id=region_id,
        rank=1,
        reference_start=1,
        reference_end=100,
        peak_start=1,
        peak_end=100,
        position_count=100,
        usable_length=100,
        usable_fraction=1.0,
        mean_conservation=1.0,
        minimum_conservation=1.0,
        mean_coverage=1.0,
        mean_gap_frequency=0.0,
        mean_entropy_bits=0.0,
    )


def make_oligo(sequence: str, start: int) -> DesignedOligo:
    return DesignedOligo(
        sequence=sequence,
        reference_start=start,
        reference_end=start + len(sequence) - 1,
        length=len(sequence),
        tm=60.0,
        gc_percent=50.0,
        penalty=0.5,
        metrics=(),
    )


def make_assay(
    assay_id: str = "a1",
    *,
    region_id: str = "r1",
    primer3_index: int = 0,
    pair_penalty: float | None = 1.0,
) -> AssayCandidate:
    return AssayCandidate(
        assay_id=assay_id,
        region_id=region_id,
        primer3_index=primer3_index,
        forward_primer=make_oligo("ACGT", 1),
        probe=make_oligo("TTAA", 9),
        reverse_primer=make_oligo("AGTC", 17),
        product_size=20,
        pair_penalty=pair_penalty,
        metrics=(),
    )


def make_primer_result(
    assays: tuple[AssayCandidate, ...] | None = None,
    candidates: tuple[CandidateRegion, ...] | None = None,
    *,
    status: str = "COMPLETE",
) -> PrimerDesignResult:
    return PrimerDesignResult(
        status=status,
        reference_id="ref" if status == "COMPLETE" else None,
        candidates=candidates if candidates is not None else (make_region(),),
        assays=assays if assays is not None else (make_assay(),),
        candidate_regions_path=None,
        assays_path=None,
        primer3_input_path=None,
        primer3_output_path=None,
        report_path=Path("primer-report.json"),
    )


def make_inclusivity_result(
    primer: PrimerDesignResult,
    compatibility: dict[str, tuple[bool, ...]] | None = None,
    *,
    sequence_ids: tuple[str, ...] = ("s1",),
    proposals: tuple[DegeneracyProposal, ...] = (),
    status: str = "COMPLETE",
) -> InclusivityResult:
    if status == "SKIPPED":
        sequence_ids = ()
        rows = ()
    else:
        configured = compatibility or {
            assay.assay_id: tuple(True for _ in sequence_ids) for assay in primer.assays
        }
        rows = tuple(
            AssayInclusivity(
                assay_id=assay.assay_id,
                sequence_id=sequence_id,
                orientation="FORWARD",
                geometry_found=True,
                source_amplicon_start=1,
                source_amplicon_end=20,
                amplicon_size=20,
                forward_match=None,
                probe_match=None,
                reverse_match=None,
                original_compatible=configured[assay.assay_id][index],
                proposed_forward=None,
                proposed_probe=None,
                proposed_reverse=None,
                proposed_compatible=configured[assay.assay_id][index],
            )
            for assay in primer.assays
            for index, sequence_id in enumerate(sequence_ids)
        )
    return InclusivityResult(
        status=status,
        evaluation_sequence_ids=sequence_ids,
        oligo_matches=(),
        assay_results=rows,
        variations=(),
        proposals=proposals,
        oligo_matches_path=None,
        assay_inclusivity_path=None,
        oligo_variations_path=None,
        degeneracy_proposals_path=None,
        report_path=Path("inclusivity-report.json"),
    )


def make_specificity_result(
    primer: PrimerDesignResult,
    *,
    dataset_names: tuple[str, ...] = ("off",),
    hit_totals: dict[tuple[str, str, str], int] | None = None,
    amplicons: tuple[PlausibleAmplicon, ...] = (),
    status: str = "COMPLETE",
) -> SpecificityResult:
    if status == "SKIPPED":
        return SpecificityResult(
            status="SKIPPED",
            dataset_names=(),
            sequence_count=0,
            assay_count=0,
            hits=(),
            amplicons=(),
            retention=(),
            off_target_hits_path=None,
            plausible_amplicons_path=None,
            report_path=Path("specificity-report.json"),
        )
    totals = hit_totals or {}
    retention = tuple(
        HitRetentionSummary(
            dataset_name=dataset,
            assay_id=assay.assay_id,
            role=role,
            total_hit_count=totals.get((dataset, assay.assay_id, role), 0),
            retained_hit_count=min(totals.get((dataset, assay.assay_id, role), 0), 20),
            truncated=totals.get((dataset, assay.assay_id, role), 0) > 20,
        )
        for dataset in dataset_names
        for assay in primer.assays
        for role in ("FORWARD", "PROBE", "REVERSE")
    )
    return SpecificityResult(
        status="COMPLETE",
        dataset_names=dataset_names,
        sequence_count=1,
        assay_count=len(primer.assays),
        hits=(),
        amplicons=amplicons,
        retention=retention,
        off_target_hits_path=None,
        plausible_amplicons_path=None,
        report_path=Path("specificity-report.json"),
    )
```

- [ ] **Step 2: Write RED tests for inclusivity classes, specificity classes, missing evidence, and proposal advisories**

Create `tests/test_ranking_classification.py`. Include these explicit cases:

```python
import unittest

from qpcr_pipeline.config import RankingConfig
from qpcr_pipeline.inclusivity import DegeneracyProposal
from qpcr_pipeline.ranking import rank_assays
from qpcr_pipeline.specificity import PlausibleAmplicon
from tests.ranking_fixtures import (
    make_inclusivity_result,
    make_primer_result,
    make_specificity_result,
)


def amplicon(*, detectable: bool) -> PlausibleAmplicon:
    return PlausibleAmplicon(
        dataset_name="off",
        assay_id="a1",
        sequence_id="x",
        orientation="FORWARD",
        source_start=1,
        source_end=20,
        amplicon_size=20,
        forward_source_start=1,
        forward_source_end=4,
        reverse_source_start=17,
        reverse_source_end=20,
        probe_source_sites=((9, 12),) if detectable else (),
        forward_hit_rank=1,
        reverse_hit_rank=1,
        probe_hit_ranks=(1,) if detectable else (),
        primer_amplicon_plausible=True,
        detectable_off_target=detectable,
    )


class RankingClassificationTests(unittest.TestCase):
    def test_default_inclusivity_boundaries_are_pass_review_high_risk(self):
        primer = make_primer_result()
        cases = (
            ((True,) * 10, "IN SILICO PASS"),
            ((True,) * 9 + (False,), "REVIEW"),
            ((True,) * 8 + (False, False), "HIGH_RISK"),
        )
        sequence_ids = tuple(f"s{i}" for i in range(10))
        for compatibility, expected in cases:
            with self.subTest(expected=expected):
                result = rank_assays(
                    primer,
                    make_inclusivity_result(
                        primer,
                        {"a1": compatibility},
                        sequence_ids=sequence_ids,
                    ),
                    make_specificity_result(primer),
                    RankingConfig(enabled=True),
                )
                self.assertEqual(result[0].classification, expected)

    def test_detectable_off_target_forces_high_risk(self):
        primer = make_primer_result()
        result = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(primer, amplicons=(amplicon(detectable=True),)),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "HIGH_RISK")
        self.assertIn("DETECTABLE_OFF_TARGET", [reason.code for reason in result[0].reasons])

    def test_plausible_amplicon_without_probe_forces_review(self):
        primer = make_primer_result()
        result = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(primer, amplicons=(amplicon(detectable=False),)),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "REVIEW")
        self.assertIn(
            "PLAUSIBLE_OFF_TARGET_AMPLICON",
            [reason.code for reason in result[0].reasons],
        )

    def test_isolated_hits_are_advisory_only(self):
        primer = make_primer_result()
        result = rank_assays(
            primer,
            make_inclusivity_result(primer),
            make_specificity_result(
                primer,
                hit_totals={("off", "a1", "FORWARD"): 2},
            ),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "IN SILICO PASS")
        self.assertIn("ISOLATED_OFF_TARGET_HITS", [reason.code for reason in result[0].reasons])

    def test_skipped_sources_force_review_but_not_high_risk(self):
        primer = make_primer_result()
        result = rank_assays(
            primer,
            make_inclusivity_result(primer, status="SKIPPED"),
            make_specificity_result(primer, status="SKIPPED"),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "REVIEW")
        self.assertEqual(
            {reason.code for reason in result[0].reasons},
            {
                "EVIDENCE_INCOMPLETE",
                "INCLUSIVITY_EVIDENCE_MISSING",
                "SPECIFICITY_EVIDENCE_MISSING",
            },
        )

    def test_detectable_risk_wins_even_when_inclusivity_is_missing(self):
        primer = make_primer_result()
        result = rank_assays(
            primer,
            make_inclusivity_result(primer, status="SKIPPED"),
            make_specificity_result(primer, amplicons=(amplicon(detectable=True),)),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "HIGH_RISK")

    def test_iupac_proposal_is_advisory_and_does_not_replace_original_assay(self):
        primer = make_primer_result()
        proposal = DegeneracyProposal(
            assay_id="a1",
            role="FORWARD",
            original_sequence="ACGT",
            proposed_sequence="ARGT",
            status="ACCEPTED",
            reason="coverage improved",
            original_degeneracy=1,
            proposed_degeneracy=2,
            changed_positions=(2,),
            binding_site_count=1,
            original_exact_count=0,
            original_exact_fraction=0.0,
            proposed_exact_count=1,
            proposed_exact_fraction=1.0,
        )
        result = rank_assays(
            primer,
            make_inclusivity_result(primer, proposals=(proposal,)),
            make_specificity_result(primer),
            RankingConfig(enabled=True),
        )
        self.assertEqual(result[0].classification, "IN SILICO PASS")
        self.assertIn("IUPAC_PROPOSAL_ACCEPTED", [reason.code for reason in result[0].reasons])
```

- [ ] **Step 3: Run classification tests and verify RED**

```bash
python -m pytest tests/test_ranking_classification.py -q
```

Expected: import fails because `qpcr_pipeline.ranking` does not exist.

- [ ] **Step 4: Implement evidence validation before any classification**

In `qpcr_pipeline/ranking.py`, create `RankingError(RuntimeError)` and validators that run before reasons are generated:

```python
_ROLE_ORDER = ("FORWARD", "PROBE", "REVERSE")
_CLASS_ORDER = {"IN SILICO PASS": 0, "REVIEW": 1, "HIGH_RISK": 2}
_SEVERITY_ORDER = {"HIGH_RISK": 0, "REVIEW": 1, "ADVISORY": 2}


def _validate_primer_evidence(primer_design: PrimerDesignResult) -> dict[str, CandidateRegion]:
    if primer_design.status != "COMPLETE":
        raise RankingError("Enabled ranking requires COMPLETE primer design evidence.")
    assay_ids = tuple(assay.assay_id for assay in primer_design.assays)
    if len(set(assay_ids)) != len(assay_ids):
        raise RankingError("Primer design assay_id values must be unique for ranking.")
    regions = {candidate.region_id: candidate for candidate in primer_design.candidates}
    if len(regions) != len(primer_design.candidates):
        raise RankingError("Primer design region_id values must be unique for ranking.")
    for assay in primer_design.assays:
        if assay.region_id not in regions:
            raise RankingError(
                f"Assay {assay.assay_id} references unknown region {assay.region_id}."
            )
    return regions
```

For `InclusivityResult.status == "COMPLETE"`, require unique `evaluation_sequence_ids`, exactly one `AssayInclusivity` for every Cartesian `(assay_id, sequence_id)` pair, no unknown pair, and at most one proposal for every `(assay_id, role)`. For `SpecificityResult.status == "COMPLETE"`, require `assay_count == len(primer_design.assays)`, unique `dataset_names`, and exactly one `HitRetentionSummary` for every `(dataset_name, assay_id, FORWARD|PROBE|REVERSE)` tuple. Reject unknown assay IDs in retention and amplicons.

- [ ] **Step 5: Generate deterministic reasons before class**

Implement reason constructors that use evidence tuples sorted by key. Use source `ranking` for the single aggregate `EVIDENCE_INCOMPLETE`, so both missing stages still create only one aggregate reason.

Classification reason rules must be exactly:

```python
if inclusivity.status == "SKIPPED":
    add("INCLUSIVITY_EVIDENCE_MISSING", "REVIEW", "inclusivity", ...)

if specificity.status == "SKIPPED":
    add("SPECIFICITY_EVIDENCE_MISSING", "REVIEW", "specificity", ...)

if inclusivity_fraction is not None:
    if inclusivity_fraction < config.min_inclusivity_before_high_risk:
        add("INCLUSIVITY_BELOW_MINIMUM", "HIGH_RISK", "inclusivity", ...)
    elif inclusivity_fraction < config.min_inclusivity_for_pass:
        add("INCLUSIVITY_BELOW_PASS", "REVIEW", "inclusivity", ...)

if detectable_count:
    add("DETECTABLE_OFF_TARGET", "HIGH_RISK", "specificity", ...)
elif plausible_count:
    add("PLAUSIBLE_OFF_TARGET_AMPLICON", "REVIEW", "specificity", ...)
elif compatible_hit_count:
    add("ISOLATED_OFF_TARGET_HITS", "ADVISORY", "specificity", ...)
```

Add `IUPAC_PROPOSAL_ACCEPTED` and `IUPAC_PROPOSAL_REJECTED` advisory reasons from proposals; `UNCHANGED` adds no reason. Deduplicate by `(code, source)` and sort reasons with:

```python
key=lambda reason: (_SEVERITY_ORDER[reason.severity], reason.source, reason.code)
```

- [ ] **Step 6: Determine class only from reasons**

Implement before score computation:

```python
def _classification(reasons: tuple[RankingReason, ...]) -> AssayClassification:
    if any(reason.severity == "HIGH_RISK" for reason in reasons):
        return "HIGH_RISK"
    if any(reason.severity == "REVIEW" for reason in reasons):
        return "REVIEW"
    return "IN SILICO PASS"
```

Before calling `_classification`, detect raw score prerequisites that are already known to be unavailable: skipped inclusivity, empty Evaluation Set, skipped specificity, or `pair_penalty is None`. Add one `REVIEW / EVIDENCE_INCOMPLETE` from source `ranking` with evidence listing the unavailable component names. This preserves the required reason → class → score execution order.

- [ ] **Step 7: Add integrity regression tests**

Extend `tests/test_ranking_classification.py` with concrete malformed evidence cases: duplicate assay ID, missing region, missing/duplicate inclusivity matrix row, unknown proposal assay, mismatched specificity `assay_count`, missing retention row, duplicate retention row, and unknown amplicon assay. Each must assert `RankingError` rather than `REVIEW`.

Example:

```python
def test_complete_specificity_missing_retention_row_is_invalid(self):
    primer = make_primer_result()
    specificity = make_specificity_result(primer)
    malformed = specificity.__class__(
        status=specificity.status,
        dataset_names=specificity.dataset_names,
        sequence_count=specificity.sequence_count,
        assay_count=specificity.assay_count,
        hits=specificity.hits,
        amplicons=specificity.amplicons,
        retention=specificity.retention[:-1],
        off_target_hits_path=specificity.off_target_hits_path,
        plausible_amplicons_path=specificity.plausible_amplicons_path,
        report_path=specificity.report_path,
    )
    with self.assertRaisesRegex(RankingError, "retention"):
        rank_assays(
            primer,
            make_inclusivity_result(primer),
            malformed,
            RankingConfig(enabled=True),
        )
```

- [ ] **Step 8: Run classification and integrity tests**

```bash
python -m pytest tests/test_ranking_classification.py -q
```

Expected: PASS.

- [ ] **Step 9: Commit the classification slice**

```bash
git add qpcr_pipeline/ranking.py tests/ranking_fixtures.py tests/test_ranking_classification.py
git commit -m "feat: classify assays from ranking evidence"
```

---

### Task 3: Add decomposable scoring and deterministic ordering

**Files:**
- Modify: `qpcr_pipeline/ranking.py`
- Create: `tests/test_ranking_scoring.py`
- Create: `tests/test_ranking_ordering.py`

**Interfaces:**
- Completes `rank_assays(...)` by populating `ScoreComponents`, `score_status`, `final_score`, and contiguous final `rank`.
- Internal formulas are absolute and independent of other assays in the run.

- [ ] **Step 1: Write RED tests for all five score components**

Create `tests/test_ranking_scoring.py`. Use ten Evaluation Set sequences with nine compatible, two isolated specificity hits, perfect candidate-region conservation, `pair_penalty=1.0`, and one accepted forward proposal with degeneracy `1 -> 2`. Assert:

```python
self.assertAlmostEqual(assay.components.inclusivity, 0.9)
self.assertAlmostEqual(assay.components.specificity, 0.96)
self.assertAlmostEqual(assay.components.conservation, 1.0)
self.assertAlmostEqual(assay.components.primer3_quality, 0.5)
self.assertAlmostEqual(assay.components.robustness, (0.5 + 1.0 + 1.0) / 3.0)
self.assertAlmostEqual(assay.final_score, 88.83333333333333)
self.assertEqual(assay.score_status, "COMPLETE")
```

Also test specificity component values `1.0`, floor `0.80`, `0.40`, and `0.0`; test the conservation entropy cap with `mean_entropy_bits=2.0`; and test that changing the other assays in the same run does not change a fixed assay's score.

- [ ] **Step 2: Write RED tests for incomplete scores**

Add exact cases:

```python
def test_missing_pair_penalty_makes_score_incomplete_and_review(self):
    primer = make_primer_result(assays=(make_assay(pair_penalty=None),))
    ranked = rank_assays(
        primer,
        make_inclusivity_result(primer),
        make_specificity_result(primer),
        RankingConfig(enabled=True),
    )[0]
    self.assertEqual(ranked.score_status, "INCOMPLETE")
    self.assertIsNone(ranked.final_score)
    self.assertEqual(ranked.classification, "REVIEW")
    self.assertIn("EVIDENCE_INCOMPLETE", [reason.code for reason in ranked.reasons])
```

Create a second case with `RankingWeights(primer3_quality=0.0, ...)` summing to `1.0`; missing pair penalty must still be incomplete.

Use a complete inclusivity result with `evaluation_sequence_ids=()` to verify no division by zero, inclusivity and robustness are `None`, score is incomplete, and classification is REVIEW.

- [ ] **Step 3: Run scoring tests and verify RED**

```bash
python -m pytest tests/test_ranking_scoring.py -q
```

Expected: failures on component and final score assertions until formulas are implemented.

- [ ] **Step 4: Implement component formulas exactly**

Add helpers to `ranking.py`:

```python
def _inclusivity_component(compatible_count: int, total: int) -> float | None:
    return None if total == 0 else compatible_count / total


def _specificity_component(
    *,
    hit_count: int,
    plausible_count: int,
    detectable_count: int,
) -> float:
    if detectable_count:
        return 0.0
    if plausible_count:
        return 0.40
    if hit_count:
        return max(0.80, 1.0 - 0.02 * hit_count)
    return 1.0


def _conservation_component(region: CandidateRegion) -> float:
    return (
        region.mean_conservation
        + region.minimum_conservation
        + region.mean_coverage
        + (1.0 - region.mean_gap_frequency)
        + (1.0 - min(region.mean_entropy_bits / 2.0, 1.0))
    ) / 5.0


def _primer3_quality_component(pair_penalty: float | None) -> float | None:
    return None if pair_penalty is None else 1.0 / (1.0 + pair_penalty)
```

Validate every present fraction as finite and in `[0, 1]`, entropy as finite and non-negative, and pair penalty as finite and non-negative. Invalid present metrics raise `RankingError`; only explicit absence becomes incomplete.

For robustness, use `1.0` for every role without an accepted proposal and `original_degeneracy / proposed_degeneracy` for an accepted proposal. Validate both degeneracies are positive and the proposed value is not lower than the original.

Compute `final_score` only when all five components are present:

```python
final_score = 100.0 * (
    config.weights.inclusivity * components.inclusivity
    + config.weights.specificity * components.specificity
    + config.weights.conservation * components.conservation
    + config.weights.primer3_quality * components.primer3_quality
    + config.weights.robustness * components.robustness
)
```

Keep full precision; do not round in the model.

- [ ] **Step 5: Run scoring tests and verify GREEN**

```bash
python -m pytest tests/test_ranking_scoring.py -q
```

Expected: PASS.

- [ ] **Step 6: Write RED ordering tests**

Create `tests/test_ranking_ordering.py` with at least four assays arranged so their numeric scores contradict their safety class. Assert final IDs are still class-first. Add a same-class pair where one score is incomplete, then explicit tie cases for inclusivity, pair penalty, Primer3 index, and assay ID.

Core assertion example:

```python
self.assertEqual(
    [item.classification for item in ranked],
    ["IN SILICO PASS", "REVIEW", "HIGH_RISK"],
)
self.assertEqual([item.rank for item in ranked], [1, 2, 3])
```

For an intentionally high-scoring HIGH_RISK assay and lower-scoring REVIEW assay:

```python
self.assertLess(
    next(item.rank for item in ranked if item.assay_id == "review"),
    next(item.rank for item in ranked if item.assay_id == "risky"),
)
```

- [ ] **Step 7: Implement deterministic sort key and contiguous ranks**

Use present/missing sentinel tuples rather than fake numeric defaults:

```python
def _optional_desc(value: float | None) -> tuple[int, float]:
    return (1, 0.0) if value is None else (0, -value)


def _optional_asc(value: float | None) -> tuple[int, float]:
    return (1, 0.0) if value is None else (0, value)


def _sort_key(item: RankedAssay) -> tuple[object, ...]:
    return (
        _CLASS_ORDER[item.classification],
        0 if item.score_status == "COMPLETE" else 1,
        *_optional_desc(item.final_score),
        *_optional_desc(item.components.inclusivity),
        *_optional_asc(item.pair_penalty),
        item.primer3_index,
        item.assay_id,
    )
```

Sort once, then create the final frozen objects with `dataclasses.replace(item, rank=index)` for `index` starting at 1.

- [ ] **Step 8: Run scoring, ordering, and classification suites together**

```bash
python -m pytest \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py -q
```

Expected: PASS.

- [ ] **Step 9: Commit scoring and ordering**

```bash
git add qpcr_pipeline/ranking.py tests/test_ranking_scoring.py tests/test_ranking_ordering.py
git commit -m "feat: score and order final assays"
```

---

### Task 4: Render the final self-contained assay report

**Files:**
- Create: `qpcr_pipeline/assay_report_html.py`
- Create: `tests/test_assay_report_html.py`

**Interfaces:**
- Produces:

```python
def render_assay_report_html(
    *,
    target_name: str,
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    assays: tuple[RankedAssay, ...],
) -> str:
    ...
```

- Does not read artifact files and does not execute JavaScript.

- [ ] **Step 1: Write RED renderer tests**

Create `tests/test_assay_report_html.py` with a PASS assay, an accepted IUPAC proposal, and specificity evidence. Assert that the returned document:

```python
self.assertIn("<!doctype html>", html_text.lower())
self.assertIn("IN SILICO PASS", html_text)
self.assertIn("ACGT", html_text)  # original forward
self.assertIn("ARGT", html_text)  # proposal displayed as context
self.assertIn("inclusivity", html_text.lower())
self.assertIn("specificity", html_text.lower())
self.assertIn("conservation", html_text.lower())
self.assertIn("primer3 quality", html_text.lower())
self.assertIn("robustness", html_text.lower())
self.assertNotIn("https://", html_text.lower())
self.assertNotIn("http://", html_text.lower())
self.assertNotIn("<script", html_text.lower())
```

Pass `target_name='<img src=x onerror="boom">'` and assert the literal unsafe tag is absent while escaped text is present:

```python
self.assertNotIn('<img src=x onerror="boom">', html_text)
self.assertIn("&lt;img", html_text)
```

Add zero-assay input and assert a visible `No assay candidates` empty state.

- [ ] **Step 2: Run the renderer test and verify RED**

```bash
python -m pytest tests/test_assay_report_html.py -q
```

Expected: import fails because `assay_report_html.py` does not exist.

- [ ] **Step 3: Implement static HTML helpers with escaping at every dynamic boundary**

Use stdlib `html.escape` and no JSON-in-script payload. The renderer should construct a summary table plus one `<details>` section per ranked assay. Add deterministic helper maps keyed by `assay_id` and `(assay_id, role)` so output follows `assays` order and roles always follow `FORWARD`, `PROBE`, `REVERSE`.

Use a complete document skeleton:

```python
from html import escape


def _text(value: object) -> str:
    return escape("" if value is None else str(value), quote=True)


def _number(value: float | None, digits: int = 2) -> str:
    return "Not available" if value is None else f"{value:.{digits}f}"


def render_assay_report_html(...):
    summary_rows = "".join(_summary_row(item) for item in assays)
    detail_blocks = "".join(
        _detail_block(item, primer_by_id[item.assay_id], region_by_id[primer_by_id[item.assay_id].region_id], ...)
        for item in assays
    )
    empty_state = "<p>No assay candidates are available.</p>" if not assays else ""
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Geison assay report</title>
<style>
body{{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;max-width:80rem;margin:0 auto;padding:1.5rem;color:#172033}}
table{{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums}}
th,td{{padding:.5rem;border-bottom:1px solid #d8deea;text-align:left;vertical-align:top}}
details{{border:1px solid #d8deea;border-radius:.5rem;padding:.75rem;margin:.75rem 0}}
code{{overflow-wrap:anywhere}}
.grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(14rem,1fr));gap:.75rem}}
.panel{{border:1px solid #d8deea;border-radius:.4rem;padding:.65rem}}
</style>
</head>
<body>
<main>
<h1>Final qPCR assay ranking</h1>
<p><strong>Target:</strong> {_text(target_name)}</p>
{empty_state}
<table><thead><tr><th>Rank</th><th>Assay</th><th>Class</th><th>Score</th><th>Reasons</th></tr></thead>
<tbody>{summary_rows}</tbody></table>
{detail_blocks}
</main>
</body>
</html>
"""
```

Each detail block must include original F/Probe/R sequence, coordinates, Tm, GC, oligo penalty, product size, pair penalty, candidate-region metrics, original inclusivity count/fraction, IUPAC proposal status/sequence/degeneracy, specificity counts by dataset, reason codes, and all five score components. Never substitute a proposed IUPAC sequence into the original oligo row.

- [ ] **Step 4: Run renderer tests and verify GREEN**

```bash
python -m pytest tests/test_assay_report_html.py -q
```

Expected: PASS.

- [ ] **Step 5: Commit the renderer**

```bash
git add qpcr_pipeline/assay_report_html.py tests/test_assay_report_html.py
git commit -m "feat: render final assay report"
```

---

### Task 5: Publish ranking TSV/JSON/HTML with correct stale-artifact ownership

**Files:**
- Modify: `qpcr_pipeline/ranking.py`
- Create: `tests/test_ranking_artifacts.py`

**Interfaces:**
- Adds public stage result:

```python
@dataclass(frozen=True, slots=True)
class RankingResult:
    status: Literal["SKIPPED", "COMPLETE"]
    assays: tuple[RankedAssay, ...]
    ranking_tsv_path: Path | None
    ranking_report_path: Path
    html_report_path: Path | None
```

- Adds public stage function:

```python
def evaluate_ranking(
    primer_design: PrimerDesignResult,
    inclusivity: InclusivityResult,
    specificity: SpecificityResult,
    config: RankingConfig,
    output_dir: Path,
    *,
    target_name: str,
) -> RankingResult:
    ...
```

- [ ] **Step 1: Write RED tests for disabled ownership and stale cleanup**

Create `tests/test_ranking_artifacts.py`. First prove the approved ownership correction:

```python
def test_disabled_ranking_preserves_existing_root_report(self):
    with tempfile.TemporaryDirectory() as tmpdir:
        output = Path(tmpdir)
        (output / "ranking").mkdir()
        (output / "ranking" / "assay_ranking.tsv").write_text("stale", encoding="utf-8")
        (output / "report.html").write_text("conservation report", encoding="utf-8")
        result = evaluate_ranking(
            object(),
            object(),
            object(),
            RankingConfig(enabled=False),
            output,
            target_name="target",
        )
        self.assertEqual(result.status, "SKIPPED")
        self.assertFalse((output / "ranking" / "assay_ranking.tsv").exists())
        self.assertEqual(
            (output / "report.html").read_text(encoding="utf-8"),
            "conservation report",
        )
        report = json.loads(result.ranking_report_path.read_text(encoding="utf-8"))
        self.assertEqual(report["status"], "SKIPPED")
```

The use of `object()` is intentional: disabled ranking must not inspect upstream scientific evidence.

Second, create stale TSV/JSON/root HTML, call enabled ranking with `make_primer_result(status="SKIPPED")`, assert `RankingError`, and then assert all three stale files are gone. This proves cleanup happens after config validation but before enabled evidence validation.

- [ ] **Step 2: Write RED tests for complete and zero-assay artifacts**

For one valid assay, assert:

- `ranking/assay_ranking.tsv`, `ranking/ranking_report.json`, and root `report.html` exist;
- TSV header includes `rank`, `assay_id`, `region_id`, `classification`, `score_status`, `final_score`, all five component columns, inclusivity counts, compatible-hit count, plausible/detectable counts, pair penalty, and `reason_codes`;
- JSON has `schema_version: 1`, effective config, class/score-state counts, ordered assays, components, structured reasons/evidence, and artifact paths;
- JSON assay ordering matches TSV ranking ordering;
- final HTML equals the newly rendered assay report, not a stale conservation report.

For a `PrimerDesignResult(status="COMPLETE", assays=(), candidates=())`, assert a header-only TSV, zero-count JSON, and HTML containing `No assay candidates`.

- [ ] **Step 3: Run artifact tests and verify RED**

```bash
python -m pytest tests/test_ranking_artifacts.py -q
```

Expected: failures because `RankingResult` and `evaluate_ranking` are not implemented.

- [ ] **Step 4: Implement artifact paths, atomic text writes, TSV, and JSON serialization**

In `ranking.py` add:

```python
def _artifact_paths(output_dir: Path) -> dict[str, Path]:
    return {
        "tsv": output_dir / "ranking" / "assay_ranking.tsv",
        "report": output_dir / "ranking" / "ranking_report.json",
        "html": output_dir / "report.html",
    }


def _atomic_write_text(destination: Path, content: str) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.parent / f".{destination.name}.{uuid.uuid4().hex}.tmp"
    try:
        with temporary.open("x", encoding="utf-8", newline="\n") as handle:
            handle.write(content)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
```

Use a fixed TSV column tuple and `json.dumps(..., indent=2, sort_keys=True) + "\n"`. Serialize `RankingReason.evidence` as a JSON object with deterministic key order. Round only presentation values in TSV/HTML; preserve full float precision in JSON.

- [ ] **Step 5: Implement disabled and enabled publication order exactly**

Start `evaluate_ranking()` with `validate_ranking_config(config)`.

Disabled path:

```python
paths["report"].parent.mkdir(parents=True, exist_ok=True)
paths["tsv"].unlink(missing_ok=True)
# Intentionally do not touch paths["html"].
_atomic_write_text(paths["report"], skipped_json)
return RankingResult("SKIPPED", (), None, paths["report"], None)
```

Enabled path, immediately after config validation:

```python
paths["report"].parent.mkdir(parents=True, exist_ok=True)
for key in ("tsv", "report", "html"):
    paths[key].unlink(missing_ok=True)
```

Then call `rank_assays`, render TSV text, JSON text, and `render_assay_report_html(...)` fully in memory. Only after all three strings are ready, atomically write TSV, HTML, then JSON report last. Return `RankingResult("COMPLETE", assays, paths["tsv"], paths["report"], paths["html"])`.

- [ ] **Step 6: Add deterministic repeated-run and full-retention tests**

Run the same valid ranking inputs into two clean output directories and assert TSV and JSON file contents are byte-identical after replacing no path-dependent values. Since artifact paths are relative, they should already match.

Create specificity evidence with `total_hit_count=25`, `retained_hit_count=1`, `truncated=True`, no amplicon, and assert the specificity score uses `25`, reaching the `0.80` floor. Do not derive isolated-hit score from `SpecificityResult.hits` length.

- [ ] **Step 7: Run all ranking stage and renderer tests**

```bash
python -m pytest \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py \
  tests/test_assay_report_html.py \
  tests/test_ranking_artifacts.py -q
```

Expected: PASS.

- [ ] **Step 8: Commit artifact publication**

```bash
git add qpcr_pipeline/ranking.py tests/test_ranking_artifacts.py
git commit -m "feat: publish ranking artifacts"
```

---

### Task 6: Integrate ranking into the pipeline and document the final contract

**Files:**
- Modify: `qpcr_pipeline/pipeline.py`
- Modify: `README.md`
- Create: `tests/test_pipeline_ranking.py`

**Interfaces:**
- `run_pipeline()` calls `evaluate_ranking(...)` after `evaluate_specificity(...)` and before `qc_report.json` is created.
- Adds exact `qc_report.json["ranking"]` keys:

```text
status
assay_count
in_silico_pass_count
review_count
high_risk_count
complete_score_count
incomplete_score_count
top_recommended_assay_id
```

- [ ] **Step 1: Write RED pipeline summary tests**

Create `tests/test_pipeline_ranking.py` following the patch style in `tests/test_pipeline_specificity.py`.

Default config assertion:

```python
self.assertEqual(
    qc["ranking"],
    {
        "status": "SKIPPED",
        "assay_count": 0,
        "in_silico_pass_count": 0,
        "review_count": 0,
        "high_risk_count": 0,
        "complete_score_count": 0,
        "incomplete_score_count": 0,
        "top_recommended_assay_id": None,
    },
)
```

For enabled ranking, patch clustering/alignment/conservation/primer design/inclusivity/specificity to return valid in-memory public result objects from `tests/ranking_fixtures.py`, leave ranking itself unpatched, run the pipeline, and assert one PASS assay is reported and `top_recommended_assay_id == "a1"`.

Add a second enabled case where the only assay is REVIEW and assert `top_recommended_assay_id is None`; REVIEW must never be silently promoted.

- [ ] **Step 2: Run pipeline ranking tests and verify RED**

```bash
python -m pytest tests/test_pipeline_ranking.py -q
```

Expected: `qc_report.json` lacks `ranking` and/or ranking is never invoked.

- [ ] **Step 3: Call ranking after specificity and construct the QC summary**

Modify `qpcr_pipeline/pipeline.py`:

```python
from qpcr_pipeline.ranking import evaluate_ranking
```

Immediately after specificity:

```python
ranking = evaluate_ranking(
    primer_design,
    inclusivity,
    specificity,
    config.ranking,
    output_dir,
    target_name=config.target_name,
)
```

Add summary generation using the final ordered `ranking.assays`:

```python
"ranking": {
    "status": ranking.status,
    "assay_count": len(ranking.assays),
    "in_silico_pass_count": sum(
        item.classification == "IN SILICO PASS" for item in ranking.assays
    ),
    "review_count": sum(item.classification == "REVIEW" for item in ranking.assays),
    "high_risk_count": sum(
        item.classification == "HIGH_RISK" for item in ranking.assays
    ),
    "complete_score_count": sum(
        item.score_status == "COMPLETE" for item in ranking.assays
    ),
    "incomplete_score_count": sum(
        item.score_status == "INCOMPLETE" for item in ranking.assays
    ),
    "top_recommended_assay_id": next(
        (
            item.assay_id
            for item in ranking.assays
            if item.classification == "IN SILICO PASS"
        ),
        None,
    ),
},
```

- [ ] **Step 4: Run pipeline and prior-stage regressions**

```bash
python -m pytest \
  tests/test_pipeline_ranking.py \
  tests/test_pipeline_specificity.py \
  tests/test_specificity_artifacts.py -q
```

Expected: PASS.

- [ ] **Step 5: Document ranking semantics and report ownership in README**

Add a new `## Classificação e ranking final dos assays` section after specificity. Include the exact YAML defaults, class precedence, the three specificity risk cases, inclusivity thresholds, five score components and default weights, the rule that score only sorts within class, structured reason codes, incomplete-score behavior, and artifacts:

```text
ranking/assay_ranking.tsv
ranking/ranking_report.json
report.html
```

State explicitly:

```text
Com ranking desabilitado, o estágio não altera report.html; portanto o relatório
publicado pela conservação continua disponível. Com ranking habilitado, o ranking
assume a propriedade de report.html e substitui o relatório de conservação pelo
relatório final consolidado dos assays.
```

Also update the conservation `report.html` paragraph to note that enabled final ranking later replaces that root file during the same run.

Document that accepted IUPAC proposals remain contextual evidence, that no proposed degenerate assay is ranked or re-tested for specificity, and that `IN SILICO PASS` is an in-silico classification rather than experimental validation.

- [ ] **Step 6: Commit pipeline integration and docs**

```bash
git add qpcr_pipeline/pipeline.py tests/test_pipeline_ranking.py README.md
git commit -m "feat: integrate final assay ranking"
```

---

### Task 7: Final regression, direct review, and completion evidence

**Files:**
- Review all files changed by Tasks 1-6.
- Do not change `.circleci/config.yml` unless a discovered defect specifically requires it; the current design does not.

**Interfaces:**
- Produces a feature branch ready for `superpowers:finishing-a-development-branch`.
- Completion evidence must include fresh local test output and later the normal CircleCI run on merged `develop`.

- [ ] **Step 1: Run the complete normal test suite**

```bash
python -m pytest -q
```

Expected: all tests PASS.

- [ ] **Step 2: Run repository whitespace/diff sanity checks**

```bash
git diff --check develop...HEAD
git status --short
git diff --stat develop...HEAD
```

Expected: `git diff --check` exits 0; status contains only intentionally preserved unrelated local files plus no uncommitted ranking implementation files.

- [ ] **Step 3: Perform a direct issue-#10 diff review**

Inspect:

```bash
git diff develop...HEAD -- \
  qpcr_pipeline/config.py \
  qpcr_pipeline/ranking.py \
  qpcr_pipeline/assay_report_html.py \
  qpcr_pipeline/pipeline.py \
  tests/test_ranking_config.py \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py \
  tests/test_ranking_artifacts.py \
  tests/test_assay_report_html.py \
  tests/test_pipeline_ranking.py \
  README.md
```

Check every acceptance requirement explicitly:

1. three classes are generated from explicit configurable rules;
2. reasons/classification are finalized before score calculation;
3. class-first ordering prevents HIGH_RISK from outranking PASS/REVIEW;
4. final score is decomposed into five named components;
5. all non-recommended candidates stay present with structured reason codes;
6. root final HTML contains F/Probe/R, degeneracy, inclusivity, specificity, class, score, and rank;
7. tests exercise class rules and ordering;
8. disabled ranking does not erase conservation `report.html`.

There is no independent subagent reviewer in this environment, so record that this was a direct diff review rather than claiming an independent review.

- [ ] **Step 4: Re-run focused issue-#10 tests after review**

```bash
python -m pytest \
  tests/test_ranking_config.py \
  tests/test_ranking_classification.py \
  tests/test_ranking_scoring.py \
  tests/test_ranking_ordering.py \
  tests/test_ranking_artifacts.py \
  tests/test_assay_report_html.py \
  tests/test_pipeline_ranking.py -q
```

Expected: PASS from fresh execution after review.

- [ ] **Step 5: Commit any review fix only after a failing regression test demonstrates it**

If review finds a defect, follow RED → GREEN for that defect, then commit only the regression test and minimal fix. If review finds no defect, make no empty commit.

- [ ] **Step 6: Hand off to branch finishing**

Invoke `superpowers:verification-before-completion`, then `superpowers:finishing-a-development-branch`. Do not claim issue #10 complete until the merged `develop` SHA receives a successful normal CircleCI run. After that evidence exists, update issue #10 acceptance checkboxes, add the verification SHA/run reference, and close the issue as completed.
