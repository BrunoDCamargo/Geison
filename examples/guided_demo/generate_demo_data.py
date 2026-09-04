#!/usr/bin/env python3
"""Generate deterministic synthetic inputs for the guided Geison Colab demo.

The data are intentionally synthetic. The generator creates files and configs
only; all scientific analysis remains in the Geison package/CLI.
"""

from __future__ import annotations

import random
import sys
from pathlib import Path

import yaml


SEED = 20260904
REFERENCE_LENGTH = 1200
UNIQUE_WINDOW_SIZE = 100
SHARED_START = 201
SHARED_END = 300
DISCRIMINANT_START = 601
DISCRIMINANT_END = 700
CONSOLIDATED_START = 61
CONSOLIDATED_END = 160

TARGET_NAME = "West Nile-like synthetic target"
CHALLENGES = (
    (
        "Related flavivirus-like synthetic A",
        "challenge-related-a.fasta",
        "CRITICAL",
        "synthetic close differential",
    ),
    (
        "Related flavivirus-like synthetic B",
        "challenge-related-b.fasta",
        "CRITICAL",
        "synthetic close differential",
    ),
    (
        "Context arbovirus-like synthetic A",
        "challenge-context-a.fasta",
        "IMPORTANT",
        "synthetic contextual differential",
    ),
)

_BASES = "ACGT"
_MUTATION_OFFSETS = (1, 2, 3)


def _mutated_base(base: str, offset: int) -> str:
    return _BASES[(_BASES.index(base) + offset) % len(_BASES)]


def _unique_reference() -> str:
    rng = random.Random(SEED)
    for _ in range(100):
        reference = "".join(rng.choice(_BASES) for _ in range(REFERENCE_LENGTH))
        windows = {
            reference[index : index + UNIQUE_WINDOW_SIZE]
            for index in range(REFERENCE_LENGTH - UNIQUE_WINDOW_SIZE + 1)
        }
        if len(windows) == REFERENCE_LENGTH - UNIQUE_WINDOW_SIZE + 1:
            return reference
    raise RuntimeError("Could not construct a unique-window synthetic reference.")


def _protected_target_position(index: int) -> bool:
    # Preserve target-side stability around both discriminant zones and the
    # shared-conserved zone used by the demonstration.
    return 0 <= index < 350 or 500 <= index < 800


def _target_variant(reference: str, variant_index: int) -> str:
    sequence = list(reference)
    allowed = [
        index
        for index in range(len(sequence))
        if not _protected_target_position(index)
    ]
    rng = random.Random(SEED + 100 * variant_index)
    for index in rng.sample(allowed, 28):
        sequence[index] = _mutated_base(sequence[index], variant_index)
    return "".join(sequence)


def _challenge_sequence(reference: str, challenge_index: int) -> str:
    sequence = list(reference)
    offset = _MUTATION_OFFSETS[challenge_index]

    # Two deterministic discriminant zones are altered base-by-base. The
    # central zone supports the focused real-engine regression. The near-start
    # zone makes several 100-base windows expand to the same 1-300 candidate
    # interval, exercising deterministic candidate consolidation end-to-end.
    for start, end in (
        (DISCRIMINANT_START, DISCRIMINANT_END),
        (CONSOLIDATED_START, CONSOLIDATED_END),
    ):
        for index in range(start - 1, end):
            sequence[index] = _mutated_base(sequence[index], offset)

    # Add deterministic background diversity without changing the shared
    # conserved demonstration interval or either deliberate discriminant zone.
    excluded = (
        set(range(SHARED_START - 1, SHARED_END))
        | set(range(DISCRIMINANT_START - 1, DISCRIMINANT_END))
        | set(range(CONSOLIDATED_START - 1, CONSOLIDATED_END))
    )
    allowed = [index for index in range(len(sequence)) if index not in excluded]
    rng = random.Random(SEED + 1000 + challenge_index)
    background_count = (18, 30, 54)[challenge_index]
    for index in rng.sample(allowed, background_count):
        sequence[index] = _mutated_base(sequence[index], offset)
    return "".join(sequence)


def _fasta_text(records: tuple[tuple[str, str], ...]) -> str:
    lines: list[str] = []
    for sequence_id, sequence in records:
        lines.append(f">{sequence_id}")
        lines.extend(
            sequence[index : index + 80]
            for index in range(0, len(sequence), 80)
        )
    return "\n".join(lines) + "\n"


def _panel_definition() -> dict[str, object]:
    return {
        "target": {
            "name": TARGET_NAME,
            "taxid": None,
            "mode": "broad_detection",
            "subtype": None,
            "groups": [
                {
                    "name": "synthetic target diversity",
                    "required": True,
                    "dataset_roles": ["DESIGN"],
                    "reasons": ["synthetic target diversity demonstration"],
                    "proposed_by": ["guided demo generator"],
                    "sequence_selection": [],
                }
            ],
        },
        "non_targets": [
            {
                "name": name,
                "taxid": None,
                "criticality": criticality,
                "dataset_roles": ["CHALLENGE"],
                "reasons": [reason],
                "proposed_by": ["guided demo generator"],
                "sequence_selection": [],
            }
            for name, _, criticality, reason in CHALLENGES
        ],
        "diagnostic_context": {
            "syndrome": "synthetic arboviral febrile scenario",
            "geography": "synthetic setting",
            "sample_type": "synthetic serum-like material",
            "vector": "synthetic mosquito-like context",
        },
    }


def _base_config() -> dict[str, object]:
    return {
        "target": {"name": TARGET_NAME},
        "input": {"fasta": "target.fasta"},
        "qc": {
            "min_length": 1000,
            "max_ambiguous_fraction": 0.0,
        },
        "alignment": {"enabled": True, "threads": 2},
        "conservation": {
            "enabled": True,
            "window_size": 100,
            "step_size": 20,
        },
        "primer_design": {
            "enabled": True,
            "max_candidate_regions": 6,
            "assays_per_region": 3,
            "candidate_region_length": 300,
            "product_size_min": 70,
            "product_size_max": 180,
            "min_mean_conservation": 0.95,
            "min_minimum_conservation": 0.90,
            "min_mean_coverage": 0.95,
            "max_mean_gap_frequency": 0.05,
            "max_mean_entropy_bits": 0.40,
            "min_usable_fraction": 0.90,
        },
        "inclusivity": {"enabled": True},
        "off_targets": [
            {"name": name, "fasta": filename}
            for name, filename, _, _ in CHALLENGES
        ],
        "specificity": {"enabled": True},
        "ranking": {"enabled": True},
    }


def _yaml_text(payload: object) -> str:
    return yaml.safe_dump(
        payload,
        sort_keys=False,
        allow_unicode=False,
        default_flow_style=False,
    )


def generate(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    reference = _unique_reference()

    target_records = (
        ("synthetic-target-reference", reference),
        ("synthetic-target-variant-1", _target_variant(reference, 1)),
        ("synthetic-target-variant-2", _target_variant(reference, 2)),
        ("synthetic-target-variant-3", _target_variant(reference, 3)),
    )
    (output_dir / "target.fasta").write_text(
        _fasta_text(target_records), encoding="utf-8"
    )

    for index, (_, filename, _, _) in enumerate(CHALLENGES):
        sequence = _challenge_sequence(reference, index)
        (output_dir / filename).write_text(
            _fasta_text(((f"synthetic-challenge-{index + 1}", sequence),)),
            encoding="utf-8",
        )

    definition = _panel_definition()
    panel_proposal = {
        "schema_version": 1,
        "status": "PROPOSED",
        "definition": definition,
    }
    (output_dir / "panel-proposal.yaml").write_text(
        _yaml_text(panel_proposal), encoding="utf-8"
    )

    base = _base_config()
    proposal = dict(base)
    proposal["panel"] = {"proposal": definition}
    proposal["contrastive_conservation"] = {"enabled": False}
    (output_dir / "config-proposal.yaml").write_text(
        _yaml_text(proposal), encoding="utf-8"
    )

    approved = dict(base)
    approved["panel"] = {"frozen_manifest": "__APPROVED_PANEL__"}
    approved["contrastive_conservation"] = {"enabled": True}
    (output_dir / "config-approved-template.yaml").write_text(
        _yaml_text(approved), encoding="utf-8"
    )


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(
            "Usage: python examples/guided_demo/generate_demo_data.py OUTPUT_DIR",
            file=sys.stderr,
        )
        return 2
    generate(Path(argv[1]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
