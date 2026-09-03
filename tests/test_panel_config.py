from copy import deepcopy
from pathlib import Path

import pytest
import yaml

from qpcr_pipeline.config import load_config


def _base_payload(tmp_path):
    fasta = tmp_path / "target.fasta"
    fasta.write_text(">s1\nACGT\n", encoding="utf-8")
    return {
        "target": {"name": "West Nile virus"},
        "input": {"fasta": fasta.as_posix()},
        "panel": {
            "proposal": {
                "target": {
                    "name": "West Nile virus",
                    "taxid": None,
                    "mode": "broad_detection",
                    "subtype": None,
                    "groups": [
                        {
                            "name": "lineage_1",
                            "required": True,
                            "dataset_roles": ["DESIGN", "CHALLENGE"],
                            "reasons": ["target_diversity"],
                            "proposed_by": ["manual"],
                            "sequence_selection": [],
                        }
                    ],
                },
                "non_targets": [
                    {
                        "name": "Usutu virus",
                        "taxid": None,
                        "criticality": "CRITICAL",
                        "dataset_roles": ["DESIGN", "CHALLENGE"],
                        "reasons": ["phylogenetic_neighbor"],
                        "proposed_by": ["manual"],
                        "sequence_selection": [],
                    }
                ],
                "diagnostic_context": {
                    "syndrome": "arboviral febrile disease",
                    "geography": "Brazil",
                    "sample_type": "human serum",
                    "vector": "mosquito",
                },
            }
        },
    }


def _load_payload(tmp_path, payload):
    path = tmp_path / "config.yaml"
    path.write_text(
        yaml.safe_dump(payload, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    return load_config(path)


def test_loads_inline_panel_proposal(tmp_path):
    config = _load_payload(tmp_path, _base_payload(tmp_path))
    assert config.panel is not None
    assert config.panel.proposal is not None
    assert config.panel.frozen_manifest is None
    assert config.panel.proposal.target.mode == "broad_detection"


def test_primer_design_requires_panel(tmp_path):
    payload = _base_payload(tmp_path)
    del payload["panel"]
    payload["alignment"] = {"enabled": True}
    payload["conservation"] = {"enabled": True}
    payload["primer_design"] = {"enabled": True}
    with pytest.raises(
        ValueError,
        match="Enabled primer design requires a panel proposal or frozen manifest",
    ):
        _load_payload(tmp_path, payload)


def test_frozen_manifest_parses_to_path(tmp_path):
    payload = _base_payload(tmp_path)
    payload["panel"] = {"frozen_manifest": "reviewed/approved-panel.json"}
    config = _load_payload(tmp_path, payload)
    assert config.panel is not None
    assert config.panel.proposal is None
    assert config.panel.frozen_manifest == Path("reviewed/approved-panel.json")


def test_panel_rejects_proposal_and_frozen_manifest_together(tmp_path):
    payload = _base_payload(tmp_path)
    payload["panel"]["frozen_manifest"] = "approved.json"
    with pytest.raises(ValueError, match="exactly one proposal or frozen_manifest"):
        _load_payload(tmp_path, payload)


def test_panel_section_rejects_neither_mode(tmp_path):
    payload = _base_payload(tmp_path)
    payload["panel"] = {}
    with pytest.raises(ValueError, match="exactly one proposal or frozen_manifest"):
        _load_payload(tmp_path, payload)


def test_panel_section_must_be_a_mapping(tmp_path):
    payload = _base_payload(tmp_path)
    payload["panel"] = True
    with pytest.raises(ValueError, match="section 'panel'.*mapping"):
        _load_payload(tmp_path, payload)


@pytest.mark.parametrize(
    ("location", "message"),
    [
        ("panel", "section 'panel'.*extra.*unrecognized"),
        ("proposal", "panel proposal.*extra.*unrecognized"),
        ("target", "panel target.*extra.*unrecognized"),
        ("target_group", "target group 1.*extra.*unrecognized"),
        ("non_target", "non-target 1.*extra.*unrecognized"),
        ("context", "diagnostic_context.*extra.*unrecognized"),
        ("selection", "sequence_selection entry 1.*extra.*unrecognized"),
    ],
)
def test_panel_config_rejects_unknown_fields_at_every_level(
    tmp_path,
    location,
    message,
):
    payload = deepcopy(_base_payload(tmp_path))
    proposal = payload["panel"]["proposal"]
    locations = {
        "panel": payload["panel"],
        "proposal": proposal,
        "target": proposal["target"],
        "target_group": proposal["target"]["groups"][0],
        "non_target": proposal["non_targets"][0],
        "context": proposal["diagnostic_context"],
    }
    if location == "selection":
        selection = {
            "dataset_role": "DESIGN",
            "method": "manual_fixture",
            "source": "unit-test",
            "details": [],
            "extra": True,
        }
        proposal["target"]["groups"][0]["sequence_selection"] = [selection]
    else:
        locations[location]["extra"] = True
    with pytest.raises(ValueError, match=message):
        _load_payload(tmp_path, payload)


def test_inline_target_name_must_match_pipeline_target_after_stripping(tmp_path):
    payload = _base_payload(tmp_path)
    payload["target"]["name"] = "Zika virus"
    with pytest.raises(ValueError, match="Inline panel target name must match"):
        _load_payload(tmp_path, payload)


def test_inline_target_name_allows_only_surrounding_whitespace_difference(tmp_path):
    payload = _base_payload(tmp_path)
    payload["target"]["name"] = "  West Nile virus  "
    config = _load_payload(tmp_path, payload)
    assert config.panel is not None


def test_panel_free_primer_design_disabled_config_remains_valid(tmp_path):
    payload = _base_payload(tmp_path)
    del payload["panel"]
    config = _load_payload(tmp_path, payload)
    assert config.panel is None
    assert config.primer_design.enabled is False


def test_panel_list_fields_reject_non_lists(tmp_path):
    payload = _base_payload(tmp_path)
    payload["panel"]["proposal"]["target"]["groups"][0][
        "dataset_roles"
    ] = "DESIGN"
    with pytest.raises(ValueError, match="dataset_roles.*list"):
        _load_payload(tmp_path, payload)
