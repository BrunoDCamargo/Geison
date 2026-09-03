import hashlib
import json
from dataclasses import replace
from pathlib import Path

import pytest
import yaml

from qpcr_pipeline.panel_manifest import (
    approve_panel_proposal,
    load_approved_panel_manifest,
    load_panel_proposal,
    materialize_approved_panel,
    proposal_semantic_sha256,
    write_approved_panel_manifest,
)

FIXTURE = Path("tests/fixtures/panels/west_nile_proposal.yaml")


def _proposal_payload():
    return yaml.safe_load(FIXTURE.read_text(encoding="utf-8"))


def _write_yaml(path, payload):
    path.write_text(
        yaml.safe_dump(payload, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )


def test_load_proposal_validates_schema_and_domain():
    proposal = load_panel_proposal(FIXTURE)
    assert proposal.schema_version == 1
    assert proposal.status == "PROPOSED"
    assert proposal.definition.target.name == "West Nile virus"
    assert proposal.definition.non_targets[0].criticality == "CRITICAL"


def test_approval_is_byte_deterministic(tmp_path):
    first = tmp_path / "first.json"
    second = tmp_path / "second.json"
    approve_panel_proposal(FIXTURE, first)
    approve_panel_proposal(FIXTURE, second)
    assert first.read_bytes() == second.read_bytes()


def test_semantic_hash_ignores_yaml_comments_and_formatting(tmp_path):
    original = load_panel_proposal(FIXTURE)
    reformatted = tmp_path / "reformatted.yaml"
    reformatted.write_text(
        "# review comment\n" + FIXTURE.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    assert proposal_semantic_sha256(original) == proposal_semantic_sha256(
        load_panel_proposal(reformatted)
    )


def test_approved_manifest_rejects_extra_fields(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["unexpected"] = True
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved panel manifest fields"):
        load_approved_panel_manifest(approved)


def test_approved_manifest_rejects_invalid_proposal_hash(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["proposal_sha256"] = "sha256:not-a-digest"
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="proposal_sha256"):
        load_approved_panel_manifest(approved)


@pytest.mark.parametrize("schema_version", [0, 2, True, "1"])
def test_proposal_rejects_wrong_schema_version(tmp_path, schema_version):
    payload = _proposal_payload()
    payload["schema_version"] = schema_version
    path = tmp_path / "proposal.yaml"
    _write_yaml(path, payload)
    with pytest.raises(ValueError, match="proposal schema_version"):
        load_panel_proposal(path)


def test_proposal_rejects_wrong_status(tmp_path):
    payload = _proposal_payload()
    payload["status"] = "APPROVED"
    path = tmp_path / "proposal.yaml"
    _write_yaml(path, payload)
    with pytest.raises(ValueError, match="proposal status"):
        load_panel_proposal(path)


def test_proposal_rejects_extra_nested_fields(tmp_path):
    payload = _proposal_payload()
    payload["definition"]["target"]["unexpected"] = True
    path = tmp_path / "proposal.yaml"
    _write_yaml(path, payload)
    with pytest.raises(ValueError, match="panel target fields"):
        load_panel_proposal(path)


@pytest.mark.parametrize("schema_version", [0, 2, True, "1"])
def test_approved_manifest_rejects_wrong_schema_version(
    tmp_path,
    schema_version,
):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["schema_version"] = schema_version
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved panel schema_version"):
        load_approved_panel_manifest(approved)


def test_approved_manifest_rejects_wrong_status(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["status"] = "PROPOSED"
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved panel status"):
        load_approved_panel_manifest(approved)


@pytest.mark.parametrize("approved_by_user", [False, 1, "true", None])
def test_approved_by_user_must_be_exactly_true(tmp_path, approved_by_user):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["approved_by_user"] = approved_by_user
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved_by_user.*true"):
        load_approved_panel_manifest(approved)


@pytest.mark.parametrize(
    "field",
    [
        "schema_version",
        "status",
        "approved_by_user",
        "proposal_sha256",
        "definition",
    ],
)
def test_approved_manifest_rejects_missing_fields(tmp_path, field):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    del payload[field]
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="approved panel manifest fields"):
        load_approved_panel_manifest(approved)


def test_approved_manifest_rejects_invalid_nested_panel_content(tmp_path):
    approved = tmp_path / "approved.json"
    approve_panel_proposal(FIXTURE, approved)
    payload = json.loads(approved.read_text(encoding="utf-8"))
    payload["definition"]["target"]["groups"] = []
    approved.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="target groups.*non-empty"):
        load_approved_panel_manifest(approved)


def test_failed_serialization_preserves_existing_output(tmp_path):
    output = tmp_path / "approved.json"
    manifest = approve_panel_proposal(FIXTURE, tmp_path / "valid.json")
    output.write_bytes(b"existing\n")

    with pytest.raises(ValueError, match="approved_by_user.*true"):
        write_approved_panel_manifest(
            replace(manifest, approved_by_user=False),
            output,
        )

    assert output.read_bytes() == b"existing\n"


def test_materialize_approved_panel_writes_canonical_copy_and_result(tmp_path):
    source = tmp_path / "source.json"
    approve_panel_proposal(FIXTURE, source)

    result = materialize_approved_panel(source, tmp_path / "run")

    expected_path = tmp_path / "run" / "panel" / "approved_panel.json"
    expected_digest = hashlib.sha256(expected_path.read_bytes()).hexdigest()
    assert result.status == "APPROVED"
    assert result.manifest_path == expected_path
    assert result.manifest_sha256 == f"sha256:{expected_digest}"
    assert result.target_mode == "broad_detection"
    assert result.non_target_count == 3
    assert expected_path.read_bytes() == source.read_bytes()
