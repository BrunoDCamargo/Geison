import json
from pathlib import Path

import pytest

from qpcr_pipeline.checkpointing import (
    CheckpointInvalidity,
    CheckpointManager,
    CheckpointRequest,
    OutputIdentity,
    canonical_sha256,
    causal_fingerprint,
    result_fingerprint,
)


class DictCodec:
    def encode(self, value, outdir):
        del outdir
        return value

    def decode(self, payload, outdir):
        del outdir
        if not isinstance(payload, dict) or "value" not in payload:
            raise ValueError("invalid test state")
        return payload


def request(stage="alignment", dependency="sha256:dep"):
    return CheckpointRequest(
        stage=stage,
        dependencies={"upstream": dependency},
        inputs={"source": "sha256:input"},
        parameters={"enabled": True},
        software={"geison": "0.1.0.dev0"},
        tools={"mafft": {"version": "7.526"}},
    )


def test_canonical_sha256_is_mapping_order_independent():
    assert canonical_sha256({"b": 2, "a": 1}) == canonical_sha256({"a": 1, "b": 2})


def test_causal_fingerprint_uses_dependency_result_identity():
    first = request(dependency="sha256:a")
    second = request(dependency="sha256:b")
    assert causal_fingerprint(first) != causal_fingerprint(second)


def test_result_fingerprint_uses_state_and_sorted_outputs():
    first = result_fingerprint(
        "sha256:causal",
        "sha256:state-a",
        (OutputIdentity("a.tsv", "sha256:a"), OutputIdentity("b.json", "sha256:b")),
    )
    reordered = result_fingerprint(
        "sha256:causal",
        "sha256:state-a",
        (OutputIdentity("b.json", "sha256:b"), OutputIdentity("a.tsv", "sha256:a")),
    )
    changed = result_fingerprint(
        "sha256:causal",
        "sha256:state-b",
        (OutputIdentity("a.tsv", "sha256:a"), OutputIdentity("b.json", "sha256:b")),
    )
    assert first == reordered
    assert first != changed


def _completed_checkpoint(tmp_path: Path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    output = outdir / "alignment" / "report.json"
    output.parent.mkdir()
    output.write_text('{"status":"COMPLETE"}\n', encoding="utf-8")
    manager = CheckpointManager(outdir)
    req = request()
    codec = DictCodec()
    manager.begin(req)
    manifest = manager.complete(req, {"value": 42}, codec, (output,))
    return outdir, output, manager, req, codec, manifest


def test_complete_checkpoint_round_trips_and_validates(tmp_path):
    outdir, output, manager, req, codec, manifest = _completed_checkpoint(tmp_path)
    validation = manager.validate(req, codec)
    assert validation.valid is True
    assert validation.invalidity is None
    assert validation.loaded is not None
    assert validation.loaded.state == {"value": 42}
    assert validation.loaded.manifest == manifest
    assert manifest.status == "COMPLETE"
    assert manifest.outputs[0].path == output.relative_to(outdir).as_posix()


def test_running_checkpoint_is_not_reusable(tmp_path):
    outdir = tmp_path / "run"
    manager = CheckpointManager(outdir)
    req = request()
    manager.begin(req)
    validation = manager.validate(req, DictCodec())
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.NON_COMPLETE_STATUS


def test_failed_checkpoint_is_not_reusable(tmp_path):
    outdir = tmp_path / "run"
    manager = CheckpointManager(outdir)
    req = request()
    manager.begin(req)
    manager.fail(req, RuntimeError("boom"))
    validation = manager.validate(req, DictCodec())
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.NON_COMPLETE_STATUS


def test_missing_output_is_invalid(tmp_path):
    _, output, manager, req, codec, _ = _completed_checkpoint(tmp_path)
    output.unlink()
    validation = manager.validate(req, codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.MISSING_OUTPUT


def test_modified_output_hash_is_invalid(tmp_path):
    _, output, manager, req, codec, _ = _completed_checkpoint(tmp_path)
    output.write_text("changed\n", encoding="utf-8")
    validation = manager.validate(req, codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.OUTPUT_HASH_MISMATCH


def test_modified_state_hash_is_invalid(tmp_path):
    outdir, _, manager, req, codec, _ = _completed_checkpoint(tmp_path)
    state_path = outdir / ".checkpoints" / req.stage / "state.json"
    state_path.write_text('{"value":43}\n', encoding="utf-8")
    validation = manager.validate(req, codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.STATE_HASH_MISMATCH


def test_malformed_state_is_invalid(tmp_path):
    outdir, _, manager, req, codec, _ = _completed_checkpoint(tmp_path)
    state_path = outdir / ".checkpoints" / req.stage / "state.json"
    malformed = json.dumps({"wrong": True}) + "\n"
    state_path.write_text(malformed, encoding="utf-8")
    manifest_path = outdir / ".checkpoints" / req.stage / "manifest.json"
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    from qpcr_pipeline.checkpointing import file_sha256
    raw["state"]["sha256"] = file_sha256(state_path)
    raw["result_fingerprint"] = result_fingerprint(
        raw["fingerprint"], raw["state"]["sha256"], tuple(OutputIdentity(**item) for item in raw["outputs"])
    )
    manifest_path.write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validation = manager.validate(req, codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.INVALID_STATE


def test_modified_result_fingerprint_is_invalid(tmp_path):
    outdir, _, manager, req, codec, _ = _completed_checkpoint(tmp_path)
    manifest_path = outdir / ".checkpoints" / req.stage / "manifest.json"
    raw = json.loads(manifest_path.read_text(encoding="utf-8"))
    raw["result_fingerprint"] = "sha256:" + "0" * 64
    manifest_path.write_text(json.dumps(raw, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    validation = manager.validate(req, codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.RESULT_FINGERPRINT_MISMATCH


def test_changed_request_fingerprint_is_invalid(tmp_path):
    _, _, manager, _, codec, _ = _completed_checkpoint(tmp_path)
    validation = manager.validate(request(dependency="sha256:new"), codec)
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.FINGERPRINT_MISMATCH


def test_missing_manifest_is_invalid(tmp_path):
    manager = CheckpointManager(tmp_path / "run")
    validation = manager.validate(request(), DictCodec())
    assert validation.valid is False
    assert validation.invalidity == CheckpointInvalidity.MISSING_MANIFEST


def test_output_outside_outdir_is_rejected(tmp_path):
    outdir = tmp_path / "run"
    outdir.mkdir()
    external = tmp_path / "outside.txt"
    external.write_text("outside\n", encoding="utf-8")
    manager = CheckpointManager(outdir)
    req = request()
    manager.begin(req)
    with pytest.raises(ValueError, match="inside the pipeline output directory"):
        manager.complete(req, {"value": 1}, DictCodec(), (external,))
