from __future__ import annotations

from pathlib import Path

from qpcr_pipeline import checkpoint_stages
from qpcr_pipeline.config import PipelineConfig


class NoToolIdentityProvider:
    def identity(self, tool_name: str):
        raise AssertionError(f"unexpected tool identity request: {tool_name}")


def test_stage_request_records_version_and_source_hash(tmp_path: Path, monkeypatch) -> None:
    config = PipelineConfig(
        target_name="target",
        input_fasta=tmp_path / "target.fasta",
    )
    monkeypatch.setattr(checkpoint_stages, "geison_version", lambda: "0.1.0.dev0")

    request = checkpoint_stages.stage_request(
        "panel",
        config,
        {},
        {},
        NoToolIdentityProvider(),
    )

    identity = request.software["geison"]
    assert isinstance(identity, dict)
    assert identity["version"] == "0.1.0.dev0"
    assert isinstance(identity["source_sha256"], str)
    assert identity["source_sha256"].startswith("sha256:")
    assert len(identity["source_sha256"]) == len("sha256:") + 64


def test_source_tree_hash_changes_when_python_source_changes(tmp_path: Path) -> None:
    assert hasattr(checkpoint_stages, "_source_tree_sha256")
    source_tree_sha256 = getattr(checkpoint_stages, "_source_tree_sha256")

    package = tmp_path / "qpcr_pipeline"
    package.mkdir()
    first = package / "first.py"
    second = package / "nested" / "second.py"
    second.parent.mkdir()
    first.write_text("VALUE = 1\n", encoding="utf-8")
    second.write_text("VALUE = 2\n", encoding="utf-8")

    initial = source_tree_sha256(package)
    first.write_text("VALUE = 3\n", encoding="utf-8")
    changed = source_tree_sha256(package)

    assert initial.startswith("sha256:")
    assert len(initial) == len("sha256:") + 64
    assert initial != changed
