from __future__ import annotations

from pathlib import Path

from qpcr_pipeline import checkpoint_stages
from qpcr_pipeline.checkpointing import causal_fingerprint
from qpcr_pipeline.config import PipelineConfig


class NoToolIdentityProvider:
    def identity(self, tool_name: str):
        raise AssertionError(f"unexpected tool identity request: {tool_name}")


def test_stage_fingerprint_changes_when_geison_source_changes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    config = PipelineConfig(
        target_name="target",
        input_fasta=tmp_path / "target.fasta",
    )
    monkeypatch.setattr(checkpoint_stages, "geison_version", lambda: "0.1.0.dev0")
    monkeypatch.setattr(
        checkpoint_stages,
        "geison_source_sha256",
        lambda: "sha256:" + "a" * 64,
    )
    first = checkpoint_stages.stage_request(
        "panel",
        config,
        {},
        {},
        NoToolIdentityProvider(),
    )

    monkeypatch.setattr(
        checkpoint_stages,
        "geison_source_sha256",
        lambda: "sha256:" + "b" * 64,
    )
    second = checkpoint_stages.stage_request(
        "panel",
        config,
        {},
        {},
        NoToolIdentityProvider(),
    )

    assert first.software["geison"] == {
        "version": "0.1.0.dev0",
        "source_sha256": "sha256:" + "a" * 64,
    }
    assert second.software["geison"] == {
        "version": "0.1.0.dev0",
        "source_sha256": "sha256:" + "b" * 64,
    }
    assert causal_fingerprint(first) != causal_fingerprint(second)
