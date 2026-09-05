import hashlib
import json
import tempfile
import unittest
from pathlib import Path

import yaml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from qpcr_pipeline.config import load_config
from qpcr_pipeline.guided import (
    KNOWLEDGE_VERSION,
    build_guided_proposal_config,
    finalize_guided_project,
    supported_guided_targets,
)
from qpcr_pipeline.ncbi import NcbiFetchedRecord, ResolvedNcbiQuery
from qpcr_pipeline.panel_manifest import approve_panel_proposal, prepare_panel_preflight


class FakeNcbiClient:
    def __init__(self, query_uids, records_by_uid):
        self.query_uids = query_uids
        self.records_by_uid = records_by_uid
        self.resolve_calls = []
        self.fetch_calls = []

    def resolve_query(self, query, max_records):
        self.resolve_calls.append((query, max_records))
        uids = tuple(self.query_uids[query][:max_records])
        return ResolvedNcbiQuery(
            uids=uids,
            reported_count=len(self.query_uids[query]),
            query_translation=query,
        )

    def fetch_records(self, identifiers, *, identifier_kind):
        self.fetch_calls.append((tuple(identifiers), identifier_kind))
        return tuple(
            NcbiFetchedRecord(request_id=identifier, record=self.records_by_uid[identifier])
            for identifier in identifiers
        )


def record(accession: str) -> SeqRecord:
    value = SeqRecord(Seq("ATGC" * 2800), id=accession, description=accession)
    value.annotations["molecule_type"] = "RNA"
    return value


class GuidedConfigTests(unittest.TestCase):
    def test_west_nile_proposal_is_ncbi_first_and_contains_reviewable_panel(self):
        proposal = build_guided_proposal_config("West Nile virus")

        self.assertEqual(supported_guided_targets(), ("West Nile virus",))
        self.assertEqual(KNOWLEDGE_VERSION, "2026-09-05")
        self.assertEqual(proposal["target"]["name"], "West Nile virus")
        self.assertEqual(
            proposal["input"]["ncbi"]["query"],
            '"West Nile virus"[Organism] AND complete genome[Title]',
        )
        self.assertEqual(proposal["input"]["ncbi"]["max_records"], 50)
        self.assertFalse(proposal["contrastive_conservation"]["enabled"])
        self.assertNotIn("fasta", proposal["input"])

        non_targets = proposal["panel"]["proposal"]["non_targets"]
        self.assertEqual(
            [(item["name"], item["criticality"]) for item in non_targets],
            [
                ("Usutu virus", "CRITICAL"),
                ("Japanese encephalitis virus", "CRITICAL"),
                ("Dengue virus", "IMPORTANT"),
            ],
        )
        self.assertTrue(
            all("geison_guided_knowledge@2026-09-05" in item["proposed_by"] for item in non_targets)
        )
        self.assertEqual(
            [item["name"] for item in proposal["off_targets"]],
            ["Usutu virus", "Japanese encephalitis virus", "Dengue virus"],
        )
        self.assertTrue(
            all("frozen_dataset" in item and "fasta" not in item for item in proposal["off_targets"])
        )

    def test_unsupported_target_fails_without_inventing_a_panel(self):
        with self.assertRaisesRegex(ValueError, "Supported guided targets: West Nile virus"):
            build_guided_proposal_config("Unknown virus")

    def test_finalize_freezes_approved_challenges_and_writes_standard_config(self):
        queries = {
            '"Usutu virus"[Organism] AND complete genome[Title]': ("u1", "u2"),
            '"Japanese encephalitis virus"[Organism] AND complete genome[Title]': ("j1",),
            '"Dengue virus"[Organism] AND complete genome[Title]': ("d1",),
        }
        client = FakeNcbiClient(
            queries,
            {
                "u1": record("USUTU_1.1"),
                "u2": record("USUTU_2.1"),
                "j1": record("JEV_1.1"),
                "d1": record("DENV_1.1"),
            },
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            workspace = Path(tmpdir)
            proposal = build_guided_proposal_config("West Nile virus")
            config_path = workspace / "config-proposal.yaml"
            config_path.write_text(yaml.safe_dump(proposal, sort_keys=False), encoding="utf-8")
            parsed = load_config(config_path)
            preflight = prepare_panel_preflight(
                parsed.panel,
                workspace / "output",
                target_name=parsed.target_name,
            )
            self.assertEqual(preflight.status, "ACTION_REQUIRED")
            assert preflight.proposal_path is not None
            approved_path = workspace / "approved_panel.json"
            approve_panel_proposal(preflight.proposal_path, approved_path)

            approved_config_path = finalize_guided_project(
                "West Nile virus",
                approved_path,
                workspace,
                ncbi_client=client,
            )

            self.assertEqual(approved_config_path, workspace / "config-approved.yaml")
            approved = yaml.safe_load(approved_config_path.read_text(encoding="utf-8"))
            self.assertEqual(
                approved["input"]["ncbi"]["query"],
                '"West Nile virus"[Organism] AND complete genome[Title]',
            )
            self.assertTrue(approved["contrastive_conservation"]["enabled"])
            self.assertEqual(approved["panel"]["frozen_manifest"], str(approved_path))
            self.assertEqual(
                [item["name"] for item in approved["off_targets"]],
                ["Usutu virus", "Japanese encephalitis virus", "Dengue virus"],
            )
            for index, item in enumerate(approved["off_targets"], 1):
                frozen = Path(item["frozen_dataset"])
                self.assertTrue(frozen.is_dir())
                self.assertTrue((frozen / "dataset_manifest.json").is_file())
                self.assertTrue((frozen / "records.gb").is_file())
                self.assertIn(f"guided_challenges/{index:03d}-", frozen.as_posix())

            acquisition_path = workspace / "guided_acquisition_manifest.json"
            acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
            self.assertEqual(acquisition["knowledge_version"], KNOWLEDGE_VERSION)
            self.assertEqual(
                acquisition["approved_panel_sha256"],
                hashlib.sha256(approved_path.read_bytes()).hexdigest(),
            )
            self.assertEqual([x["record_count"] for x in acquisition["datasets"]], [2, 1, 1])
            serialized = json.dumps(acquisition)
            self.assertNotIn("NCBI_EMAIL", serialized)
            self.assertNotIn("NCBI_API_KEY", serialized)
            self.assertNotIn("@example", serialized)

        self.assertEqual(
            client.resolve_calls,
            [
                ('"Usutu virus"[Organism] AND complete genome[Title]', 20),
                ('"Japanese encephalitis virus"[Organism] AND complete genome[Title]', 20),
                ('"Dengue virus"[Organism] AND complete genome[Title]', 20),
            ],
        )

    def test_finalize_rejects_approved_target_mismatch_before_network(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workspace = Path(tmpdir)
            proposal = build_guided_proposal_config("West Nile virus")
            proposal["panel"]["proposal"]["target"]["name"] = "Different virus"
            config_path = workspace / "config.yaml"
            config_path.write_text(yaml.safe_dump(proposal, sort_keys=False), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "panel target name must match"):
                load_config(config_path)


if __name__ == "__main__":
    unittest.main()
