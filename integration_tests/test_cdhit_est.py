import json
import shutil
import tempfile
import unittest
from pathlib import Path

from qpcr_pipeline.clustering import cluster_sequences
from qpcr_pipeline.config import ClusteringConfig
from qpcr_pipeline.local_input import LocalSequenceRecord
from qpcr_pipeline.models import EvaluationSet


@unittest.skipUnless(shutil.which("cd-hit-est"), "cd-hit-est is not installed")
class CdHitEstIntegrationTests(unittest.TestCase):
    def test_real_cdhit_clusters_similar_sequences_and_publishes_artifacts(self):
        sequence_one = "ACGT" * 100
        records = (
            LocalSequenceRecord("seq-1", sequence_one),
            LocalSequenceRecord("seq-2", sequence_one),
            LocalSequenceRecord("seq-3", "T" * 400),
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            result = cluster_sequences(
                records,
                EvaluationSet(("seq-1", "seq-2", "seq-3")),
                ClusteringConfig(enabled=True, identity=0.95),
                Path(tmpdir),
            )

            self.assertEqual(result.discovery_set.sequence_ids, ("seq-1", "seq-3"))
            self.assertEqual(
                {member.sequence_id for member in result.clusters[0].members},
                {"seq-1", "seq-2"},
            )
            self.assertEqual(
                json.loads(result.report_path.read_text())["counts"],
                {"evaluation": 3, "discovery": 2},
            )
