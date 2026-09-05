import io
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml

from qpcr_pipeline.cli import main


class GuidedCliTests(unittest.TestCase):
    def test_guided_prepare_writes_loadable_proposal_configuration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workspace = Path(tmpdir)
            stdout = io.StringIO()
            with patch(
                "sys.argv",
                [
                    "qpcr-pipeline",
                    "guided",
                    "prepare",
                    "--target",
                    "West Nile virus",
                    "--workspace",
                    str(workspace),
                ],
            ), patch("sys.stdout", stdout):
                code = main()

            self.assertEqual(code, 0)
            proposal_path = workspace / "config-proposal.yaml"
            self.assertTrue(proposal_path.is_file())
            config = yaml.safe_load(proposal_path.read_text(encoding="utf-8"))
            self.assertEqual(config["target"]["name"], "West Nile virus")
            self.assertIn("ncbi", config["input"])
            self.assertNotIn("fasta", config["input"])
            self.assertIn(str(proposal_path), stdout.getvalue())

    def test_guided_finalize_dispatches_approved_panel_to_guided_finalizer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            workspace = Path(tmpdir)
            approved_panel = workspace / "approved_panel.json"
            approved_panel.write_text("{}", encoding="utf-8")
            approved_config = workspace / "config-approved.yaml"
            approved_config.write_text("target: {name: West Nile virus}\n", encoding="utf-8")
            stdout = io.StringIO()

            with patch(
                "qpcr_pipeline.cli.finalize_guided_project",
                return_value=approved_config,
            ) as finalize, patch(
                "sys.argv",
                [
                    "qpcr-pipeline",
                    "guided",
                    "finalize",
                    "--target",
                    "West Nile virus",
                    "--workspace",
                    str(workspace),
                    "--approved-panel",
                    str(approved_panel),
                ],
            ), patch("sys.stdout", stdout):
                code = main()

            self.assertEqual(code, 0)
            finalize.assert_called_once_with(
                "West Nile virus",
                approved_panel,
                workspace,
            )
            self.assertIn(str(approved_config), stdout.getvalue())


if __name__ == "__main__":
    unittest.main()
