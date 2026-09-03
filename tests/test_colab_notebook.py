import json
import unittest
from pathlib import Path


_REPO_ROOT = Path(__file__).resolve().parents[1]
_NOTEBOOK_PATH = _REPO_ROOT / "notebooks" / "geison_colab.ipynb"


class ColabNotebookBranchTests(unittest.TestCase):
    def test_official_notebook_tracks_main_branch(self):
        notebook = json.loads(_NOTEBOOK_PATH.read_text(encoding="utf-8"))
        text = "".join(
            source_line
            for cell in notebook["cells"]
            for source_line in cell.get("source", [])
        )

        self.assertNotIn("develop", text)
        self.assertIn("git -C Geison checkout main", text)
        self.assertIn("git -C Geison pull --ff-only origin main", text)
        self.assertIn("git clone --branch main", text)


if __name__ == "__main__":
    unittest.main()
