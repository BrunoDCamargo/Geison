import json
from pathlib import Path


NOTEBOOK_PATH = Path(__file__).parents[1] / "notebooks" / "geison_colab.ipynb"


def _load_notebook():
    assert NOTEBOOK_PATH.is_file(), "Official Colab notebook is missing"
    return json.loads(NOTEBOOK_PATH.read_text(encoding="utf-8"))


def _cell_text(notebook, cell_type):
    return "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook.get("cells", [])
        if cell.get("cell_type") == cell_type
    )


def test_official_colab_notebook_covers_issue_13_flow():
    notebook = _load_notebook()
    code = _cell_text(notebook, "code")
    markdown = _cell_text(notebook, "markdown")

    assert notebook["nbformat"] == 4
    assert "git clone" in code
    assert "--branch develop" in code
    assert "pull --ff-only origin develop" in code
    assert "python -m pip install -e" in code

    assert "apt-get" in code
    assert "cd-hit" in code
    assert "mafft" in code
    assert "primer3" in code

    assert "qpcr-pipeline doctor" in code
    assert "config.yaml" in code
    assert "qpcr-pipeline run" in code
    assert "--outdir" in code
    assert "--resume" in code
    assert "report.html" in code

    assert "git pull" in markdown.lower()
    assert "resume" in markdown.lower()
    assert "report.html" in markdown.lower()


def test_colab_notebook_delegates_scientific_work_to_geison_cli():
    notebook = _load_notebook()
    code = _cell_text(notebook, "code")

    assert "from qpcr_pipeline" not in code
    assert "import qpcr_pipeline" not in code
    assert "from Bio" not in code
    assert "import Bio" not in code
    assert "def " not in code
    assert "class " not in code
