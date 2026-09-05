from __future__ import annotations

import json
from pathlib import Path


NOTEBOOK = Path("notebooks/geison_guided_colab.ipynb")
GUIDE = Path("docs/guided-colab.md")


def _notebook_text():
    notebook = json.loads(NOTEBOOK.read_text(encoding="utf-8"))
    markdown = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
        if cell.get("cell_type") == "markdown"
    )
    code = "\n".join(
        "".join(cell.get("source", []))
        for cell in notebook["cells"]
        if cell.get("cell_type") == "code"
    )
    return notebook, markdown, code


def test_guided_notebook_has_researcher_facing_sections_and_cli_flow():
    notebook, markdown, code = _notebook_text()
    assert notebook["nbformat"] == 4

    required_markdown = [
        "Geison - Assay discovery workbench",
        "Project and panel",
        "Data readiness",
        "Target conservation",
        "Target vs non-target contrast",
        "Assay design",
        "Target coverage",
        "Specificity",
        "Final candidates",
        "Reproducibility",
        "11. Researcher report",
    ]
    for marker in required_markdown:
        assert marker in markdown

    for marker in [
        "#@param",
        "Guided (NCBI)",
        "Demo (synthetic)",
        "Advanced (local sequences)",
        "qpcr-pipeline doctor",
        "qpcr-pipeline",
        "run",
        "panel",
        "approve",
        "guided",
        "prepare",
        "finalize",
        "--resume",
    ]:
        assert marker in code


def test_guided_notebook_defaults_to_low_intrusion_ncbi_mode():
    _, _, code = _notebook_text()

    assert 'mode = "Guided (NCBI)"' in code
    assert 'target_name = "West Nile virus"' in code
    assert 'os.environ["NCBI_EMAIL"]' in code
    assert '"guided", "prepare"' in code
    assert '"guided", "finalize"' in code

    # The active kernel must see the cloned source immediately; editable install
    # path hooks created by a shell cell are not guaranteed to affect it.
    assert 'geison_repo = Path("/content/Geison")' in code
    assert "sys.path.insert(0, str(geison_repo))" in code


def test_guided_notebook_keeps_local_sequence_paths_in_advanced_mode_only():
    _, _, code = _notebook_text()

    assert 'mode == "Advanced (local sequences)"' in code
    assert "project_target_fasta" in code
    assert "challenge_1_fasta" in code
    assert "challenge_2_fasta" in code

    # Guided mode must be built by Geison, rather than assembling a local FASTA
    # config in the default branch of notebook logic.
    guided_start = code.index('mode == "Guided (NCBI)"')
    advanced_start = code.index('mode == "Advanced (local sequences)"')
    guided_block = code[guided_start:advanced_start]
    assert "project_target_fasta" not in guided_block
    assert "challenge_1_fasta" not in guided_block
    assert "challenge_2_fasta" not in guided_block


def test_guided_notebook_does_not_duplicate_scientific_implementation():
    _, _, code = _notebook_text()
    for forbidden in [
        "from qpcr_pipeline.conservation",
        "from qpcr_pipeline.contrastive_conservation",
        "from qpcr_pipeline.primer_design",
        "from qpcr_pipeline.inclusivity",
        "from qpcr_pipeline.specificity",
        "from qpcr_pipeline.ranking",
        "import qpcr_pipeline.conservation",
        "import qpcr_pipeline.contrastive_conservation",
        "import qpcr_pipeline.primer_design",
        "import qpcr_pipeline.inclusivity",
        "import qpcr_pipeline.specificity",
        "import qpcr_pipeline.ranking",
        "from Bio",
        "import Bio",
        "PairwiseAligner",
        "find_plausible_amplicons",
        "design_primers(",
        "analyze_conservation(",
        "analyze_contrastive_conservation(",
    ]:
        assert forbidden not in code

    assert (
        "from qpcr_pipeline.evidence_bundle import create_evidence_bundle" in code
    )


def test_guided_notebook_exposes_states_configs_and_published_artifacts():
    _, markdown, code = _notebook_text()
    combined = markdown + "\n" + code
    for marker in [
        "run_manifest.json",
        "config-proposal.yaml",
        "config-approved.yaml",
        "ACTION_REQUIRED",
        "PARTIAL",
        "FAILED",
        "COMPLETED",
        "contrastive_conservation/window_metrics.tsv",
        "contrastive_conservation/dataset_metrics.tsv",
        "contrastive_conservation/candidate_regions.tsv",
    ]:
        assert marker in combined


def test_guided_notebook_keeps_explicit_human_approval_gate():
    _, markdown, code = _notebook_text()
    combined = markdown + "\n" + code
    assert "APROVAR" in combined
    assert "PANEL_APPROVAL_REQUIRED" in combined
    assert "approved_panel.json" in combined


def test_guided_notebook_exposes_researcher_report_and_evidence_downloads():
    _, markdown, code = _notebook_text()
    combined = markdown + "\n" + code

    for marker in [
        "Researcher report",
        "View report",
        "Download report.html",
        "Download evidence bundle.zip",
        "report.html",
        "evidence_bundle.zip",
        "create_evidence_bundle",
        "files.download",
    ]:
        assert marker in combined

    assert "google.colab import files" in code
    assert "approved_config" in code
    assert "approved_panel_path" in code


def test_guided_colab_guide_documents_normal_and_advanced_use():
    text = GUIDE.read_text(encoding="utf-8")
    for marker in [
        "Guided (NCBI)",
        "Demo (synthetic)",
        "Advanced (local sequences)",
        "NCBI_EMAIL",
        "PANEL_APPROVAL_REQUIRED",
        "APROVAR",
        "contrastive_conservation",
        "specificity",
        "run_manifest.json",
        "synthetic",
        "experimental",
        "Researcher report",
        "report.html",
        "evidence bundle",
    ]:
        assert marker in text
