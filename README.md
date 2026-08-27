# Geison

Geison is a pipeline for the **in silico design and evaluation of qPCR and RT-qPCR assays**.

The project is under active development. The current implementation lives on the [`develop`](https://github.com/BrunoDCamargo/Geison/tree/develop) branch; `main` serves as the stable project landing page until the development work is integrated.

## Current development scope

The pipeline currently covers the path from a quality-controlled sequence population to auditable primer and probe candidates:

```text
Evaluation Set
    ↓
optional clustering (CD-HIT)
    ↓
Discovery Set
    ↓
multiple-sequence alignment (MAFFT)
    ↓
genomic conservation analysis
    ↓
candidate-region selection
    ↓
primer/probe design (Primer3)
```

The development branch preserves intermediate artifacts and reports in formats such as FASTA, TSV, and JSON. It also generates a self-contained HTML report for conservation analysis.

## Design choices

- The **Evaluation Set** retains the full sequence population approved by quality control.
- The **Discovery Set** may use clustering to reduce redundancy during discovery steps without replacing the Evaluation Set.
- Alignment, conservation analysis, and primer design are configurable stages.
- Published coordinates are traceable to the selected reference sequence.
- Primer3 input and output can be retained for auditability when assay design runs.

## What is next

Two major evaluation stages remain open in the current roadmap:

- [Inclusivity and degeneracy evaluation against the full target diversity](https://github.com/BrunoDCamargo/Geison/issues/8)
- [Specificity evaluation against off-target datasets](https://github.com/BrunoDCamargo/Geison/issues/9)

The current Primer3 stage therefore produces **assay candidates**, not a final assay-risk decision.

## Technology

The current development package targets **Python 3.10+** and uses Biopython and PyYAML. Optional pipeline stages call external tools when enabled:

- `cd-hit-est` for clustering
- `mafft` for multiple-sequence alignment
- `primer3_core` for primer and probe design

See the [`develop` README](https://github.com/BrunoDCamargo/Geison/blob/develop/README.md) for the current configuration model, generated artifacts, and implementation details.
