# coevo_rnase_16s

**A modular, reproducible pipeline to investigate co-occurrence of a protein family and a structural motif in bacterial 16S rRNA sequences.**

---

## What is this pipeline?

`coevo_rnase_16s` searches for organisms that contain **both**:

1. Homologs of a given protein family (detected via BLASTP against a local `nr` database).
2. A specific structural motif in their 16S rRNA (defined by residue patterns at alignment positions, detected via BLASTN + MAFFT alignment).

The pipeline then quantifies the **statistical co-occurrence** of these two features across organisms using Jaccard similarity and Fisher's exact test, and groups results by phylum.

This is useful for exploring whether a particular protein family and an RNA structural feature are found together more often than expected by chance — an initial step towards studying potential **co-evolutionary signals**.

---

## Quick Start

### Prerequisites

- Python ≥ 3.10
- BLAST+ suite (`blastp`, `blastn`) in `$PATH`
- MAFFT in `$PATH`
- Local BLAST databases: `nr` (protein), `core_nt` (nucleotide)

### Install

```bash
git clone https://github.com/biofe/coevo_pipeline.git
cd coevo_pipeline
pip install -e .
```

### Run the pipeline on example data

=== "CLI"

    ```bash
    coevo run-all examples/protein_query.fasta examples/small_16s.fasta \
        --config config.yaml
    ```

=== "Snakemake"

    ```bash
    snakemake --cores 8 --configfile config.yaml
    ```

Results will be written under `results/`.

---

## Pipeline steps at a glance

| Step | Description | Output |
|------|-------------|--------|
| 1 — BLASTP | Find protein family homologs in `nr` | `results/blast/protein_hits.tsv` |
| 2 — BLASTN | Find 16S rRNA sequences in `core_nt` | `results/blast/rna_hits.tsv` |
| 3 — Extract taxa | Pull taxids from BLAST hits | `results/taxa/protein_taxa.txt`, `results/taxa/rna_taxa.txt` |
| 4 — Prepare alignment | Build reference + hit FASTA | `results/alignment/16s_input.fasta`, `results/alignment/16s_metadata.tsv` |
| 5 — MAFFT alignment | Align 16S sequences | `results/alignment/16s_alignment.fasta` |
| 6 — Motif detection | Score each sequence for the motif | `results/analysis/motif_results.tsv` |
| 7–8 — Co-occurrence & stats | Jaccard, Fisher's test, phylum summary | `results/analysis/cooccurrence.tsv`, `results/analysis/phylum_summary.tsv` |

---

## Documentation structure

This site follows the [Diátaxis](https://diataxis.fr/) framework:

- **[Tutorials](tutorials/end-to-end.md)** — learning-oriented: run the pipeline end-to-end.
- **[How-to Guides](how-to/blast-only.md)** — task-oriented: accomplish specific goals.
- **[Reference](reference/cli.md)** — information-oriented: CLI commands, config schema, output files.
- **[Explanation](explanation/16s-context.md)** — understanding-oriented: background concepts, caveats.

---

## Repository

Source code: [github.com/biofe/coevo_pipeline](https://github.com/biofe/coevo_pipeline)

License: MIT
