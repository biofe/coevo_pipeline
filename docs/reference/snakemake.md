# Snakemake rules reference

The Snakemake workflow is defined in `Snakefile` (which includes three rule modules) and produces three final targets.

---

## Entry point: `Snakefile`

```python
configfile: "config.yaml"

include: "workflows/rules_blast.smk"
include: "workflows/rules_alignment.smk"
include: "workflows/rules_analysis.smk"

rule all:
    input:
        expand("{results}/analysis/cooccurrence.tsv",   results=RESULTS),
        expand("{results}/analysis/phylum_summary.tsv", results=RESULTS),
        expand("{results}/analysis/motif_results.tsv",  results=RESULTS),
```

---

## Rule dependency graph (high-level)

```
all
├── analyse           (rules_analysis.smk)
│   ├── extract_taxa  (rules_blast.smk)
│   │   ├── blast_protein  (rules_blast.smk)
│   │   └── blast_16s      (rules_blast.smk)
│   └── detect_motif  (rules_alignment.smk)
│       └── align_16s      (rules_alignment.smk)
│           └── prepare_alignment  (rules_alignment.smk)
│               └── blast_16s      (rules_blast.smk)
└── (motif_results.tsv already produced by detect_motif above)
```

---

## `workflows/rules_blast.smk`

### `rule blast_protein`

| | |
|---|---|
| **Input** | `config["blast"]["protein_query"]` (protein query FASTA) |
| **Output** | `{RESULTS}/blast/protein_hits.tsv` |
| **Calls** | `coevo blast-protein {input} --config {CONFIG_FILE}` |

Runs BLASTP against `blast.protein_db` with the configured e-value, threads and max_target_seqs.

### `rule blast_16s`

| | |
|---|---|
| **Input** | `config["blast"]["rna_query"]` (16S rRNA query FASTA) |
| **Output** | `{RESULTS}/blast/rna_hits.tsv` |
| **Calls** | `coevo blast-16s {input} --config {CONFIG_FILE}` |

Runs BLASTN against `blast.nucleotide_db`.

### `rule extract_taxa`

| | |
|---|---|
| **Input** | `{RESULTS}/blast/protein_hits.tsv`, `{RESULTS}/blast/rna_hits.tsv` |
| **Output** | `{RESULTS}/taxa/protein_taxa.txt`, `{RESULTS}/taxa/rna_taxa.txt` |
| **Calls** | `coevo extract-taxa --config {CONFIG_FILE}` |

---

## `workflows/rules_alignment.smk`

### `rule prepare_alignment`

| | |
|---|---|
| **Input** | `config["blast"]["rna_query"]`, `{RESULTS}/blast/rna_hits.tsv` |
| **Output** | `{RESULTS}/alignment/16s_input.fasta`, `{RESULTS}/alignment/16s_metadata.tsv` |
| **Calls** | `coevo prepare-alignment {input.query} --config {CONFIG_FILE}` |

### `rule align_16s`

| | |
|---|---|
| **Input** | `{RESULTS}/alignment/16s_input.fasta` |
| **Output** | `{RESULTS}/alignment/16s_alignment.fasta` |
| **Calls** | `coevo align-16s {input} --config {CONFIG_FILE}` |

### `rule detect_motif`

| | |
|---|---|
| **Input** | `{RESULTS}/alignment/16s_alignment.fasta` |
| **Output** | `{RESULTS}/analysis/motif_results.tsv` |
| **Calls** | `coevo detect-motif --config {CONFIG_FILE}` |

---

## `workflows/rules_analysis.smk`

### `rule analyse`

| | |
|---|---|
| **Input** | `{RESULTS}/taxa/protein_taxa.txt`, `{RESULTS}/taxa/rna_taxa.txt`, `{RESULTS}/analysis/motif_results.tsv` |
| **Output** | `{RESULTS}/analysis/cooccurrence.tsv`, `{RESULTS}/analysis/phylum_summary.tsv` |
| **Calls** | `coevo analyse --config {CONFIG_FILE}` |

---

## Useful Snakemake invocations

```bash
# Preview the workflow without executing
snakemake --dry-run --cores 8 --configfile config.yaml

# Execute the full pipeline
snakemake --cores 8 --configfile config.yaml

# Force re-run a specific rule
snakemake --cores 8 --configfile config.yaml --forcerun detect_motif

# Print rule execution order
snakemake --dag --configfile config.yaml | dot -Tsvg > dag.svg

# Use conda environments
snakemake --cores 8 --configfile config.yaml --use-conda
```
