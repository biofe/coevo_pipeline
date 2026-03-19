# CLI reference

The `coevo` command provides individual pipeline steps and an end-to-end runner.

All commands accept these global options:

| Option | Default | Description |
|--------|---------|-------------|
| `--config` / `-c` | `config.yaml` | Path to the YAML configuration file |
| `--log-level` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |
| `--help` | — | Show help and exit |

---

## `coevo blast-protein`

Run BLASTP to identify protein family homologs.

```
coevo blast-protein QUERY [OPTIONS]
```

| Argument / Option | Required | Description |
|-------------------|----------|-------------|
| `QUERY` | yes | Path to protein query FASTA |
| `--config` | no | Path to config.yaml (default: `config.yaml`) |
| `--log-level` | no | Logging level (default: `INFO`) |

**Output:** `{results_dir}/blast/protein_hits.tsv`

---

## `coevo blast-16s`

Run BLASTN to find 16S rRNA sequences for relevant organisms.

```
coevo blast-16s QUERY [OPTIONS]
```

| Argument / Option | Required | Description |
|-------------------|----------|-------------|
| `QUERY` | yes | Path to 16S rRNA query FASTA |
| `--config` | no | Path to config.yaml (default: `config.yaml`) |
| `--log-level` | no | Logging level (default: `INFO`) |

**Output:** `{results_dir}/blast/rna_hits.tsv`

---

## `coevo extract-taxa`

Extract taxids from BLAST results and write taxid lists.

```
coevo extract-taxa [OPTIONS]
```

Reads `{results_dir}/blast/protein_hits.tsv` and `{results_dir}/blast/rna_hits.tsv`.

**Outputs:**

- `{results_dir}/taxa/protein_taxa.txt`
- `{results_dir}/taxa/rna_taxa.txt`

---

## `coevo prepare-alignment`

Build alignment input FASTA and sequence metadata from BLAST results.

```
coevo prepare-alignment QUERY [OPTIONS]
```

| Argument / Option | Required | Description |
|-------------------|----------|-------------|
| `QUERY` | yes | Path to reference 16S rRNA FASTA (same file used for `blast-16s`) |
| `--config` | no | Path to config.yaml (default: `config.yaml`) |

Reads `{results_dir}/blast/rna_hits.tsv`, filters by `filters.min_identity`, deduplicates hit sequences, and prepends the reference sequence.

**Outputs:**

- `{results_dir}/alignment/16s_input.fasta`
- `{results_dir}/alignment/16s_metadata.tsv`

---

## `coevo align-16s`

Align 16S rRNA sequences using MAFFT.

```
coevo align-16s INPUT_FASTA [OPTIONS]
```

| Argument / Option | Required | Description |
|-------------------|----------|-------------|
| `INPUT_FASTA` | yes | Path to 16S FASTA to align (usually `results/alignment/16s_input.fasta`) |
| `--config` | no | Path to config.yaml (default: `config.yaml`) |

**Output:** `{results_dir}/alignment/16s_alignment.fasta`

---

## `coevo detect-motif`

Detect the 16S rRNA structural motif in the alignment.

```
coevo detect-motif [OPTIONS]
```

Uses `motif.positions`, `motif.fragments`, `motif.molecule_type`, and `motif.tolerance` from config.

**Output:** `{results_dir}/analysis/motif_results.tsv`

---

## `coevo analyse`

Compute co-occurrence statistics and produce summary tables.

```
coevo analyse [OPTIONS]
```

Reads `protein_taxa.txt` and `rna_taxa.txt`, computes Jaccard co-occurrence and Fisher's exact test, and writes a phylum summary.

**Outputs:**

- `{results_dir}/analysis/cooccurrence.tsv`
- `{results_dir}/analysis/phylum_summary.tsv`

---

## `coevo run-all`

Run the full pipeline end-to-end.

```
coevo run-all PROTEIN_QUERY RNA_QUERY [OPTIONS]
```

| Argument | Required | Description |
|----------|----------|-------------|
| `PROTEIN_QUERY` | yes | Path to protein query FASTA |
| `RNA_QUERY` | yes | Path to 16S rRNA query FASTA |
| `--config` | no | Path to config.yaml (default: `config.yaml`) |

Calls all individual steps in order.
