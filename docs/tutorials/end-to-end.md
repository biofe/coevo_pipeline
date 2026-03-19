# Tutorial: Run the pipeline end-to-end

This tutorial walks you through a complete run of `coevo_rnase_16s` on the provided example data.

**Goal:** understand every step, what it produces, and where to find the results.

---

## Prerequisites

Before starting, make sure you have:

- `coevo` installed (`pip install -e .` from the repo root)
- BLAST+ suite in `$PATH` (`blastp`, `blastn`)
- MAFFT in `$PATH`
- Local BLAST databases available: `nr` (protein), `core_nt` (nucleotide)

Check your installation:

```bash
coevo --help
blastp -version
blastn -version
mafft --version
```

---

## Step 0 — Review the configuration

Open `config.yaml` and verify the database paths and motif definition match your environment:

```yaml
blast:
  protein_db: nr          # path or name of your local nr database
  nucleotide_db: core_nt  # path or name of your local core_nt database
  threads: 8
  evalue: 0.05
  max_target_seqs: 50000

filters:
  min_identity: 30
  min_alignment_coverage: 0.7   # optional; fraction of query covered (0–1)

motif:
  positions: [1408, 1409, 1410]
  fragments: ["A", "C", "G"]
  molecule_type: rna
  tolerance: 0

output:
  results_dir: results
```

!!! note "Example motif"
    Positions 1408–1410 (in the *E. coli* K-12 coordinate system) correspond to the **A-site** of the 16S rRNA — a region critical for decoding accuracy and the binding site of aminoglycoside antibiotics.  The residues `A`, `C`, `G` represent the wild-type nucleotides at these positions in most bacteria.  Change these values to define your own motif of interest.

!!! tip
    For a local development run, you can point `protein_db` and `nucleotide_db` at small test databases to avoid waiting for full NCBI searches.

---

## Step 1 — Run BLASTP (protein homolog search)

```bash
coevo blast-protein examples/protein_query.fasta --config config.yaml
```

This runs `blastp` against the configured protein database and writes tabular output (outfmt 6 + staxids) to:

```
results/blast/protein_hits.tsv
```

**What it contains:** one line per HSP (high-scoring pair), with columns `qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, `bitscore`, `staxids`.

---

## Step 2 — Run BLASTN (16S rRNA search)

```bash
coevo blast-16s examples/small_16s.fasta --config config.yaml
```

This runs `blastn` against the nucleotide database and writes:

```
results/blast/rna_hits.tsv
```

Same column schema as the protein hits table.

---

## Step 3 — Extract taxids

```bash
coevo extract-taxa --config config.yaml
```

Parses both BLAST hit tables, extracts the `staxids` column (semicolon-separated), and writes one taxid per line to:

```
results/taxa/protein_taxa.txt
results/taxa/rna_taxa.txt
```

---

## Step 4 — Prepare alignment input

```bash
coevo prepare-alignment examples/small_16s.fasta --config config.yaml
```

Reads the 16S BLAST hits, filters by minimum identity (`filters.min_identity`), deduplicates subject sequences, and assembles:

- `results/alignment/16s_input.fasta` — the reference sequence (from `small_16s.fasta`) followed by all unique BLAST hits.
- `results/alignment/16s_metadata.tsv` — mapping of local sequence index → taxid.

---

## Step 5 — Align sequences

```bash
coevo align-16s results/alignment/16s_input.fasta --config config.yaml
```

Calls MAFFT on the input FASTA and writes the multiple-sequence alignment to:

```
results/alignment/16s_alignment.fasta
```

---

## Step 6 — Detect the structural motif

```bash
coevo detect-motif --config config.yaml
```

Scans every aligned sequence for the motif defined in `config.yaml` (`motif.positions` + `motif.fragments`).

Output:

```
results/analysis/motif_results.tsv
```

Columns: `seq_id`, `taxid`, `motif_present` (True/False), `matched_fragment`, `position`.

---

## Step 7–8 — Co-occurrence & summary statistics

```bash
coevo analyse --config config.yaml
```

Computes:

- **Co-occurrence statistics** (Jaccard index, intersection, union sizes) → `results/analysis/cooccurrence.tsv`
- **Fisher's exact test** (contingency table odds ratio and p-value; logged to console)
- **Phylum summary** → `results/analysis/phylum_summary.tsv`

---

## Run everything at once

### CLI

```bash
coevo run-all examples/protein_query.fasta examples/small_16s.fasta \
    --config config.yaml
```

### Snakemake

```bash
# Preview without executing
snakemake --dry-run --cores 8 --configfile config.yaml

# Execute
snakemake --cores 8 --configfile config.yaml
```

---

## Expected outputs

After a successful run, `results/` will look like this:

```
results/
├── blast/
│   ├── protein_hits.tsv
│   └── rna_hits.tsv
├── taxa/
│   ├── protein_taxa.txt
│   └── rna_taxa.txt
├── alignment/
│   ├── 16s_input.fasta
│   ├── 16s_metadata.tsv
│   └── 16s_alignment.fasta
└── analysis/
    ├── motif_results.tsv
    ├── cooccurrence.tsv
    ├── fisher_stats.tsv
    └── phylum_summary.tsv
```

See [Output files reference](../reference/output-files.md) for a description of every column in each file.

---

## Next steps

- [How to run only the BLAST steps](../how-to/blast-only.md) — skip downstream analysis
- [How to change the motif definition](../how-to/motif-config.md)
- [Interpret your results](../how-to/interpret-outputs.md)
