# coevo_rnase_16s

A modular, reproducible Python pipeline to investigate the **co-occurrence** of a protein family and a specific structural motif in bacterial **16S rRNA sequences** across organisms.

---

## Project Goal

The pipeline detects organisms that contain **both**:

1. Homologs of a given protein family (searched via BLASTP against a local `nr` database)
2. A specific structural motif in their 16S rRNA (defined by three consecutive residues at alignment positions)

**Current version:** presence/absence correlation  
**Future work:** residue coevolution analysis

---

## Repository Structure

```
coevo_pipeline/
├── README.md
├── LICENSE
├── pyproject.toml
├── Makefile
├── config.yaml
├── Snakefile
│
├── coevo/                      # Main Python package
│   ├── __init__.py
│   ├── config.py               # Configuration loading
│   ├── cli.py                  # Typer CLI entry point
│   ├── logging_utils.py        # Loguru-based logging
│   │
│   ├── blast/
│   │   ├── blast_runner.py     # BLASTP / BLASTN via subprocess
│   │   └── blast_parser.py     # Parse tabular BLAST output
│   │
│   ├── sequences/
│   │   ├── fasta_utils.py      # Read / write FASTA (Biopython)
│   │   ├── alignment.py        # MAFFT multiple alignment
│   │   └── motif_detection.py  # Motif detection in alignments
│   │
│   ├── taxonomy/
│   │   └── taxid_utils.py      # Read / write / intersect taxid sets
│   │
│   ├── analysis/
│   │   ├── cooccurrence.py     # Jaccard / co-occurrence statistics
│   │   ├── statistics.py       # Fisher's exact test
│   │   └── phylogeny.py        # Phylum-level summary tables
│   │
│   └── io/
│       └── result_writer.py    # Write DataFrames and dicts to TSV
│
├── workflows/                  # Snakemake rule files
│   ├── rules_blast.smk
│   ├── rules_alignment.smk
│   └── rules_analysis.smk
│
├── tests/
│   ├── test_blast_parser.py
│   ├── test_motif_detection.py
│   └── test_cooccurrence.py
│
└── examples/
    ├── protein_query.fasta
    ├── small_16s.fasta
    └── example_config.yaml
```

---

## Requirements

- Python ≥ 3.10
- BLAST+ suite (`blastp`, `blastn`) available in `$PATH`
- MAFFT available in `$PATH`
- Local BLAST databases: `nr` (protein), `core_nt` (nucleotide)

---

## Installation

```bash
# Clone the repository
git clone https://github.com/biofe/coevo_pipeline.git
cd coevo_pipeline

# Create a virtual environment
python -m venv .venv
source .venv/bin/activate

# Install the package and its dependencies
pip install -e .

# For development (includes pytest)
pip install -e ".[dev]"
```

---

## Configuration

All parameters are controlled via `config.yaml`:

```yaml
blast:
  protein_db: nr
  nucleotide_db: core_nt
  threads: 8
  evalue: 0.05
  max_target_seqs: 50000

filters:
  min_identity: 30
  min_alignment_coverage: 0.7

motif:
  positions: [1408, 1409, 1410]
  residues: ["A", "C", "G"]

alignment:
  method: mafft

output:
  results_dir: results
```

See `examples/example_config.yaml` for a fully annotated example.

---

## Usage

### CLI

The `coevo` command exposes individual pipeline steps:

```bash
# Run BLASTP against the local nr database
coevo blast-protein examples/protein_query.fasta --config config.yaml

# Run BLASTN for 16S rRNA sequences
coevo blast-16s examples/small_16s.fasta --config config.yaml

# Extract taxids from BLAST results
coevo extract-taxa --config config.yaml

# Build alignment input FASTA and sequence metadata from BLAST results
# (reference sequence first, then all unique 16S hits; writes 16s_input.fasta and 16s_metadata.tsv)
coevo prepare-alignment examples/small_16s.fasta --config config.yaml

# Align 16S rRNA sequences with MAFFT
coevo align-16s results/alignment/16s_input.fasta --config config.yaml

# Detect the structural motif in the alignment
coevo detect-motif --config config.yaml

# Compute co-occurrence statistics
coevo analyse --config config.yaml

# Run everything end-to-end
coevo run-all examples/protein_query.fasta examples/small_16s.fasta
```

Use `coevo --help` or `coevo <command> --help` for full options.

### Snakemake

```bash
# Dry-run to preview the workflow
snakemake --dry-run --cores 8 --configfile config.yaml

# Execute the full pipeline
snakemake --cores 8 --configfile config.yaml

# Run with Slurm (example)
snakemake --cores 8 --configfile config.yaml \
    --executor slurm --default-resources slurm_partition=main
```

---

## Intermediate Results

The pipeline writes all intermediate files under `results/` (configurable):

```
results/
├── blast/
│   ├── protein_hits.tsv        # BLASTP tabular output
│   └── rna_hits.tsv            # BLASTN tabular output
├── taxa/
│   ├── protein_taxa.txt        # One taxid per line
│   └── rna_taxa.txt
├── alignment/
│   ├── 16s_input.fasta         # Reference + unique 16S sequences (alignment input)
│   ├── 16s_metadata.tsv        # seq_index → taxid mapping for each 16S sequence
│   └── 16s_alignment.fasta     # MAFFT-aligned 16S sequences
└── analysis/
    ├── motif_results.tsv       # Per-sequence motif presence/absence
    ├── cooccurrence.tsv        # Jaccard index and counts
    ├── fisher_stats.tsv        # Odds ratio and p-value
    └── phylum_summary.tsv      # Organisms grouped by phylum
```

---

## Running Tests

```bash
# Run all tests
make test

# With coverage
make test-cov

# Specific test module
pytest tests/test_blast_parser.py -v
```

---

## Future Extensions

- **Residue coevolution analysis** – co-variation statistics between alignment columns
- **HMM-based protein detection** – replace BLAST with `hmmsearch` for sensitive family detection
- **MMseqs2 integration** – faster sequence clustering and search
- **Phylogenetic tree colouring** – visualise co-occurrence on a species tree using ETE3
- **Multiple motif support** – scan for several motifs simultaneously
- **Entrez taxonomy integration** – automated phylum assignment via NCBI E-utilities

---

## License

MIT License – see [LICENSE](LICENSE).
