# Config reference

All pipeline behaviour is controlled by `config.yaml`.  This page describes every key.

---

## Full annotated example

```yaml
blast:
  protein_db: nr            # BLAST protein database name or path
  nucleotide_db: core_nt    # BLAST nucleotide database name or path
  threads: 8                # CPU threads passed to blastp/blastn via -num_threads
  evalue: 0.05              # E-value cutoff passed to blastp/blastn via -evalue
  max_target_seqs: 50000    # Maximum hits per query (-max_target_seqs)

filters:
  min_identity: 30          # Minimum % identity to retain a BLAST hit (0–100)
  min_alignment_coverage: 0.7   # Minimum alignment coverage fraction (0–1)

motif:
  # One-based residue positions in the ungapped 16S rRNA sequence.
  # Each fragment is a pattern string describing consecutive residues.
  # Syntax: plain char = exact; [ag] = a or g; {a} = anything except a; x = any.
  # positions:[100,101,102] fragments:["a","g","c"] is identical to
  # positions:[100] fragments:["agc"].
  positions: [1408, 1409, 1410]
  fragments: ["A", "C", "G"]
  molecule_type: rna   # protein | dna | rna — controls alphabet normalisation
  tolerance: 0         # allow ±tolerance shift around each configured position

alignment:
  method: mafft          # alignment tool (currently only mafft is supported)

output:
  results_dir: results   # root directory for all pipeline outputs

conda_env: environment.yml   # conda environment file used by Snakemake --use-conda
```

---

## Section: `blast`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `protein_db` | string | `nr` | BLAST protein database (name or absolute path) |
| `nucleotide_db` | string | `core_nt` | BLAST nucleotide database (name or absolute path) |
| `threads` | integer | `8` | Number of CPU threads for BLAST (`-num_threads`) |
| `evalue` | float | `0.05` | E-value threshold for BLAST hits |
| `max_target_seqs` | integer | `50000` | Maximum number of aligned sequences to keep per query |

---

## Section: `filters`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `min_identity` | float | `30` | Minimum percentage identity (0–100) to retain a hit |
| `min_alignment_coverage` | float | `0.7` | Minimum fraction of query covered by the alignment (0–1) |

---

## Section: `motif`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `positions` | list of int | `[1408, 1409, 1410]` | One-based ungapped positions to interrogate |
| `fragments` | list of string | `["A", "C", "G"]` | Pattern strings at each position (see pattern syntax below) |
| `molecule_type` | string | `rna` | Sequence alphabet: `rna`, `dna`, or `protein` |
| `tolerance` | integer | `0` | Allow ±N shift around each position |

### Pattern syntax

| Pattern | Matches |
|---------|---------|
| `A` | Exactly `A` (case-insensitive) |
| `[AG]` | `A` or `G` |
| `{A}` | Any residue except `A` |
| `x` | Any single residue |
| `AGC` | Three consecutive residues `A`, `G`, `C` |

### Legacy key (`residues`)

If `motif.residues` is present in the config it is automatically normalised to `motif.fragments`.  New configs should always use `fragments`.

---

## Section: `alignment`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `method` | string | `mafft` | Multiple alignment tool.  Only `mafft` is currently supported. |

---

## Section: `output`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `results_dir` | string | `results` | Root directory for all output files |

---

## `conda_env`

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `conda_env` | string | `environment.yml` | Conda environment file used when Snakemake is invoked with `--use-conda` |
