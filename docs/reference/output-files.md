# Output files reference

All output files are written under `{output.results_dir}` (default: `results/`).

---

## Directory layout

```
results/
├── blast/
│   ├── protein_hits.tsv        # BLASTP tabular output
│   └── rna_hits.tsv            # BLASTN tabular output
├── taxa/
│   ├── protein_taxa.txt        # One taxid per line (protein homolog organisms)
│   └── rna_taxa.txt            # One taxid per line (16S rRNA hit organisms)
├── alignment/
│   ├── 16s_input.fasta         # Reference + unique 16S hit sequences (alignment input)
│   ├── 16s_metadata.tsv        # seq_index → taxid mapping
│   └── 16s_alignment.fasta     # MAFFT-aligned 16S sequences
└── analysis/
    ├── motif_results.tsv       # Per-sequence motif presence/absence
    ├── cooccurrence.tsv        # Jaccard index and set sizes
    ├── fisher_stats.tsv        # Odds ratio and p-value
    └── phylum_summary.tsv      # Taxids grouped by phylum
```

---

## `results/blast/protein_hits.tsv`

**Produced by:** `coevo blast-protein` / `rule blast_protein`

BLAST tabular output (outfmt 6) with an additional `staxids` column.

| Column | Description |
|--------|-------------|
| `qseqid` | Query sequence ID |
| `sseqid` | Subject (database) sequence ID |
| `pident` | Percentage identity |
| `length` | Alignment length |
| `mismatch` | Number of mismatches |
| `gapopen` | Number of gap openings |
| `qstart` | Query start position |
| `qend` | Query end position |
| `sstart` | Subject start position |
| `send` | Subject end position |
| `evalue` | Expect value |
| `bitscore` | Bit score |
| `staxids` | Semicolon-separated NCBI taxids of the subject |

---

## `results/blast/rna_hits.tsv`

**Produced by:** `coevo blast-16s` / `rule blast_16s`

Same column schema as `protein_hits.tsv`, but from a BLASTN search against the nucleotide database.

---

## `results/taxa/protein_taxa.txt`

**Produced by:** `coevo extract-taxa` / `rule extract_taxa`

Plain text file, one NCBI taxid per line, representing organisms with a protein family homolog.

---

## `results/taxa/rna_taxa.txt`

**Produced by:** `coevo extract-taxa` / `rule extract_taxa`

Plain text file, one NCBI taxid per line, representing organisms with a 16S rRNA database hit.

---

## `results/alignment/16s_input.fasta`

**Produced by:** `coevo prepare-alignment` / `rule prepare_alignment`

FASTA file containing:

1. The reference 16S rRNA sequence (from the query FASTA passed to `blast-16s`).
2. All unique subject sequences from `rna_hits.tsv` that pass the identity filter.

Sequence headers follow the enumeration-based naming convention
``seq_{index}``, where the index matches the ``seq_index`` column in
``16s_metadata.tsv``:

- ``seq_0`` — the reference sequence (index 0).
- ``seq_1``, ``seq_2``, … — unique BLAST hit sequences in encounter order.

This naming is preserved through MAFFT alignment into ``16s_alignment.fasta``
and is used as the ``sequence_id`` in ``motif_results.tsv``, enabling
unambiguous cross-referencing across all pipeline output files.

Used as input to MAFFT.

---

## `results/alignment/16s_metadata.tsv`

**Produced by:** `coevo prepare-alignment` / `rule prepare_alignment`

| Column | Description |
|--------|-------------|
| `seq_index` | Local integer index (0-based) in `16s_input.fasta` |
| `taxid` | NCBI taxid of the sequence |

Maps each sequence in the alignment FASTA to its organism.  The ``seq_index``
value corresponds directly to the numeric part of the sequence identifier in
``16s_input.fasta`` and ``16s_alignment.fasta`` (e.g. ``seq_index = 1``
→ header ``>seq_1``).

---

## `results/alignment/16s_alignment.fasta`

**Produced by:** `coevo align-16s` / `rule align_16s`

MAFFT multiple-sequence alignment of all sequences in `16s_input.fasta`.  Gap characters (`-`) are inserted to produce a common coordinate system.  Sequence headers (``seq_0``, ``seq_1``, …) are inherited unchanged from ``16s_input.fasta`` and can be looked up in ``16s_metadata.tsv`` via the ``seq_index`` column.

---

## `results/analysis/motif_results.tsv`

**Produced by:** `coevo detect-motif` / `rule detect_motif`

| Column | Type | Description |
|--------|------|-------------|
| `sequence_id` | string | Enumeration-based sequence identifier from the alignment FASTA (e.g. `seq_0`, `seq_1`) |
| `motif_present` | bool | `True` if the motif was found |

The ``sequence_id`` values match the FASTA headers in ``16s_alignment.fasta``
and can be joined to ``16s_metadata.tsv`` via the numeric part of the identifier
(``seq_index``).

---

## `results/analysis/cooccurrence.tsv`

**Produced by:** `coevo analyse` / `rule analyse`

| Key | Description |
|-----|-------------|
| `n_protein_taxa` | Number of unique taxids with a protein homolog |
| `n_rna_taxa` | Number of unique taxids with a 16S rRNA hit |
| `n_intersection` | Taxids in **both** sets |
| `n_union` | Taxids in **either** set |
| `jaccard_index` | `n_intersection / n_union` (0–1) |

---

## `results/analysis/fisher_stats.tsv`

**Produced by:** `coevo analyse` / `rule analyse`

| Key | Description |
|-----|-------------|
| `odds_ratio` | Odds ratio from the 2×2 contingency table |
| `p_value` | Two-sided Fisher's exact test p-value |

---

## `results/analysis/phylum_summary.tsv`

**Produced by:** `coevo analyse` / `rule analyse`

| Column | Description |
|--------|-------------|
| `taxid` | NCBI taxid |
| `phylum` | Phylum name (currently `"Unknown"` — see [caveats](../explanation/caveats.md)) |
