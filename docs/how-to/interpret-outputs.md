# How to interpret motif and co-occurrence outputs

This guide explains what the pipeline output files mean and how to draw conclusions from them.

---

## motif_results.tsv

Located at `results/analysis/motif_results.tsv`.

| Column | Type | Description |
|--------|------|-------------|
| `seq_id` | string | Sequence identifier from the alignment FASTA |
| `taxid` | integer | NCBI taxonomy ID of the organism |
| `motif_present` | bool | `True` if the motif was found at the configured positions |
| `matched_fragment` | string | The actual residues found at the motif positions (or empty) |
| `position` | integer | The (ungapped) alignment position where the match was found |

**Reading the results:**

- A row with `motif_present = True` means the sequence carries the exact residue pattern defined by `motif.fragments` at the configured positions (±`motif.tolerance`).
- If many sequences from the same organism have `motif_present = True`, that organism is considered **motif-positive**.

---

## cooccurrence.tsv

Located at `results/analysis/cooccurrence.tsv`.

| Key | Description |
|-----|-------------|
| `n_protein_taxa` | Number of unique taxids with a protein homolog |
| `n_rna_taxa` | Number of unique taxids with a 16S rRNA hit |
| `n_intersection` | Taxids present in **both** sets |
| `n_union` | Taxids present in **either** set |
| `jaccard_index` | `n_intersection / n_union` (0–1; higher = more overlap) |

**Reading the results:**

- A **Jaccard index close to 1** means nearly every organism with the protein also has the RNA feature, and vice versa.
- A **Jaccard index close to 0** means the two features appear in mostly different organisms.
- The raw counts (`n_intersection`, `n_union`) tell you the scale — a Jaccard of 0.8 from 5 organisms is far less informative than one from 500.

---

## fisher_stats.tsv / console output

The Fisher's exact test result is logged to the console and (if produced) written to `results/analysis/fisher_stats.tsv`.

| Key | Description |
|-----|-------------|
| `odds_ratio` | Odds ratio from the 2×2 contingency table |
| `p_value` | Two-sided p-value |

**Reading the results:**

- **p-value < 0.05** (or your chosen threshold): the co-occurrence is unlikely to be due to chance given the observed taxid sets.
- **odds_ratio > 1**: the protein and RNA feature are positively associated.
- **odds_ratio < 1**: they are negatively associated (anti-correlated).

!!! warning "Interpret with caution"
    A significant p-value does not establish causation. See [Fisher's exact test — caveats](../explanation/fisher-test.md) for limitations specific to this pipeline.

---

## phylum_summary.tsv

Located at `results/analysis/phylum_summary.tsv`.

| Column | Description |
|--------|-------------|
| `taxid` | NCBI taxonomy ID |
| `phylum` | Phylum name (currently `"Unknown"` — taxonomy integration is a planned feature) |

!!! note "Placeholder behaviour"
    Phylum assignment currently returns `"Unknown"` for all taxids because automated NCBI taxonomy lookup has not yet been integrated. This is a known limitation — see [Caveats & assumptions](../explanation/caveats.md).
