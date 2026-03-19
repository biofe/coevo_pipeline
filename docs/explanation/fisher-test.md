# Fisher's exact test

This page explains the statistical test used to assess whether the co-occurrence of the protein and RNA motif is significant.

---

## The contingency table

Fisher's exact test operates on a **2×2 contingency table** that counts organisms in each combination of presence/absence for the two features:

|  | Motif present | Motif absent |
|--|---------------|--------------|
| **Protein present** | a | b |
| **Protein absent** | c | d |

Where:

- **a** = organisms with **both** the protein and the motif (`|P ∩ R|`)
- **b** = organisms with the protein **but not** the motif (`|P| - a`)
- **c** = organisms with the motif **but not** the protein (`|R| - a`)
- **d** = organisms with **neither** feature (estimated from the total database size, or treated as the complement)

---

## What the test computes

Fisher's exact test computes the exact probability of observing a table as extreme as the one seen, under the null hypothesis that the two features are **independent** (i.e., knowing an organism has the protein tells you nothing about whether it has the motif).

The pipeline reports:

| Statistic | Meaning |
|-----------|---------|
| `odds_ratio` | How much more likely it is to find the motif in protein-positive organisms vs. protein-negative organisms |
| `p_value` | Two-sided p-value; small values (<0.05) suggest non-random co-occurrence |

---

## Interpreting the odds ratio

- **odds_ratio > 1**: protein-positive organisms are more likely to carry the motif (positive association).
- **odds_ratio = 1**: no association.
- **odds_ratio < 1**: protein-positive organisms are less likely to carry the motif (negative association).

---

## Caveats and limitations

### The "absent" count is uncertain

The number of organisms with **neither** feature (cell **d** in the table) depends on the total number of organisms in the searched database.  Because the pipeline uses local BLAST databases of variable completeness, this count is an approximation.

### Taxid-level, not sequence-level

The test operates on **unique taxids**, not individual sequences.  One taxid with 100 16S sequences counts the same as one taxid with a single sequence.  This avoids pseudoreplication but means the power of the test depends on the number of distinct organisms, not the number of sequences.

### Multiple testing

If you run this pipeline for multiple protein families or motifs, the p-values should be corrected for multiple testing (e.g., Bonferroni or Benjamini–Hochberg correction).  The pipeline does not currently perform this correction automatically.

### The result is exploratory

A significant p-value means the co-occurrence pattern is unlikely to be due to sampling noise alone.  It does **not** establish:

- Causation or a direct functional link
- That the association holds across all environments or phylogenetic lineages
- That the statistical model is correctly specified

Significant results should be followed up with biological validation.
