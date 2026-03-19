# Co-occurrence and the Jaccard index

This page explains how the pipeline measures whether a protein family and an RNA structural motif tend to occur in the same organisms.

---

## The two organism sets

The pipeline constructs two sets of NCBI taxids:

- **Protein set (P):** organisms with at least one BLASTP hit above the configured thresholds — i.e., organisms that appear to carry a homolog of the query protein.
- **RNA set (R):** organisms with at least one BLASTN hit whose 16S rRNA carries the defined structural motif.

Each organism is identified by its NCBI **taxid** (a unique integer identifier).

---

## Intersection and union

| Concept | Meaning |
|---------|---------|
| **Intersection** \( P \cap R \) | Taxids present in **both** sets — organisms with both the protein and the motif |
| **Union** \( P \cup R \) | Taxids present in **either** set — organisms with at least one feature |

---

## Jaccard similarity index

The **Jaccard index** is a simple, normalised measure of overlap between two sets:

$$
J(P, R) = \frac{|P \cap R|}{|P \cup R|}
$$

- **J = 1**: perfect overlap — every organism with one feature also has the other.
- **J = 0**: no overlap — the two features are found in completely different organisms.
- **J ≈ 0.5**: moderate overlap.

The Jaccard index is symmetric and intuitive, but ignores the **total number of organisms** in the database.  A Jaccard of 0.9 from 10 organisms is much weaker evidence than 0.9 from 10,000 organisms — this is why the Fisher's exact test is also computed.

---

## Worked example

Suppose:

- Protein set: taxids `{1, 2, 3, 4, 5}` → 5 organisms
- RNA motif set: taxids `{3, 4, 5, 6, 7}` → 5 organisms
- Intersection: `{3, 4, 5}` → 3 organisms
- Union: `{1, 2, 3, 4, 5, 6, 7}` → 7 organisms

$$J = \frac{3}{7} \approx 0.43$$

---

## What is stored in `cooccurrence.tsv`

| Key | Value (example) |
|-----|-----------------|
| `n_protein_taxa` | 5 |
| `n_rna_taxa` | 5 |
| `n_intersection` | 3 |
| `n_union` | 7 |
| `jaccard_index` | 0.4286 |

---

## Limitations of Jaccard alone

- It does not account for the background rate of either feature in the full set of sequenced organisms.
- It treats each taxid equally, regardless of how many sequences were found for it.
- It is not a statistical test — it cannot provide a p-value.

Use the [Fisher's exact test](fisher-test.md) output alongside the Jaccard index for a more complete picture.
