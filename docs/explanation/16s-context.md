# 16S rRNA and phylogenetics — background

This page provides a brief, accessible introduction to the biological context of this pipeline.

---

## What is 16S rRNA?

Ribosomes are the molecular machines that synthesise proteins in all living cells.  They are built from RNA and protein components.  In bacteria, the **small ribosomal subunit** contains a ~1,500-nucleotide RNA called **16S rRNA** (the "16S" refers to its sedimentation coefficient).

16S rRNA is useful for evolutionary and ecological studies for two reasons:

1. **It is universally present** — every bacterium has it.
2. **It evolves slowly** — its sequence is conserved enough to be comparable across distant organisms, yet variable enough to distinguish species.

For decades, sequencing 16S rRNA has been the gold standard for **microbial taxonomy** and **phylogenetic reconstruction** (the "tree of life" for bacteria).

---

## Structural positions and motifs

The 16S rRNA sequence folds into a complex secondary structure of loops, stems, and bulges.  **Specific nucleotide positions** within this structure are critical for ribosome function:

- Some positions interact directly with transfer RNA (tRNA) during translation.
- Others contact the large ribosomal subunit (23S rRNA) or ribosome-associated proteins.
- Mutations at key positions can confer antibiotic resistance (e.g., the **A-site** around positions 1408–1410 is targeted by aminoglycoside antibiotics).

A **motif** in this context is a pattern of residues at defined alignment positions.  Because the 16S gene is **highly conserved**, the same position numbers (in the *E. coli* coordinate system) are broadly comparable across bacteria.

---

## Why co-occurrence analysis?

If a protein family and a specific 16S rRNA structural motif are found together in the same organisms significantly more (or less) often than expected by chance, it suggests a **functional or evolutionary link**:

- The protein might interact with that region of the ribosome.
- Both features might be under the same selective pressure.
- The protein might modify the RNA at those positions.

This pipeline provides a **first-pass** statistical test of co-occurrence.  A significant result motivates deeper investigation (e.g., structural modelling, experimental validation) but does not by itself establish a mechanism.

---

## Multiple sequence alignment

Before detecting a motif across many organisms, all 16S sequences must be **aligned** — gaps are inserted so that homologous positions are in the same column.  This pipeline uses **MAFFT**, a fast and accurate multiple-sequence aligner, to produce a gapped alignment.

Motif positions are then interpreted relative to the **ungapped** (reference) coordinate system, so that a position number always refers to the same location regardless of how many gaps the aligner inserted.
