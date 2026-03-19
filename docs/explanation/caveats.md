# Caveats and assumptions

This page documents known limitations and assumptions of the pipeline.  Understanding them is important for correctly interpreting results.

---

## Local BLAST databases

The pipeline runs BLAST against **local copies** of NCBI databases (`nr`, `core_nt`).  This means:

- Results depend on the **version and completeness** of the local database.
- Databases go out of date.  Run BLAST against an old `nr` and you will miss recently deposited sequences.
- Database paths must be configured manually in `config.yaml` and are not validated before the pipeline starts.

**Mitigation:** document the database version and download date in your `config.yaml` or lab notebook. Use `blastdbcmd -db <path> -info` to check the database build date.

---

## Taxid assignment

Taxids are extracted from the `staxids` field of BLAST tabular output.  This field:

- Is populated only if the database was built with taxonomy information (e.g., NCBI `nr` and `core_nt` include it; custom databases may not).
- Can be **semicolon-separated** (one hit may map to multiple taxids).  The pipeline splits on `;` and takes all taxids.
- Does not distinguish **strain-level** from **species-level** taxids.  Two strains of the same species will appear as different taxids if NCBI assigned them separately.

---

## "Phylum summary" is a placeholder

The `phylum_summary.tsv` output currently assigns `"Unknown"` to all taxids because automated NCBI taxonomy lookup (via Entrez or a local taxonomy database) has not yet been integrated into the pipeline.

**Impact:** the phylum summary table does not currently provide meaningful groupings.

**Future work:** integrate NCBI `ete3` or `taxonkit` to map taxids to phylum-level classifications.

---

## Motif positions use the *E. coli* coordinate system

Motif positions are specified in **one-based, ungapped coordinates** following the *E. coli* K-12 16S rRNA numbering convention.  If your reference sequence uses a different numbering system, you must convert positions before configuring the pipeline.

---

## MAFFT alignment quality

Multiple sequence alignment is inherently uncertain for divergent sequences.  If the 16S sequences in your dataset are very diverse (e.g., spanning phyla with >30% sequence divergence), the alignment quality may be reduced and motif position detection may be unreliable.

Consider:

- Inspecting the alignment manually before proceeding.
- Using a curated reference alignment (e.g., from SILVA) as a structural scaffold.

---

## Presence/absence only

The current version of the pipeline measures **presence or absence** of the protein family and the motif.  It does not:

- Quantify copy number or expression level.
- Distinguish between active and inactive/pseudogene copies.
- Model the phylogenetic non-independence of organisms (related organisms are more likely to share features simply because of shared ancestry, not because of co-evolution).

**Future work:** phylogenetic comparative methods (e.g., ancestral state reconstruction) would be needed to properly test for co-evolution beyond simple co-occurrence.

---

## No multiple-testing correction

If the pipeline is used to screen many protein families or motif definitions, the Fisher's exact test p-values are **not automatically corrected** for multiple comparisons.  Apply Bonferroni, Benjamini–Hochberg, or another appropriate correction manually.
