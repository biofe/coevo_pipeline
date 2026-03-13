"""Snakemake rules for co-occurrence analysis and summary statistics."""


rule cooccurrence_analysis:
    """Compute co-occurrence statistics from protein and RNA taxid sets."""
    input:
        protein_taxa=f"{config['output']['results_dir']}/taxa/protein_taxa.txt",
        rna_taxa=f"{config['output']['results_dir']}/taxa/rna_taxa.txt",
    output:
        cooccurrence=f"{config['output']['results_dir']}/analysis/cooccurrence.tsv",
        stats=f"{config['output']['results_dir']}/analysis/fisher_stats.tsv",
    run:
        from coevo.taxonomy.taxid_utils import read_taxids
        from coevo.analysis.cooccurrence import compute_cooccurrence
        from coevo.analysis.statistics import contingency_table, fisher_exact_test
        from coevo.io.result_writer import write_dict
        import os

        protein_taxa = read_taxids(input.protein_taxa)
        rna_taxa = read_taxids(input.rna_taxa)

        cooccurrence = compute_cooccurrence(protein_taxa, rna_taxa)
        os.makedirs(os.path.dirname(output.cooccurrence), exist_ok=True)
        write_dict(cooccurrence, output.cooccurrence)

        table = contingency_table(protein_taxa, rna_taxa)
        stats = fisher_exact_test(table)
        write_dict(stats, output.stats)


rule phylum_summary:
    """Group organisms by phylum and produce a summary table."""
    input:
        protein_taxa=f"{config['output']['results_dir']}/taxa/protein_taxa.txt",
    output:
        summary=f"{config['output']['results_dir']}/analysis/phylum_summary.tsv",
    run:
        from coevo.taxonomy.taxid_utils import read_taxids
        from coevo.analysis.phylogeny import phylum_summary
        from coevo.io.result_writer import write_dataframe
        import os

        protein_taxa = read_taxids(input.protein_taxa)
        summary = phylum_summary(protein_taxa)
        os.makedirs(os.path.dirname(output.summary), exist_ok=True)
        write_dataframe(summary, output.summary)
