"""Snakemake rules for co-occurrence analysis and summary statistics."""


rule analyse:
    """Compute co-occurrence statistics and phylum summary from taxid sets."""
    input:
        protein_taxa=f"{config['output']['results_dir']}/taxa/protein_taxa.txt",
        rna_taxa=f"{config['output']['results_dir']}/taxa/rna_taxa.txt",
    output:
        cooccurrence=f"{config['output']['results_dir']}/analysis/cooccurrence.tsv",
        phylum_summary=f"{config['output']['results_dir']}/analysis/phylum_summary.tsv",
    params:
        config_file=CONFIG_FILE,
    log:
        f"{config['output']['results_dir']}/logs/analyse.log",
    conda:
        CONDA_ENV
    shell:
        "coevo analyse --config {params.config_file} > {log} 2>&1"
