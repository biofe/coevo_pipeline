"""Snakemake rules for BLAST searches and taxid extraction."""


rule blast_protein:
    """Run BLASTP to identify protein family homologs."""
    input:
        query=config.get("protein_query", "examples/protein_query.fasta"),
    output:
        hits=f"{config['output']['results_dir']}/blast/protein_hits.tsv",
    threads:
        config["blast"]["threads"]
    params:
        db=config["blast"]["protein_db"],
        evalue=config["blast"]["evalue"],
        max_target_seqs=config["blast"]["max_target_seqs"],
    log:
        f"{config['output']['results_dir']}/logs/blast_protein.log",
    conda:
        CONDA_ENV
    shell:
        """
        blastp \
            -query {input.query} \
            -db {params.db} \
            -out {output.hits} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids" \
            > {log} 2>&1
        """


rule blast_16s:
    """Run BLASTN to find 16S rRNA sequences for relevant organisms."""
    input:
        query=config.get("rna_query", "examples/small_16s.fasta"),
    output:
        hits=f"{config['output']['results_dir']}/blast/rna_hits.tsv",
    threads:
        config["blast"]["threads"]
    params:
        db=config["blast"]["nucleotide_db"],
        evalue=config["blast"]["evalue"],
        max_target_seqs=config["blast"]["max_target_seqs"],
    log:
        f"{config['output']['results_dir']}/logs/blast_16s.log",
    conda:
        CONDA_ENV
    shell:
        """
        blastn \
            -query {input.query} \
            -db {params.db} \
            -out {output.hits} \
            -num_threads {threads} \
            -evalue {params.evalue} \
            -max_target_seqs {params.max_target_seqs} \
            -outfmt "6 qseqid sseqid pident length evalue bitscore staxids" \
            > {log} 2>&1
        """


rule extract_taxa:
    """Parse BLAST hits and extract unique taxids for protein and RNA results."""
    input:
        protein_hits=f"{config['output']['results_dir']}/blast/protein_hits.tsv",
        rna_hits=f"{config['output']['results_dir']}/blast/rna_hits.tsv",
    output:
        protein_taxids=f"{config['output']['results_dir']}/taxa/protein_taxa.txt",
        rna_taxids=f"{config['output']['results_dir']}/taxa/rna_taxa.txt",
    params:
        config_file=CONFIG_FILE,
    log:
        f"{config['output']['results_dir']}/logs/extract_taxa.log",
    conda:
        CONDA_ENV
    shell:
        "coevo extract-taxa --config {params.config_file} > {log} 2>&1"
