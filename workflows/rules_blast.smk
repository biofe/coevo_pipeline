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


rule extract_protein_taxids:
    """Parse BLASTP hits and extract unique taxids."""
    input:
        hits=f"{config['output']['results_dir']}/blast/protein_hits.tsv",
    output:
        taxids=f"{config['output']['results_dir']}/taxa/protein_taxa.txt",
    run:
        from coevo.blast.blast_parser import parse_blast_tabular, extract_taxids
        from coevo.taxonomy.taxid_utils import write_taxids

        df = parse_blast_tabular(input.hits)
        taxids = extract_taxids(df)
        write_taxids(taxids, output.taxids)


rule extract_rna_taxids:
    """Parse BLASTN hits and extract unique taxids."""
    input:
        hits=f"{config['output']['results_dir']}/blast/rna_hits.tsv",
    output:
        taxids=f"{config['output']['results_dir']}/taxa/rna_taxa.txt",
    run:
        from coevo.blast.blast_parser import parse_blast_tabular, extract_taxids
        from coevo.taxonomy.taxid_utils import write_taxids

        df = parse_blast_tabular(input.hits)
        taxids = extract_taxids(df)
        write_taxids(taxids, output.taxids)
