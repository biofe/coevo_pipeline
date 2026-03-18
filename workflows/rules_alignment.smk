"""Snakemake rules for 16S rRNA alignment and motif detection."""


rule prepare_alignment:
    """Build alignment input FASTA and metadata from reference + BLAST hits."""
    input:
        query=config.get("rna_query", "examples/small_16s.fasta"),
        blast_hits=f"{config['output']['results_dir']}/blast/rna_hits.tsv",
    output:
        fasta=f"{config['output']['results_dir']}/alignment/16s_input.fasta",
        metadata=f"{config['output']['results_dir']}/alignment/16s_metadata.tsv",
    params:
        config_file=CONFIG_FILE,
    log:
        f"{config['output']['results_dir']}/logs/prepare_alignment.log",
    conda:
        CONDA_ENV
    shell:
        "coevo prepare-alignment {input.query} --config {params.config_file} > {log} 2>&1"


rule align_16s:
    """Align 16S rRNA sequences with MAFFT."""
    input:
        fasta=f"{config['output']['results_dir']}/alignment/16s_input.fasta",
    output:
        alignment=f"{config['output']['results_dir']}/alignment/16s_alignment.fasta",
    log:
        f"{config['output']['results_dir']}/logs/align_16s.log",
    conda:
        CONDA_ENV
    shell:
        """
        mafft --auto {input.fasta} > {output.alignment} 2> {log}
        """


rule detect_motif:
    """Detect the 16S rRNA structural motif in the alignment."""
    input:
        alignment=f"{config['output']['results_dir']}/alignment/16s_alignment.fasta",
    output:
        results=f"{config['output']['results_dir']}/analysis/motif_results.tsv",
    params:
        config_file=CONFIG_FILE,
    log:
        f"{config['output']['results_dir']}/logs/detect_motif.log",
    conda:
        CONDA_ENV
    shell:
        "coevo detect-motif --config {params.config_file} > {log} 2>&1"
