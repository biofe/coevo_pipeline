"""Snakemake rules for 16S rRNA alignment and motif detection."""


rule align_16s:
    """Align 16S rRNA sequences with MAFFT."""
    input:
        fasta=config.get("rna_query", "examples/small_16s.fasta"),
    output:
        alignment=f"{config['output']['results_dir']}/alignment/16s_alignment.fasta",
    log:
        f"{config['output']['results_dir']}/logs/align_16s.log",
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
        positions=config["motif"]["positions"],
        residues=config["motif"]["residues"],
    run:
        from coevo.sequences.motif_detection import detect_motif_in_alignment
        from coevo.io.result_writer import write_dataframe

        df = detect_motif_in_alignment(
            alignment_file=input.alignment,
            positions=params.positions,
            residues=params.residues,
        )
        import os
        os.makedirs(os.path.dirname(output.results), exist_ok=True)
        write_dataframe(df, output.results)
