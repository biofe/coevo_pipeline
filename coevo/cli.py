"""Command-line interface for the coevo pipeline using Typer."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from coevo.config import load_config, get_results_dir
from coevo.logging_utils import setup_logging

app = typer.Typer(
    name="coevo",
    help="Co-occurrence pipeline: protein family vs 16S rRNA structural motif.",
    add_completion=False,
)

CONFIG_OPTION = typer.Option("config.yaml", "--config", "-c", help="Path to config.yaml")
LOG_LEVEL_OPTION = typer.Option("INFO", "--log-level", help="Logging level")


def _load(config_path: str, log_level: str) -> dict:
    """Helper: set up logging and return config dict."""
    setup_logging(log_level=log_level)
    return load_config(config_path)


@app.command("blast-protein")
def blast_protein(
    query: str = typer.Argument(..., help="Path to protein query FASTA"),
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Run BLASTP to identify protein family homologs."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.blast.blast_runner import run_blastp

    output_file = results_dir / "blast" / "protein_hits.tsv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    run_blastp(
        query_fasta=query,
        db=cfg["blast"]["protein_db"],
        output_file=str(output_file),
        threads=cfg["blast"]["threads"],
        evalue=cfg["blast"]["evalue"],
        max_target_seqs=cfg["blast"]["max_target_seqs"],
    )
    logger.info(f"BLASTP results written to {output_file}")


@app.command("blast-16s")
def blast_16s(
    query: str = typer.Argument(..., help="Path to 16S rRNA query FASTA"),
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Run BLASTN to find 16S rRNA sequences for relevant organisms."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.blast.blast_runner import run_blastn

    output_file = results_dir / "blast" / "rna_hits.tsv"
    output_file.parent.mkdir(parents=True, exist_ok=True)

    run_blastn(
        query_fasta=query,
        db=cfg["blast"]["nucleotide_db"],
        output_file=str(output_file),
        threads=cfg["blast"]["threads"],
        evalue=cfg["blast"]["evalue"],
        max_target_seqs=cfg["blast"]["max_target_seqs"],
    )
    logger.info(f"BLASTN results written to {output_file}")


@app.command("extract-taxa")
def extract_taxa(
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Extract taxids from BLAST results and write to taxa files."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.blast.blast_parser import parse_blast_tabular, extract_taxids
    from coevo.taxonomy.taxid_utils import write_taxids

    taxa_dir = results_dir / "taxa"
    taxa_dir.mkdir(parents=True, exist_ok=True)

    for name, blast_file, out_file in [
        ("protein", results_dir / "blast" / "protein_hits.tsv", taxa_dir / "protein_taxa.txt"),
        ("rna", results_dir / "blast" / "rna_hits.tsv", taxa_dir / "rna_taxa.txt"),
    ]:
        df = parse_blast_tabular(blast_file)
        taxids = extract_taxids(df)
        write_taxids(taxids, out_file)
        logger.info(f"Wrote {len(taxids)} {name} taxids to {out_file}")


@app.command("prepare-alignment")
def prepare_alignment(
    query: str = typer.Argument(
        ..., help="Path to reference 16S rRNA FASTA (used as blast-16s input)"
    ),
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Build alignment input FASTA and sequence metadata from BLAST results.

    Reads the BLAST 16S hits produced by the ``blast-16s`` step, deduplicates
    subject sequences, and assembles an input FASTA that begins with the
    reference sequence followed by all unique BLAST hits.  A metadata TSV
    mapping each local sequence index to its organism (taxid) is written
    alongside the FASTA.
    """
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.blast.blast_parser import parse_blast_tabular, filter_blast_hits
    from coevo.sequences.fasta_utils import read_fasta, write_fasta, build_alignment_fasta
    from coevo.io.result_writer import write_dataframe

    blast_file = results_dir / "blast" / "rna_hits.tsv"
    alignment_dir = results_dir / "alignment"
    alignment_dir.mkdir(parents=True, exist_ok=True)

    reference_records = read_fasta(query)
    blast_df = parse_blast_tabular(blast_file)

    min_identity = cfg.get("filters", {}).get("min_identity", 0.0)
    blast_df = filter_blast_hits(blast_df, min_identity=min_identity)

    records, metadata = build_alignment_fasta(reference_records, blast_df)

    output_fasta = alignment_dir / "16s_input.fasta"
    write_fasta(records, output_fasta)
    logger.info(f"Alignment input FASTA written to {output_fasta}")

    metadata_file = alignment_dir / "16s_metadata.tsv"
    write_dataframe(metadata, metadata_file)
    logger.info(f"Sequence metadata written to {metadata_file}")


@app.command("align-16s")
def align_16s(
    input_fasta: str = typer.Argument(..., help="Path to 16S rRNA FASTA to align"),
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Align 16S rRNA sequences using MAFFT."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.sequences.alignment import align_sequences_mafft

    output_fasta = results_dir / "alignment" / "16s_alignment.fasta"
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    align_sequences_mafft(input_fasta=input_fasta, output_fasta=str(output_fasta))
    logger.info(f"Alignment written to {output_fasta}")


@app.command("detect-motif")
def detect_motif(
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Detect the 16S rRNA structural motif in the alignment."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.sequences.motif_detection import detect_motif_in_alignment
    from coevo.io.result_writer import write_dataframe

    alignment_file = results_dir / "alignment" / "16s_alignment.fasta"
    positions = cfg["motif"]["positions"]
    residues = cfg["motif"]["residues"]

    results = detect_motif_in_alignment(
        alignment_file=str(alignment_file),
        positions=positions,
        residues=residues,
    )

    out_file = results_dir / "analysis" / "motif_results.tsv"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    write_dataframe(results, out_file)
    logger.info(f"Motif detection results written to {out_file}")


@app.command("analyse")
def analyse(
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Compute co-occurrence statistics and produce summary tables."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.taxonomy.taxid_utils import read_taxids
    from coevo.analysis.cooccurrence import compute_cooccurrence
    from coevo.analysis.statistics import contingency_table, fisher_exact_test
    from coevo.analysis.phylogeny import phylum_summary
    from coevo.io.result_writer import write_dataframe, write_dict

    taxa_dir = results_dir / "taxa"
    analysis_dir = results_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    protein_taxa = read_taxids(taxa_dir / "protein_taxa.txt")
    rna_taxa = read_taxids(taxa_dir / "rna_taxa.txt")

    cooccurrence = compute_cooccurrence(protein_taxa, rna_taxa)
    write_dict(cooccurrence, analysis_dir / "cooccurrence.tsv")
    logger.info(f"Co-occurrence stats: {cooccurrence}")

    table = contingency_table(protein_taxa, rna_taxa)
    stats = fisher_exact_test(table)
    logger.info(f"Fisher exact test: {stats}")

    summary = phylum_summary(protein_taxa)
    write_dataframe(summary, analysis_dir / "phylum_summary.tsv")
    logger.info(f"Phylum summary written to {analysis_dir / 'phylum_summary.tsv'}")


@app.command("run-all")
def run_all(
    protein_query: str = typer.Argument(..., help="Path to protein query FASTA"),
    rna_query: str = typer.Argument(..., help="Path to 16S rRNA query FASTA"),
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
) -> None:
    """Run the full pipeline end-to-end."""
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    logger.info("Starting full pipeline run")
    blast_protein(query=protein_query, config_path=config_path, log_level=log_level)
    blast_16s(query=rna_query, config_path=config_path, log_level=log_level)
    extract_taxa(config_path=config_path, log_level=log_level)
    prepare_alignment(query=rna_query, config_path=config_path, log_level=log_level)
    prepared_fasta = str(results_dir / "alignment" / "16s_input.fasta")
    align_16s(input_fasta=prepared_fasta, config_path=config_path, log_level=log_level)
    detect_motif(config_path=config_path, log_level=log_level)
    analyse(config_path=config_path, log_level=log_level)
    logger.info("Full pipeline run complete")


if __name__ == "__main__":
    app()
