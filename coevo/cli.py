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
    fragments = cfg["motif"]["fragments"]
    molecule_type = cfg["motif"].get("molecule_type", "rna")
    tolerance = cfg["motif"].get("tolerance", 0)

    results = detect_motif_in_alignment(
        alignment_file=str(alignment_file),
        positions=positions,
        fragments=fragments,
        molecule_type=molecule_type,
        tolerance=tolerance,
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
    from coevo.analysis.phylogeny import enterobacteriaceae_summary, motif_position_histogram
    from coevo.io.result_writer import write_dict, write_dataframe

    import pandas as pd

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

    entero_summary = enterobacteriaceae_summary(protein_taxa, rna_taxa)
    logger.info(f"Enterobacteriaceae summary: {len(entero_summary)} rows")

    # Load motif results and compute position histogram when available
    motif_file = analysis_dir / "motif_results.tsv"
    if motif_file.exists():
        motif_df = pd.read_csv(motif_file, sep="\t")
        motif_histogram = motif_position_histogram(motif_df)
        logger.info(f"Motif position histogram: {len(motif_histogram)} offset(s)")
    else:
        motif_histogram = pd.DataFrame(columns=["offset", "count"])
        logger.info("Motif results not found; skipping 16S position histogram")

    write_dataframe(entero_summary, analysis_dir / "enterobacteriaceae_summary.tsv")
    logger.info(
        f"Enterobacteriaceae summary written to {analysis_dir / 'enterobacteriaceae_summary.tsv'}"
    )

    write_dataframe(motif_histogram, analysis_dir / "motif_histogram.tsv")
    logger.info(f"Motif position histogram written to {analysis_dir / 'motif_histogram.tsv'}")



@app.command("draw-tree")
def draw_tree(
    config_path: str = CONFIG_OPTION,
    log_level: str = LOG_LEVEL_OPTION,
    output: Optional[str] = typer.Option(
        None,
        "--output",
        "-o",
        help=(
            "Path to save the rendered tree image (PNG or SVG). "
            "If omitted the tree is displayed interactively (requires a graphical display)."
        ),
    ),
    max_nodes: int = typer.Option(200, "--max-nodes", help="Maximum number of visible leaf nodes."),
    collapse_threshold: float = typer.Option(
        0.9,
        "--collapse-threshold",
        help=(
            "Fraction of leaf descendants sharing a category required to collapse a node (0–1). "
            "Higher values mean less collapsing and more nodes visible. Ignored with --show-all."
        ),
    ),
    show_all: bool = typer.Option(
        False,
        "--show-all",
        help="Display every leaf node without any collapsing, coloured by category.",
    ),
) -> None:
    """Draw a circular phylogenetic tree from BLAST taxonomy results.

    Reads ``protein_taxa.txt`` and ``rna_taxa.txt`` produced by the
    ``extract-taxa`` step.  When ``motif_results.tsv`` and
    ``16s_metadata.tsv`` are present, leaf nodes of organisms that carry the
    16S motif are highlighted with a gold background.

    Requires the optional ``ete4`` package::

        pip install ete4
    """
    cfg = _load(config_path, log_level)
    results_dir = get_results_dir(cfg)

    from coevo.taxonomy.taxid_utils import read_taxids
    from coevo.analysis.phylogeny import draw_circular_tree

    taxa_dir = results_dir / "taxa"
    protein_taxa = read_taxids(taxa_dir / "protein_taxa.txt")
    rna_taxa = read_taxids(taxa_dir / "rna_taxa.txt")

    # Optionally load motif taxids from motif results + alignment metadata
    motif_taxids: set[int] | None = None
    motif_file = results_dir / "analysis" / "motif_results.tsv"
    metadata_file = results_dir / "alignment" / "16s_metadata.tsv"
    if motif_file.exists() and metadata_file.exists():
        import pandas as pd

        motif_df = pd.read_csv(motif_file, sep="\t")
        metadata_df = pd.read_csv(metadata_file, sep="\t")
        # sequence_id in motif_results is "seq_N"; extract integer N to match
        # the seq_index column in the metadata table.
        motif_seq_ids = motif_df.loc[motif_df["motif_present"], "sequence_id"]
        motif_seq_indices: set[int] = set()
        for sid in motif_seq_ids:
            try:
                motif_seq_indices.add(int(str(sid).removeprefix("seq_")))
            except ValueError:
                logger.warning(f"Could not parse seq_index from sequence_id {sid!r}")
        motif_taxids = set(
            metadata_df.loc[
                metadata_df["seq_index"].isin(motif_seq_indices), "taxid"
            ]
        )
        logger.info(f"Loaded {len(motif_taxids)} motif taxids")
    else:
        logger.info("Motif results or metadata not found; proceeding without motif highlights")

    draw_circular_tree(
        protein_taxids=protein_taxa,
        rna_taxids=rna_taxa,
        motif_taxids=motif_taxids,
        output_file=output,
        max_nodes=max_nodes,
        collapse_threshold=collapse_threshold,
        show_all=show_all,
    )
    if output:
        logger.info(f"Phylogenetic tree saved to {output}")
    else:
        logger.info("Phylogenetic tree displayed interactively")


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
