"""BLAST runner: execute blastp and blastn via subprocess."""

from __future__ import annotations

import subprocess
from pathlib import Path

from loguru import logger


def run_blastp(
    query_fasta: str | Path,
    db: str,
    output_file: str | Path,
    threads: int = 8,
    evalue: float = 0.05,
    max_target_seqs: int = 50000,
) -> None:
    """Run BLASTP against a local protein database.

    Parameters
    ----------
    query_fasta:
        Path to the query protein FASTA file.
    db:
        Name or path of the local BLAST protein database (e.g. ``nr``).
    output_file:
        Path where tabular BLAST output will be written.
    threads:
        Number of CPU threads to use.
    evalue:
        E-value threshold for reporting hits.
    max_target_seqs:
        Maximum number of aligned sequences to keep.

    Raises
    ------
    subprocess.CalledProcessError
        If the blastp command exits with a non-zero status.
    """
    cmd = [
        "blastp",
        "-query", str(query_fasta),
        "-db", db,
        "-out", str(output_file),
        "-num_threads", str(threads),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore staxids sseq",
    ]
    _run_blast(cmd, label="BLASTP")


def run_blastn(
    query_fasta: str | Path,
    db: str,
    output_file: str | Path,
    threads: int = 8,
    evalue: float = 0.05,
    max_target_seqs: int = 50000,
) -> None:
    """Run BLASTN against a local nucleotide database.

    Parameters
    ----------
    query_fasta:
        Path to the query nucleotide FASTA file.
    db:
        Name or path of the local BLAST nucleotide database (e.g. ``core_nt``).
    output_file:
        Path where tabular BLAST output will be written.
    threads:
        Number of CPU threads to use.
    evalue:
        E-value threshold for reporting hits.
    max_target_seqs:
        Maximum number of aligned sequences to keep.

    Raises
    ------
    subprocess.CalledProcessError
        If the blastn command exits with a non-zero status.
    """
    cmd = [
        "blastn",
        "-query", str(query_fasta),
        "-db", db,
        "-out", str(output_file),
        "-num_threads", str(threads),
        "-evalue", str(evalue),
        "-max_target_seqs", str(max_target_seqs),
        "-outfmt", "6 qseqid sseqid pident length evalue bitscore staxids sseq",
    ]
    _run_blast(cmd, label="BLASTN")


def _run_blast(cmd: list[str], label: str) -> None:
    """Execute a BLAST command and raise on failure.

    Parameters
    ----------
    cmd:
        Command and arguments as a list of strings.
    label:
        Human-readable label for log messages (e.g. ``"BLASTP"``).
    """
    logger.info(f"Running {label}: {' '.join(cmd)}")
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        logger.error(f"{label} failed (exit {result.returncode}):\n{result.stderr}")
        raise subprocess.CalledProcessError(
            returncode=result.returncode,
            cmd=cmd,
            output=result.stdout,
            stderr=result.stderr,
        )
    logger.debug(f"{label} stdout: {result.stdout}")
    if result.stderr:
        logger.debug(f"{label} stderr: {result.stderr}")
    logger.info(f"{label} completed successfully")
