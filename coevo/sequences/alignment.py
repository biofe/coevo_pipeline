"""Multiple sequence alignment via external tools (MAFFT)."""

from __future__ import annotations

import subprocess
from pathlib import Path

from loguru import logger


def align_sequences_mafft(
    input_fasta: str | Path,
    output_fasta: str | Path,
    extra_args: list[str] | None = None,
) -> None:
    """Align sequences in *input_fasta* using MAFFT and write to *output_fasta*.

    MAFFT is called via :func:`subprocess.run`.  The ``--auto`` flag is used
    by default so that MAFFT selects the appropriate strategy based on the
    number and length of sequences.

    Parameters
    ----------
    input_fasta:
        Path to the unaligned FASTA file.
    output_fasta:
        Path where the aligned FASTA will be written.
    extra_args:
        Additional command-line arguments passed directly to MAFFT.

    Raises
    ------
    subprocess.CalledProcessError
        If MAFFT exits with a non-zero status.
    """
    input_fasta = Path(input_fasta)
    output_fasta = Path(output_fasta)
    output_fasta.parent.mkdir(parents=True, exist_ok=True)

    cmd = ["mafft", "--auto"]
    if extra_args:
        cmd.extend(extra_args)
    cmd.append(str(input_fasta))

    logger.info(f"Running MAFFT: {' '.join(cmd)}")

    with output_fasta.open("w") as out_fh:
        result = subprocess.run(
            cmd,
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True,
        )

    if result.returncode != 0:
        logger.error(f"MAFFT failed (exit {result.returncode}):\n{result.stderr}")
        raise subprocess.CalledProcessError(
            returncode=result.returncode,
            cmd=cmd,
            stderr=result.stderr,
        )

    logger.info(f"MAFFT alignment written to {output_fasta}")
    if result.stderr:
        logger.debug(f"MAFFT stderr: {result.stderr}")
