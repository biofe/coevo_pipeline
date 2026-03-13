"""Motif detection in multiple sequence alignments."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from loguru import logger


def detect_motif_in_alignment(
    alignment_file: str | Path,
    positions: list[int],
    residues: list[str],
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Detect a nucleotide/amino-acid motif in each sequence of an alignment.

    The motif is defined by a list of alignment column indices (*positions*)
    and the expected characters at those columns (*residues*).  A sequence is
    considered to have the motif when ALL specified positions match their
    expected residues (case-insensitive).

    Parameters
    ----------
    alignment_file:
        Path to the aligned FASTA file.
    positions:
        Zero-based column indices within the alignment to inspect.
    residues:
        Expected character at each position (same order as *positions*).
    n_jobs:
        Number of parallel worker processes.  Use ``1`` to disable
        multiprocessing (useful in tests).

    Returns
    -------
    pandas.DataFrame
        Table with columns ``sequence_id`` and ``motif_present`` (bool).

    Raises
    ------
    ValueError
        If *positions* and *residues* have different lengths.
    """
    if len(positions) != len(residues):
        raise ValueError(
            f"positions and residues must have the same length, "
            f"got {len(positions)} and {len(residues)}"
        )

    alignment_file = Path(alignment_file)
    logger.info(f"Loading alignment from {alignment_file}")
    alignment: MultipleSeqAlignment = AlignIO.read(str(alignment_file), "fasta")
    logger.info(f"Alignment has {len(alignment)} sequences, {alignment.get_alignment_length()} columns")

    records = list(alignment)

    args_list = [(str(rec.id), str(rec.seq), positions, residues) for rec in records]

    if n_jobs == 1:
        rows = [_check_motif(*args) for args in args_list]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            rows = list(pool.map(_check_motif_star, args_list))

    df = pd.DataFrame(rows, columns=["sequence_id", "motif_present"])
    n_present = df["motif_present"].sum()
    logger.info(
        f"Motif detected in {n_present}/{len(df)} sequences "
        f"(positions={positions}, residues={residues})"
    )
    return df


def _check_motif(
    seq_id: str,
    sequence: str,
    positions: list[int],
    residues: list[str],
) -> dict[str, Any]:
    """Check whether a single sequence contains the motif.

    Parameters
    ----------
    seq_id:
        Identifier of the sequence.
    sequence:
        Full aligned sequence string (may contain gap characters).
    positions:
        Alignment column indices to inspect.
    residues:
        Expected characters at each position.

    Returns
    -------
    dict
        ``{"sequence_id": str, "motif_present": bool}``
    """
    motif_present = all(
        pos < len(sequence) and sequence[pos].upper() == res.upper()
        for pos, res in zip(positions, residues)
    )
    return {"sequence_id": seq_id, "motif_present": motif_present}


def _check_motif_star(args: tuple) -> dict[str, Any]:
    """Unpack tuple and call :func:`_check_motif` (for multiprocessing)."""
    return _check_motif(*args)
