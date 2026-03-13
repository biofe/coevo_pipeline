"""FASTA file utilities using Biopython."""

from __future__ import annotations

from pathlib import Path
from typing import Generator

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from loguru import logger


def read_fasta(file_path: str | Path) -> list[SeqRecord]:
    """Read a FASTA file and return a list of SeqRecord objects.

    Parameters
    ----------
    file_path:
        Path to the FASTA file.

    Returns
    -------
    list[SeqRecord]
        Parsed sequence records.

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

    records = list(SeqIO.parse(str(file_path), "fasta"))
    logger.info(f"Read {len(records)} sequences from {file_path}")
    return records


def write_fasta(records: list[SeqRecord], file_path: str | Path) -> None:
    """Write a list of SeqRecord objects to a FASTA file.

    Parameters
    ----------
    records:
        Sequence records to write.
    file_path:
        Destination FASTA file path.
    """
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    count = SeqIO.write(records, str(file_path), "fasta")
    logger.info(f"Wrote {count} sequences to {file_path}")


def iter_fasta(file_path: str | Path) -> Generator[SeqRecord, None, None]:
    """Iterate over FASTA records without loading all into memory.

    Parameters
    ----------
    file_path:
        Path to the FASTA file.

    Yields
    ------
    SeqRecord
        Individual sequence records.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {file_path}")
    yield from SeqIO.parse(str(file_path), "fasta")


def filter_fasta_by_ids(
    file_path: str | Path,
    ids: set[str],
) -> list[SeqRecord]:
    """Return only those records whose IDs are in *ids*.

    Parameters
    ----------
    file_path:
        Path to the source FASTA file.
    ids:
        Set of sequence IDs to retain.

    Returns
    -------
    list[SeqRecord]
        Filtered sequence records.
    """
    selected = [rec for rec in iter_fasta(file_path) if rec.id in ids]
    logger.info(f"Filtered FASTA: retained {len(selected)}/{len(ids)} requested IDs")
    return selected
