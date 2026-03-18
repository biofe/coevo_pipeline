"""FASTA file utilities using Biopython."""

from __future__ import annotations

from pathlib import Path
from typing import Generator

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
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


def build_alignment_fasta(
    reference_records: list[SeqRecord],
    blast_df: pd.DataFrame,
) -> tuple[list[SeqRecord], pd.DataFrame]:
    """Build alignment input FASTA from a reference sequence and BLAST hits.

    The reference sequence is placed first at local index 0.  Unique
    sequences extracted from *blast_df* (deduplicated by normalised sequence
    content) are appended starting at index 1.  Any BLAST sequence whose
    content is identical to the reference is not added as a separate record
    (so the reference is never duplicated), but its taxids are recorded under
    index 0 in the metadata.

    Parameters
    ----------
    reference_records:
        Sequence records from the reference FASTA file.  Only the **first**
        record is used (local index 0).
    blast_df:
        DataFrame returned by
        :func:`~coevo.blast.blast_parser.parse_blast_tabular`.  Must contain
        the ``sseq`` and ``staxids`` columns.

    Returns
    -------
    tuple[list[SeqRecord], pandas.DataFrame]
        - Ordered list of :class:`~Bio.SeqRecord.SeqRecord` objects: the
          reference followed by the unique BLAST sequences.
        - Metadata :class:`~pandas.DataFrame` with columns
          ``seq_index`` (int) and ``taxid`` (int), one row per
          *(sequence, organism)* pair.  When a BLAST hit carries a sequence
          identical to the reference, its taxids are linked to index 0 so
          that all organisms harbouring the reference sequence are
          represented in the metadata.

    Raises
    ------
    ValueError
        If *reference_records* is empty.
    """
    if not reference_records:
        raise ValueError("reference_records must not be empty")

    ref = reference_records[0]
    # Normalise the reference sequence the same way as BLAST hits so that a
    # BLAST result identical to the reference (modulo case and gap characters)
    # is correctly identified as a duplicate.
    ref_seq_norm = str(ref.seq).replace("-", "").upper()

    records: list[SeqRecord] = [
        SeqRecord(ref.seq, id=ref.id, description=ref.description)
    ]
    metadata_rows: list[dict] = []
    # Maps normalised sequence → its assigned local index so that when the
    # same sequence appears in multiple BLAST hits (different organisms) the
    # taxids from every hit are all linked to the single stored sequence.
    seen_seqs: dict[str, int] = {ref_seq_norm: 0}
    seq_idx = 1

    for _, row in blast_df.iterrows():
        # BLAST's sseq column contains the aligned subject sequence, which may
        # include gap characters.  Strip gaps before storing so that the output
        # FASTA contains raw (un-gapped) sequences suitable for re-alignment.
        seq_norm = str(row["sseq"]).replace("-", "").upper()

        if seq_norm in seen_seqs:
            existing_idx = seen_seqs[seq_norm]
            # Same sequence was already stored under *existing_idx* (which may
            # be 0 for the reference).  Still collect the taxids from this hit
            # so that all organisms sharing the sequence are represented in the
            # metadata.
            for raw_taxid in str(row["staxids"]).split(";"):
                raw_taxid = raw_taxid.strip()
                if raw_taxid and raw_taxid != "N/A":
                    try:
                        metadata_rows.append({"seq_index": existing_idx, "taxid": int(raw_taxid)})
                    except ValueError:
                        logger.warning(f"Could not parse taxid: {raw_taxid!r}")
            continue

        seen_seqs[seq_norm] = seq_idx
        records.append(SeqRecord(Seq(seq_norm), id=str(row["sseqid"]), description=""))

        for raw_taxid in str(row["staxids"]).split(";"):
            raw_taxid = raw_taxid.strip()
            if raw_taxid and raw_taxid != "N/A":
                try:
                    metadata_rows.append({"seq_index": seq_idx, "taxid": int(raw_taxid)})
                except ValueError:
                    logger.warning(f"Could not parse taxid: {raw_taxid!r}")

        seq_idx += 1

    metadata_df = (
        pd.DataFrame(metadata_rows, columns=["seq_index", "taxid"])
        .drop_duplicates()
        .reset_index(drop=True)
    )
    logger.info(
        f"Built alignment FASTA: 1 reference + {seq_idx - 1} unique blast sequences"
    )
    return records, metadata_df
