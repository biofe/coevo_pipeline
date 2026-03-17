"""BLAST output parser: parse tabular results and extract taxids."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
from loguru import logger

# Column order as specified by the BLAST -outfmt string used in blast_runner.py:
# "6 qseqid sseqid pident length evalue bitscore staxids sseq"
BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "evalue",
    "bitscore",
    "staxids",
    "sseq",
]

BLAST_DTYPES: dict[str, type] = {
    "qseqid": str,
    "sseqid": str,
    "pident": float,
    "length": int,
    "evalue": float,
    "bitscore": float,
    "staxids": str,
    "sseq": str,
}


def parse_blast_tabular(file_path: str | Path) -> pd.DataFrame:
    """Parse a BLAST tabular output file (format 6) into a DataFrame.

    The expected column order is::

        qseqid  sseqid  pident  length  evalue  bitscore  staxids

    Parameters
    ----------
    file_path:
        Path to the BLAST tabular output file.

    Returns
    -------
    pandas.DataFrame
        Parsed BLAST hits with typed columns.

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"BLAST output file not found: {file_path}")

    logger.info(f"Parsing BLAST tabular file: {file_path}")

    df = pd.read_csv(
        file_path,
        sep="\t",
        header=None,
        names=BLAST_COLUMNS,
        dtype=BLAST_DTYPES,
        comment="#",
    )
    logger.info(f"Parsed {len(df)} BLAST hits from {file_path}")
    return df


def extract_taxids(df: pd.DataFrame) -> set[int]:
    """Extract unique taxids from a parsed BLAST DataFrame.

    BLAST may report multiple taxids per hit separated by semicolons.
    All taxids are expanded and returned as a flat set of integers.

    Parameters
    ----------
    df:
        DataFrame returned by :func:`parse_blast_tabular`.

    Returns
    -------
    set[int]
        All unique taxids present in the BLAST results.
    """
    taxids: set[int] = set()
    for raw in df["staxids"].dropna():
        for part in str(raw).split(";"):
            part = part.strip()
            if part and part != "N/A":
                try:
                    taxids.add(int(part))
                except ValueError:
                    logger.warning(f"Could not parse taxid: {part!r}")
    logger.info(f"Extracted {len(taxids)} unique taxids")
    return taxids


def filter_blast_hits(
    df: pd.DataFrame,
    min_identity: float = 30.0,
    min_coverage: float = 0.0,
    query_length: int | None = None,
) -> pd.DataFrame:
    """Filter BLAST hits by percent identity and optional alignment coverage.

    Parameters
    ----------
    df:
        DataFrame returned by :func:`parse_blast_tabular`.
    min_identity:
        Minimum percent identity to retain a hit.
    min_coverage:
        Minimum alignment coverage (``length / query_length``) to retain a
        hit.  Ignored when *query_length* is ``None`` or ``0``.
    query_length:
        Length of the query sequence used to compute alignment coverage.

    Returns
    -------
    pandas.DataFrame
        Filtered DataFrame.
    """
    mask = df["pident"] >= min_identity
    if query_length and min_coverage > 0:
        mask &= (df["length"] / query_length) >= min_coverage
    filtered = df[mask].copy()
    logger.info(
        f"Filtered BLAST hits: {len(filtered)}/{len(df)} passed "
        f"(min_identity={min_identity}, min_coverage={min_coverage})"
    )
    return filtered


def deduplicate_blast_hits(df: pd.DataFrame) -> pd.DataFrame:
    """Deduplicate BLAST hits by accession ID and subject sequence.

    Rows sharing both the same subject accession (``sseqid``) and the same
    aligned subject sequence (``sseq``) are treated as the same hit.  For
    each such group, the row with the highest bitscore is retained.

    Deduplication based solely on taxonomy ID or sequence ID is intentionally
    avoided because multiple distinct sequences can share the same identifier
    (e.g. multiple 16S rRNA entries with identical seqid but different
    sequences).

    If the ``sseq`` column is absent (e.g. older output files), deduplication
    falls back to ``sseqid`` alone.

    Parameters
    ----------
    df:
        DataFrame returned by :func:`parse_blast_tabular`.

    Returns
    -------
    pandas.DataFrame
        Deduplicated DataFrame with reset index.
    """
    if df.empty:
        return df.copy()

    if "sseq" not in df.columns:
        dedup = (
            df.sort_values("bitscore", ascending=False)
            .drop_duplicates(subset=["sseqid"])
            .reset_index(drop=True)
        )
        logger.info(
            f"Deduplicated BLAST hits (by sseqid only): {len(dedup)}/{len(df)} retained"
        )
        return dedup

    dedup = (
        df.sort_values("bitscore", ascending=False)
        .drop_duplicates(subset=["sseqid", "sseq"])
        .reset_index(drop=True)
    )
    logger.info(
        f"Deduplicated BLAST hits (by sseqid + sseq): {len(dedup)}/{len(df)} retained"
    )
    return dedup
