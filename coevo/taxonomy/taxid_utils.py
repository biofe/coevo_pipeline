"""Taxonomy utilities: reading and writing taxid sets."""

from __future__ import annotations

from pathlib import Path

from loguru import logger


def write_taxids(taxids: set[int], file_path: str | Path) -> None:
    """Write a set of integer taxids to a plain-text file (one per line).

    Parameters
    ----------
    taxids:
        Set of NCBI taxids to write.
    file_path:
        Destination file path.
    """
    file_path = Path(file_path)
    file_path.parent.mkdir(parents=True, exist_ok=True)
    with file_path.open("w") as fh:
        for taxid in sorted(taxids):
            fh.write(f"{taxid}\n")
    logger.info(f"Wrote {len(taxids)} taxids to {file_path}")


def read_taxids(file_path: str | Path) -> set[int]:
    """Read integer taxids from a plain-text file (one per line).

    Parameters
    ----------
    file_path:
        Path to the taxid file written by :func:`write_taxids`.

    Returns
    -------
    set[int]
        Taxids parsed from the file.

    Raises
    ------
    FileNotFoundError
        If *file_path* does not exist.
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Taxid file not found: {file_path}")

    taxids: set[int] = set()
    with file_path.open() as fh:
        for line in fh:
            line = line.strip()
            if line:
                try:
                    taxids.add(int(line))
                except ValueError:
                    logger.warning(f"Could not parse taxid line: {line!r}")
    logger.info(f"Read {len(taxids)} taxids from {file_path}")
    return taxids


def intersect_taxids(set_a: set[int], set_b: set[int]) -> set[int]:
    """Return the intersection of two taxid sets.

    Parameters
    ----------
    set_a:
        First set of taxids.
    set_b:
        Second set of taxids.

    Returns
    -------
    set[int]
        Taxids present in both sets.
    """
    result = set_a & set_b
    logger.debug(f"Taxid intersection: {len(result)} out of {len(set_a)} / {len(set_b)}")
    return result
