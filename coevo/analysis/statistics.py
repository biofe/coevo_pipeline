"""Statistical tests for co-occurrence analysis."""

from __future__ import annotations

from typing import TypedDict

import numpy as np
from scipy.stats import fisher_exact

from loguru import logger


class ContingencyTable(TypedDict):
    """2×2 contingency table for Fisher's exact test."""

    both: int          # protein AND motif
    protein_only: int  # protein BUT NOT motif
    motif_only: int    # motif BUT NOT protein
    neither: int       # NEITHER protein NOR motif (requires total_taxa)


class FisherResult(TypedDict):
    """Result of Fisher's exact test."""

    odds_ratio: float
    p_value: float


def contingency_table(
    protein_taxa: set[int],
    rna_taxa: set[int],
    all_taxa: set[int] | None = None,
) -> ContingencyTable:
    """Build a 2×2 contingency table from taxid sets.

    Parameters
    ----------
    protein_taxa:
        Taxids of organisms that carry the protein family.
    rna_taxa:
        Taxids of organisms that carry the 16S rRNA motif.
    all_taxa:
        Complete universe of taxids.  If ``None``, the union of both sets is
        used as the universe, so ``neither`` will be ``0``.

    Returns
    -------
    ContingencyTable
        Counts for each cell of the 2×2 table.
    """
    if all_taxa is None:
        all_taxa = protein_taxa | rna_taxa

    both = len(protein_taxa & rna_taxa)
    protein_only = len(protein_taxa - rna_taxa)
    motif_only = len(rna_taxa - protein_taxa)
    neither = len(all_taxa - protein_taxa - rna_taxa)

    table: ContingencyTable = {
        "both": both,
        "protein_only": protein_only,
        "motif_only": motif_only,
        "neither": neither,
    }
    logger.debug(f"Contingency table: {table}")
    return table


def fisher_exact_test(table: ContingencyTable) -> FisherResult:
    """Apply Fisher's exact test to a 2×2 contingency table.

    The table is arranged as::

        [[both, motif_only],
         [protein_only, neither]]

    Parameters
    ----------
    table:
        Contingency table returned by :func:`contingency_table`.

    Returns
    -------
    FisherResult
        Dictionary with ``odds_ratio`` and ``p_value`` keys.
    """
    matrix = np.array(
        [[table["both"], table["motif_only"]],
         [table["protein_only"], table["neither"]]],
        dtype=int,
    )
    odds_ratio, p_value = fisher_exact(matrix)
    result: FisherResult = {"odds_ratio": float(odds_ratio), "p_value": float(p_value)}
    logger.info(f"Fisher exact test: OR={odds_ratio:.4f}, p={p_value:.4e}")
    return result
