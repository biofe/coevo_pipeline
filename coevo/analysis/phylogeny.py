"""Phylogenetic overview utilities: group organisms by phylum."""

from __future__ import annotations

from typing import Any

import pandas as pd
from loguru import logger

# TODO: phylogenetic tree colouring with ete3 (future extension)
# TODO: integrate ETE3 for proper tree-based grouping when available


def phylum_summary(
    taxids: set[int],
    phylum_map: dict[int, str] | None = None,
) -> pd.DataFrame:
    """Group organisms by phylum and return a summary table.

    For the MVP, phylum assignment is done via an optional *phylum_map*
    lookup.  When no map is provided, all organisms are placed in the
    ``"Unknown"`` phylum (a placeholder until NCBI Entrez or a local taxonomy
    database is integrated).

    Parameters
    ----------
    taxids:
        Set of NCBI taxids to summarise.
    phylum_map:
        Optional mapping from taxid (int) to phylum name (str).

    Returns
    -------
    pandas.DataFrame
        Table with columns ``phylum`` and ``count``, sorted descending by
        count.
    """
    if phylum_map is None:
        phylum_map = {}

    rows: list[dict[str, Any]] = []
    for taxid in taxids:
        phylum = phylum_map.get(taxid, "Unknown")
        rows.append({"taxid": taxid, "phylum": phylum})

    if not rows:
        logger.warning("phylum_summary called with empty taxid set")
        return pd.DataFrame(columns=["phylum", "count"])

    df = pd.DataFrame(rows)
    summary = (
        df.groupby("phylum")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    logger.info(f"Phylum summary: {len(summary)} phyla from {len(taxids)} taxids")
    return summary


def top_phyla(
    taxids: set[int],
    phylum_map: dict[int, str] | None = None,
    n: int = 10,
) -> pd.DataFrame:
    """Return the *n* most common phyla for a set of taxids.

    Parameters
    ----------
    taxids:
        Set of NCBI taxids.
    phylum_map:
        Optional mapping from taxid to phylum name.
    n:
        Number of top phyla to return.

    Returns
    -------
    pandas.DataFrame
        Top-*n* rows from :func:`phylum_summary`.
    """
    summary = phylum_summary(taxids, phylum_map=phylum_map)
    return summary.head(n)
