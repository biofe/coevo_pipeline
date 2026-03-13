"""Co-occurrence analysis: compute organism-level statistics."""

from __future__ import annotations

from typing import TypedDict

from loguru import logger


class CooccurrenceResult(TypedDict):
    """Typed dict for co-occurrence statistics."""

    protein_only: int
    rna_only: int
    intersection: int
    union: int
    jaccard_index: float
    protein_total: int
    rna_total: int


def compute_cooccurrence(
    protein_taxa: set[int],
    rna_taxa: set[int],
) -> CooccurrenceResult:
    """Compute co-occurrence statistics between two sets of taxids.

    A taxid is considered to show *co-occurrence* when it appears in BOTH
    *protein_taxa* (has a protein-family homolog) and *rna_taxa* (has the
    16S rRNA motif).

    Parameters
    ----------
    protein_taxa:
        Set of NCBI taxids from the protein BLAST search.
    rna_taxa:
        Set of NCBI taxids from the 16S rRNA motif detection step.

    Returns
    -------
    CooccurrenceResult
        Dictionary with the following keys:

        - ``protein_total`` – total protein taxids
        - ``rna_total`` – total RNA taxids
        - ``intersection`` – taxids in both sets
        - ``union`` – taxids in either set
        - ``jaccard_index`` – ``|intersection| / |union|`` (0 if union is empty)
        - ``protein_only`` – taxids only in protein set
        - ``rna_only`` – taxids only in RNA set
    """
    intersection = protein_taxa & rna_taxa
    union = protein_taxa | rna_taxa
    jaccard = len(intersection) / len(union) if union else 0.0

    result: CooccurrenceResult = {
        "protein_total": len(protein_taxa),
        "rna_total": len(rna_taxa),
        "intersection": len(intersection),
        "union": len(union),
        "jaccard_index": jaccard,
        "protein_only": len(protein_taxa - rna_taxa),
        "rna_only": len(rna_taxa - protein_taxa),
    }

    logger.info(
        f"Co-occurrence: intersection={result['intersection']}, "
        f"jaccard={result['jaccard_index']:.4f}"
    )
    return result
