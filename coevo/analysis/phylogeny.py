"""Phylogenetic overview utilities: group organisms by phylum and draw trees.

Tree visualisation requires the optional ``ete4`` package (NCBI's fork of the
ETE toolkit).  Install it with::

    pip install ete4

or via the package's *phylo* extra::

    pip install 'coevo-rnase-16s[phylo]'
"""

from __future__ import annotations

from typing import Any

import pandas as pd
from loguru import logger

# ---------------------------------------------------------------------------
# Optional ete4 import
# ---------------------------------------------------------------------------

try:
    import ete4 as _ete4  # noqa: F401

    HAS_ETE4 = True
except ImportError:  # pragma: no cover
    HAS_ETE4 = False

# ---------------------------------------------------------------------------
# Category constants and colours
# ---------------------------------------------------------------------------

#: Taxid found only in the protein BLAST search.
CATEGORY_PROTEIN_ONLY: str = "protein_only"
#: Taxid found only in the rRNA BLAST search.
CATEGORY_RNA_ONLY: str = "rna_only"
#: Taxid found in both searches.
CATEGORY_BOTH: str = "both"

#: Hex colour assigned to each category for tree node styling.
CATEGORY_COLORS: dict[str, str] = {
    CATEGORY_PROTEIN_ONLY: "#4CAF50",  # green
    CATEGORY_RNA_ONLY: "#2196F3",       # blue
    CATEGORY_BOTH: "#FF9800",           # orange
}

#: Default (uncategorised) node colour.
DEFAULT_NODE_COLOR: str = "#999999"  # grey


# ---------------------------------------------------------------------------
# Public API – category classification
# ---------------------------------------------------------------------------


def classify_taxids(
    protein_taxids: set[int],
    rna_taxids: set[int],
) -> dict[int, str]:
    """Classify taxids by their presence in protein and rRNA BLAST results.

    Parameters
    ----------
    protein_taxids:
        Set of taxids identified in the protein BLAST search.
    rna_taxids:
        Set of taxids identified in the rRNA BLAST search.

    Returns
    -------
    dict[int, str]
        Mapping from taxid to one of:

        - ``"protein_only"`` – found only in the protein search.
        - ``"rna_only"`` – found only in the rRNA search.
        - ``"both"`` – found in both searches.
    """
    categories: dict[int, str] = {}
    for taxid in protein_taxids | rna_taxids:
        in_protein = taxid in protein_taxids
        in_rna = taxid in rna_taxids
        if in_protein and in_rna:
            categories[taxid] = CATEGORY_BOTH
        elif in_protein:
            categories[taxid] = CATEGORY_PROTEIN_ONLY
        else:
            categories[taxid] = CATEGORY_RNA_ONLY
    return categories


# ---------------------------------------------------------------------------
# Public API – phylum summary
# ---------------------------------------------------------------------------


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
        Table with columns ``phylum``, ``count``, and ``fraction``, sorted
        descending by count.
    """
    if phylum_map is None:
        phylum_map = {}

    rows: list[dict[str, Any]] = []
    for taxid in taxids:
        phylum = phylum_map.get(taxid, "Unknown")
        rows.append({"taxid": taxid, "phylum": phylum})

    if not rows:
        logger.warning("phylum_summary called with empty taxid set")
        return pd.DataFrame(columns=["phylum", "count", "fraction"])

    df = pd.DataFrame(rows)
    summary = (
        df.groupby("phylum")
        .size()
        .reset_index(name="count")
        .sort_values("count", ascending=False)
        .reset_index(drop=True)
    )
    total = summary["count"].sum()
    summary["fraction"] = summary["count"] / total
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


# ---------------------------------------------------------------------------
# Public API – circular tree visualisation
# ---------------------------------------------------------------------------


def draw_circular_tree(
    protein_taxids: set[int],
    rna_taxids: set[int],
    motif_taxids: set[int] | None = None,
    output_file: str | None = None,
    max_nodes: int = 200,
    collapse_threshold: float = 0.9,
) -> Any:
    """Draw a circular phylogenetic tree from BLAST taxonomy IDs.

    Builds a tree of life for the union of *protein_taxids* and *rna_taxids*
    using :class:`ete4.NCBITaxa`.  Leaf nodes are colour-coded by their
    category of origin:

    - **green** (``#4CAF50``) – found only in the protein BLAST search.
    - **blue** (``#2196F3``) – found only in the rRNA BLAST search.
    - **orange** (``#FF9800``) – found in both searches.

    If more than *collapse_threshold* (default: 90 %) of the direct children of
    an internal node share the same category, that node is collapsed and
    coloured accordingly.  The tree is limited to *max_nodes* visible leaves by
    collapsing nodes starting from those closest to the root (deepest subtrees
    first).

    Nodes whose taxid appears in *motif_taxids* receive a gold (``#FFD700``)
    background on their label to indicate motif presence.

    Parameters
    ----------
    protein_taxids:
        Taxids identified in the protein BLAST search.
    rna_taxids:
        Taxids identified in the rRNA BLAST search.
    motif_taxids:
        Optional set of taxids for which the 16S motif was detected.
        These nodes receive a highlighted label background.
    output_file:
        If provided, the tree is rendered to this image file (PNG or SVG).
        When *None*, the tree is displayed interactively via ``tree.show()``,
        which requires a graphical display.
    max_nodes:
        Maximum number of visible leaf nodes after collapsing.  Default 200.
    collapse_threshold:
        Fraction of direct children that must share a category for a node to
        be collapsed into a single coloured leaf.  Default 0.9 (90 %).

    Returns
    -------
    ete4.Tree
        The pruned and annotated ete4 :class:`~ete4.Tree` object.

    Raises
    ------
    ImportError
        If ``ete4`` is not installed.
    ValueError
        If both *protein_taxids* and *rna_taxids* are empty.
    """
    if not HAS_ETE4:
        raise ImportError(
            "ete4 is required for tree visualisation. "
            "Install it with: pip install ete4"
        )

    from ete4 import NCBITaxa  # type: ignore[import]
    from ete4.treeview import (  # type: ignore[import]
        TreeStyle,
        TextFace,
        RectFace,
        faces,
    )

    if motif_taxids is None:
        motif_taxids = set()

    union_taxids = protein_taxids | rna_taxids
    if not union_taxids:
        raise ValueError("At least one taxid must be provided")

    categories = classify_taxids(protein_taxids, rna_taxids)

    ncbi = NCBITaxa()
    taxid_list = sorted(union_taxids)
    tree = ncbi.get_topology(taxid_list, intermediate_nodes=True)

    # Retrieve scientific names for all nodes in the tree
    all_tree_taxids = [int(n.name) for n in tree.traverse() if n.name.isdigit()]
    sci_names = ncbi.get_taxid_translator(all_tree_taxids)

    for node in tree.traverse():
        taxid = int(node.name) if node.name.isdigit() else None
        node.taxid = taxid
        node.sci_name = sci_names.get(taxid, node.name) if taxid else node.name
        node.category = categories.get(taxid) if taxid else None

    # Collapse subtrees where children share a dominant category
    _collapse_by_category(tree, threshold=collapse_threshold)

    # Further limit visible nodes to max_nodes
    _limit_visible_nodes(tree, max_nodes=max_nodes)

    # Summarise category counts for the logger
    cat_counts: dict[str, int] = {
        CATEGORY_PROTEIN_ONLY: 0,
        CATEGORY_RNA_ONLY: 0,
        CATEGORY_BOTH: 0,
    }
    for cat in categories.values():
        if cat in cat_counts:
            cat_counts[cat] += 1
    logger.info(
        "Tree category summary — "
        f"Protein Only: {cat_counts[CATEGORY_PROTEIN_ONLY]}, "
        f"Rna Only: {cat_counts[CATEGORY_RNA_ONLY]}, "
        f"Both: {cat_counts[CATEGORY_BOTH]}"
    )

    # Build tree style (circular, custom layout)
    def _layout(node: Any) -> None:
        taxid = getattr(node, "taxid", None)
        category = getattr(node, "category", None)
        sci_name = getattr(node, "sci_name", node.name)

        # Use img_style directly so that PhyloTree nodes (which store style as
        # ``img_style`` rather than the private ``_img_style`` attribute) are
        # handled correctly and the interactive viewer does not raise an
        # AttributeError when inspecting node properties.
        nstyle = node.img_style
        if category:
            color = CATEGORY_COLORS.get(category, DEFAULT_NODE_COLOR)
            nstyle["fgcolor"] = color
            nstyle["size"] = 8
        else:
            nstyle["fgcolor"] = DEFAULT_NODE_COLOR
            nstyle["size"] = 4

        if node.is_leaf:
            label_bg = "#FFD700" if (taxid and taxid in motif_taxids) else None
            tf = TextFace(sci_name, fsize=9)
            if label_bg:
                tf.background.color = label_bg
            faces.add_face_to_node(tf, node, column=0)

    ts = TreeStyle()
    ts.mode = "c"  # circular layout
    ts.show_leaf_name = False
    ts.layout_fn = _layout
    ts.title.add_face(TextFace("Circular Phylogenetic Tree", fsize=14), column=0)

    # Legend
    for category, color in CATEGORY_COLORS.items():
        label = category.replace("_", " ").title()
        ts.legend.add_face(RectFace(20, 20, color, color), column=0)
        ts.legend.add_face(TextFace(f"  {label}", fsize=10), column=1)
    ts.legend_position = 1

    if output_file:
        tree.render(str(output_file), tree_style=ts)
        logger.info(f"Circular phylogenetic tree rendered to {output_file}")
    else:
        tree.show(tree_style=ts)

    return tree


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _get_leaves(node: Any) -> list[Any]:
    """Return leaf nodes for any ete4 tree type.

    Tries ``get_leaves()`` first (available on ``PhyloTree`` and modern ete4
    trees), then falls back to iterating via ``traverse()`` filtered by
    ``is_leaf``.

    Parameters
    ----------
    node:
        Any ete4 tree node.
    """
    if hasattr(node, "get_leaves"):
        return node.get_leaves()
    return [n for n in node.traverse() if n.is_leaf]


def _dominant_category(
    child_categories: list[str],
    threshold: float = 0.9,
) -> str | None:
    """Return the dominant category if one exceeds *threshold* fraction.

    Parameters
    ----------
    child_categories:
        Category label for each child/leaf to be considered.
    threshold:
        Minimum fraction required to declare a dominant category.

    Returns
    -------
    str | None
        The dominant category name, or ``None`` if no single category
        meets the threshold.
    """
    if not child_categories:
        return None
    counts: dict[str, int] = {}
    for cat in child_categories:
        counts[cat] = counts.get(cat, 0) + 1
    total = len(child_categories)
    for cat, count in counts.items():
        if count / total >= threshold:
            return cat
    return None


def _collapse_by_category(tree: Any, threshold: float = 0.9) -> None:
    """Collapse internal nodes whose direct children share a dominant category.

    Traverses the tree in post-order so that children are evaluated before
    their parents.  When more than *threshold* of an internal node's direct
    children share the same category, all children are detached (making the
    node a leaf) and the node is annotated with the dominant category.

    Parameters
    ----------
    tree:
        Root of an ete4 tree.
    threshold:
        Fraction of direct children required to collapse.  Default 0.9.
    """
    for node in list(tree.traverse("postorder")):
        if node.is_leaf:
            continue
        # Never collapse the root – doing so would leave the tree as a single
        # node (no edges, no visible structure).
        if node.up is None:
            continue
        child_cats = [
            getattr(child, "category", None)
            for child in node.children
            if getattr(child, "category", None) is not None
        ]
        dominant = _dominant_category(child_cats, threshold=threshold)
        if dominant is not None:
            node.category = dominant
            for child in list(node.children):
                child.detach()


def _limit_visible_nodes(tree: Any, max_nodes: int) -> None:
    """Collapse nodes until the leaf count is at most *max_nodes*.

    Nodes are selected for collapsing in order of decreasing depth (deepest
    subtrees first) so that the coarsest possible summary is retained when the
    node budget is exhausted.

    Parameters
    ----------
    tree:
        Root of an ete4 tree (modified in-place).
    max_nodes:
        Target maximum number of leaf nodes.
    """
    leaves = _get_leaves(tree)
    if len(leaves) <= max_nodes:
        return

    # Sort leaves by depth descending; collapse their parents first
    leaves_sorted = sorted(
        leaves,
        key=lambda n: n.get_distance(tree, topology_only=True),
        reverse=True,
    )

    for leaf in leaves_sorted:
        if len(_get_leaves(tree)) <= max_nodes:
            break
        parent = leaf.up
        if parent is None:
            continue
        # Determine a representative category for the collapsed node
        sib_cats = [
            getattr(sib, "category", None)
            for sib in _get_leaves(parent)
            if getattr(sib, "category", None) is not None
        ]
        dominant = _dominant_category(sib_cats, threshold=0.0) or CATEGORY_BOTH
        parent.category = dominant
        for child in list(parent.children):
            child.detach()
