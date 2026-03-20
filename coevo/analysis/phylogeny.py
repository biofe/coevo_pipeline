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
    show_all: bool = False,
) -> Any:
    """Draw a circular phylogenetic tree from BLAST taxonomy IDs.

    Builds a tree of life for the union of *protein_taxids* and *rna_taxids*
    using :class:`ete4.NCBITaxa`.  Leaf nodes are colour-coded by their
    category of origin:

    - **green** (``#4CAF50``) – found only in the protein BLAST search.
    - **blue** (``#2196F3``) – found only in the rRNA BLAST search.
    - **orange** (``#FF9800``) – found in both searches.

    Unless *show_all* is ``True``, subtrees where at least *collapse_threshold*
    of all leaf descendants share the same category are collapsed into a single
    coloured leaf.  The collapse check uses the original leaf distribution so
    that a higher threshold always means less collapsing (more nodes visible).
    After category-based collapsing the tree is further limited to *max_nodes*
    visible leaves by collapsing nodes starting from the deepest subtrees first.

    When *show_all* is ``True``, no collapsing or node limiting is applied and
    every leaf is displayed, coloured by its category.

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
        When *None*, the tree is displayed interactively via ``tree.explore()``,
        which requires a graphical display.
    max_nodes:
        Maximum number of visible leaf nodes after collapsing.  Default 200.
        Ignored when *show_all* is ``True``.
    collapse_threshold:
        Fraction of leaf descendants that must share a category for a node to
        be collapsed into a single coloured leaf.  Higher values mean less
        collapsing (more nodes visible).  Default 0.9 (90 %).
        Ignored when *show_all* is ``True``.
    show_all:
        When ``True``, skip all collapsing and display every leaf node in the
        tree coloured by category.  Default ``False``.

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

    # Collapse subtrees where leaf descendants share a dominant category
    if not show_all:
        _collapse_by_category(tree, threshold=collapse_threshold)

        # Further limit visible nodes to max_nodes
        _limit_visible_nodes(tree, max_nodes=max_nodes)

    # Colour all remaining internal nodes by their subtree's dominant category
    _propagate_categories(tree)

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

    if output_file:
        # File rendering uses the PyQt6-backed treeview API (requires PyQt6).
        from ete4.treeview import (  # type: ignore[import]
            TreeStyle,
            TextFace,
            RectFace,
            faces,
        )

        def _layout(node: Any) -> None:
            taxid = getattr(node, "taxid", None)
            category = getattr(node, "category", None)
            sci_name = getattr(node, "sci_name", node.name)

            nstyle = node.img_style
            if category:
                color = CATEGORY_COLORS.get(category, DEFAULT_NODE_COLOR)
                nstyle["fgcolor"] = color
                nstyle["size"] = 5
            else:
                nstyle["fgcolor"] = DEFAULT_NODE_COLOR
                nstyle["size"] = 2

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

        for category, color in CATEGORY_COLORS.items():
            label = category.replace("_", " ").title()
            ts.legend.add_face(RectFace(20, 20, color, color), column=0)
            ts.legend.add_face(TextFace(f"  {label}", fsize=10), column=1)
        ts.legend_position = 1

        tree.render(str(output_file), tree_style=ts)
        logger.info(f"Circular phylogenetic tree rendered to {output_file}")
    else:
        # Interactive mode uses ETE4's smartview API which does not require
        # PyQt6 and serialises all data as plain JSON-compatible dicts.
        # ``keep_server=True`` makes the server thread non-daemon so the
        # process stays alive until the user terminates it.
        from ete4.smartview import Layout as _SmartLayout  # type: ignore[import]
        from ete4.smartview.faces import (  # type: ignore[import]
            LegendFace as _SmartLegendFace,
            TextFace as _SmartTextFace,
        )

        # Annotate each node's props for the click-popup info window.
        for node in tree.traverse():
            taxid = getattr(node, "taxid", None)
            if taxid is not None:
                node.props["Taxonomy ID"] = str(taxid)
            node_leaves = _get_leaves(node) if not node.is_leaf else [node]
            leaf_cats = [
                getattr(leaf, "category", None)
                for leaf in node_leaves
                if getattr(leaf, "category", None) is not None
            ]
            total = len(leaf_cats)
            if total > 0:
                node.props["Protein"] = (
                    f"{100 * leaf_cats.count(CATEGORY_PROTEIN_ONLY) / total:.0f}%"
                )
                node.props["16S rRNA"] = (
                    f"{100 * leaf_cats.count(CATEGORY_RNA_ONLY) / total:.0f}%"
                )
                node.props["Both"] = (
                    f"{100 * leaf_cats.count(CATEGORY_BOTH) / total:.0f}%"
                )

        def _draw_tree_fn(tree: Any) -> Any:
            yield {
                "shape": "circular",
                "collapsed-shape": "outline",
                "show-popup-props": ["Taxonomy ID", "Protein", "16S rRNA", "Both"],
            }
            yield _SmartLegendFace(
                "Category",
                "discrete",
                colormap={
                    cat.replace("_", " ").title(): color
                    for cat, color in CATEGORY_COLORS.items()
                },
            )

        def _draw_node(node: Any, collapsed: bool) -> Any:
            taxid = getattr(node, "taxid", None)
            category = getattr(node, "category", None)
            sci_name = getattr(node, "sci_name", node.name)

            color = (
                CATEGORY_COLORS.get(category, DEFAULT_NODE_COLOR)
                if category
                else DEFAULT_NODE_COLOR
            )

            # Style the node dot with the category colour.
            yield {"dot": {"fill": color, "r": 5 if category else 2}}

            # Add a text label for leaf nodes.
            if node.is_leaf or collapsed:
                label_bg = "#FFD700" if (taxid and taxid in motif_taxids) else None
                face_style = f"fill:{color}"
                if label_bg:
                    face_style += f";background-color:{label_bg}"
                yield _SmartTextFace(
                    sci_name, fs_max=9, style=face_style, position="right"
                )

        layout = _SmartLayout(
            name="coevo_layout",
            draw_tree=_draw_tree_fn,
            draw_node=_draw_node,
        )

        tree.explore(layouts=[layout], keep_server=True)

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
    """Collapse internal nodes whose leaf descendants share a dominant category.

    Traverses the tree in post-order so that children are evaluated before
    their parents.  The collapse decision is based on the distribution of
    *original* leaf descendants (computed before any collapsing takes place),
    so that a higher threshold always means less collapsing and more visible
    nodes.  When at least *threshold* of a non-root internal node's original
    leaf descendants share the same category, all children are detached (making
    the node a leaf) and the node is annotated with the dominant category.

    Parameters
    ----------
    tree:
        Root of an ete4 tree.
    threshold:
        Fraction of leaf descendants required to collapse.  Default 0.9.
        A higher value means less collapsing (more nodes remain visible).
    """
    # Pre-compute the original leaf-category list for every node *before* any
    # collapsing occurs.  This prevents the cascade where a collapsed subtree
    # makes its parent appear homogeneous, leading to unwanted further collapse.
    node_leaf_cats: dict[int, list[str]] = {}
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            cat = getattr(node, "category", None)
            node_leaf_cats[id(node)] = [cat] if cat is not None else []
        else:
            cats: list[str] = []
            for child in node.children:
                cats.extend(node_leaf_cats.get(id(child), []))
            node_leaf_cats[id(node)] = cats

    for node in list(tree.traverse("postorder")):
        if node.is_leaf:
            continue
        # Never collapse the root – doing so would leave the tree as a single
        # node (no edges, no visible structure).
        if node.up is None:
            continue
        leaf_cats = node_leaf_cats.get(id(node), [])
        dominant = _dominant_category(leaf_cats, threshold=threshold)
        if dominant is not None:
            node.category = dominant
            for child in list(node.children):
                child.detach()


def _node_depth(node: Any) -> int:
    """Return the number of edges from the root to *node*.

    Computed by traversing ``node.up`` links upward.  This avoids calling
    ``get_distance()`` which may not support the ``topology_only`` keyword
    argument in all ete4 tree types.

    Parameters
    ----------
    node:
        Any ete4 tree node with an ``up`` attribute.
    """
    depth = 0
    current = node
    while current.up is not None:
        depth += 1
        current = current.up
    return depth


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
        key=_node_depth,
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


def _propagate_categories(tree: Any) -> None:
    """Propagate category from leaf descendants to uncategorised internal nodes.

    After :func:`_collapse_by_category` and :func:`_limit_visible_nodes` some
    internal nodes may still lack a ``category`` attribute.  This function
    assigns to each such node the dominant category among its remaining leaf
    descendants so that intermediate nodes are coloured consistently with
    collapsed nodes.

    Parameters
    ----------
    tree:
        Root of an ete4 tree (modified in-place).
    """
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            continue
        if getattr(node, "category", None) is not None:
            continue
        leaf_cats = [
            getattr(leaf, "category", None)
            for leaf in _get_leaves(node)
            if getattr(leaf, "category", None) is not None
        ]
        dominant = _dominant_category(leaf_cats, threshold=0.0)
        if dominant is not None:
            node.category = dominant
