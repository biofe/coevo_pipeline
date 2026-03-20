"""Unit tests for coevo.analysis.phylogeny."""

from __future__ import annotations

from typing import Any
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from coevo.analysis.phylogeny import (
    CATEGORY_BOTH,
    CATEGORY_COLORS,
    CATEGORY_PROTEIN_ONLY,
    CATEGORY_RNA_ONLY,
    HAS_ETE4,
    _collapse_by_category,
    _dominant_category,
    _limit_visible_nodes,
    classify_taxids,
    draw_circular_tree,
    phylum_summary,
    top_phyla,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

PROTEIN_TAXA = {1, 2, 3}
RNA_TAXA = {3, 4, 5}
# 1, 2 -> protein_only; 4, 5 -> rna_only; 3 -> both


# ---------------------------------------------------------------------------
# Tests for classify_taxids
# ---------------------------------------------------------------------------


class TestClassifyTaxids:
    def test_returns_dict(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, RNA_TAXA)
        assert isinstance(result, dict)

    def test_protein_only(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, RNA_TAXA)
        assert result[1] == CATEGORY_PROTEIN_ONLY
        assert result[2] == CATEGORY_PROTEIN_ONLY

    def test_rna_only(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, RNA_TAXA)
        assert result[4] == CATEGORY_RNA_ONLY
        assert result[5] == CATEGORY_RNA_ONLY

    def test_both(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, RNA_TAXA)
        assert result[3] == CATEGORY_BOTH

    def test_all_taxids_present(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, RNA_TAXA)
        assert set(result.keys()) == PROTEIN_TAXA | RNA_TAXA

    def test_empty_protein(self) -> None:
        result = classify_taxids(set(), RNA_TAXA)
        for taxid in RNA_TAXA:
            assert result[taxid] == CATEGORY_RNA_ONLY

    def test_empty_rna(self) -> None:
        result = classify_taxids(PROTEIN_TAXA, set())
        for taxid in PROTEIN_TAXA:
            assert result[taxid] == CATEGORY_PROTEIN_ONLY

    def test_both_empty(self) -> None:
        result = classify_taxids(set(), set())
        assert result == {}

    def test_identical_sets(self) -> None:
        taxa = {10, 20, 30}
        result = classify_taxids(taxa, taxa)
        for taxid in taxa:
            assert result[taxid] == CATEGORY_BOTH

    def test_disjoint_sets(self) -> None:
        result = classify_taxids({1, 2}, {3, 4})
        assert result[1] == CATEGORY_PROTEIN_ONLY
        assert result[3] == CATEGORY_RNA_ONLY


# ---------------------------------------------------------------------------
# Tests for _dominant_category
# ---------------------------------------------------------------------------


class TestDominantCategory:
    def test_clear_dominant(self) -> None:
        cats = [CATEGORY_PROTEIN_ONLY] * 9 + [CATEGORY_RNA_ONLY]
        assert _dominant_category(cats) == CATEGORY_PROTEIN_ONLY

    def test_exact_threshold(self) -> None:
        # 9 out of 10 = 0.9 -> exactly meets threshold
        cats = [CATEGORY_BOTH] * 9 + [CATEGORY_RNA_ONLY]
        assert _dominant_category(cats, threshold=0.9) == CATEGORY_BOTH

    def test_below_threshold(self) -> None:
        # 8 out of 10 = 0.8 < 0.9
        cats = [CATEGORY_PROTEIN_ONLY] * 8 + [CATEGORY_RNA_ONLY] * 2
        assert _dominant_category(cats, threshold=0.9) is None

    def test_all_same(self) -> None:
        cats = [CATEGORY_RNA_ONLY, CATEGORY_RNA_ONLY, CATEGORY_RNA_ONLY]
        assert _dominant_category(cats) == CATEGORY_RNA_ONLY

    def test_empty_list(self) -> None:
        assert _dominant_category([]) is None

    def test_threshold_zero_picks_first_majority(self) -> None:
        # With threshold=0.0, any category is dominant (first one with count >= 0)
        cats = [CATEGORY_PROTEIN_ONLY, CATEGORY_RNA_ONLY]
        result = _dominant_category(cats, threshold=0.0)
        # Both have 50%, both are >= 0%, so a result is returned
        assert result in (CATEGORY_PROTEIN_ONLY, CATEGORY_RNA_ONLY)

    def test_single_element(self) -> None:
        assert _dominant_category([CATEGORY_BOTH]) == CATEGORY_BOTH


# ---------------------------------------------------------------------------
# Tests for phylum_summary
# ---------------------------------------------------------------------------


class TestPhylumSummary:
    def test_returns_dataframe(self) -> None:
        df = phylum_summary({1, 2, 3})
        assert isinstance(df, pd.DataFrame)

    def test_columns_present(self) -> None:
        df = phylum_summary({1, 2, 3})
        assert "phylum" in df.columns
        assert "count" in df.columns
        assert "fraction" in df.columns

    def test_unknown_phylum_when_no_map(self) -> None:
        df = phylum_summary({1, 2, 3})
        assert len(df) == 1
        assert df.iloc[0]["phylum"] == "Unknown"
        assert df.iloc[0]["count"] == 3

    def test_fraction_sums_to_one(self) -> None:
        phylum_map = {1: "A", 2: "B", 3: "A"}
        df = phylum_summary({1, 2, 3}, phylum_map=phylum_map)
        assert pytest.approx(df["fraction"].sum()) == 1.0

    def test_fraction_values(self) -> None:
        phylum_map = {1: "A", 2: "A", 3: "B", 4: "B"}
        df = phylum_summary({1, 2, 3, 4}, phylum_map=phylum_map)
        for _, row in df.iterrows():
            assert row["fraction"] == pytest.approx(0.5)

    def test_phylum_map_applied(self) -> None:
        phylum_map = {1: "Firmicutes", 2: "Proteobacteria", 3: "Firmicutes"}
        df = phylum_summary({1, 2, 3}, phylum_map=phylum_map)
        phyla = set(df["phylum"])
        assert "Firmicutes" in phyla
        assert "Proteobacteria" in phyla

    def test_sorted_by_count_descending(self) -> None:
        phylum_map = {1: "A", 2: "B", 3: "A", 4: "A"}
        df = phylum_summary({1, 2, 3, 4}, phylum_map=phylum_map)
        assert df.iloc[0]["phylum"] == "A"
        assert df.iloc[0]["count"] == 3

    def test_empty_taxids(self) -> None:
        df = phylum_summary(set())
        assert df.empty
        assert list(df.columns) == ["phylum", "count", "fraction"]

    def test_partial_phylum_map(self) -> None:
        phylum_map = {1: "Firmicutes"}
        df = phylum_summary({1, 2, 3}, phylum_map=phylum_map)
        phyla = set(df["phylum"])
        assert "Firmicutes" in phyla
        assert "Unknown" in phyla


# ---------------------------------------------------------------------------
# Tests for top_phyla
# ---------------------------------------------------------------------------


class TestTopPhyla:
    def test_returns_dataframe(self) -> None:
        df = top_phyla({1, 2, 3})
        assert isinstance(df, pd.DataFrame)

    def test_limits_to_n(self) -> None:
        phylum_map = {i: f"Phylum_{i}" for i in range(1, 21)}
        df = top_phyla(set(range(1, 21)), phylum_map=phylum_map, n=5)
        assert len(df) <= 5

    def test_default_n_10(self) -> None:
        phylum_map = {i: f"Phylum_{i}" for i in range(1, 8)}
        df = top_phyla(set(range(1, 8)), phylum_map=phylum_map)
        # Only 7 phyla exist, so all are returned
        assert len(df) == 7

    def test_n_larger_than_phyla(self) -> None:
        phylum_map = {1: "A", 2: "B"}
        df = top_phyla({1, 2}, phylum_map=phylum_map, n=100)
        assert len(df) == 2


# ---------------------------------------------------------------------------
# Tests for category constants and colours
# ---------------------------------------------------------------------------


class TestCategoryConstants:
    def test_all_categories_have_colours(self) -> None:
        for cat in (CATEGORY_PROTEIN_ONLY, CATEGORY_RNA_ONLY, CATEGORY_BOTH):
            assert cat in CATEGORY_COLORS

    def test_colours_are_hex_strings(self) -> None:
        for color in CATEGORY_COLORS.values():
            assert color.startswith("#")
            assert len(color) == 7


# ---------------------------------------------------------------------------
# Tests for draw_circular_tree (mocked ete4)
# ---------------------------------------------------------------------------


class TestDrawCircularTree:
    def test_raises_import_error_when_ete4_missing(self) -> None:
        """draw_circular_tree should raise ImportError when ete4 is absent."""
        with patch("coevo.analysis.phylogeny.HAS_ETE4", False):
            with pytest.raises(ImportError, match="ete4"):
                draw_circular_tree(PROTEIN_TAXA, RNA_TAXA)

    def test_raises_value_error_on_empty_taxids(self) -> None:
        """Both empty sets should raise ValueError."""
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            # We still need ete4 importable; patch the inner import too
            mock_ncbi = MagicMock()
            mock_tree = _make_mock_tree()
            mock_ncbi.return_value.get_topology.return_value = mock_tree
            mock_ncbi.return_value.get_taxid_translator.return_value = {}

            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                with pytest.raises(ValueError, match="taxid"):
                    draw_circular_tree(set(), set())

    def test_returns_tree_object(self) -> None:
        """draw_circular_tree should return the tree object."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {
            1: "Organism1",
            3: "Organism3",
        }

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                result = draw_circular_tree(
                    {1}, {3}, output_file="/tmp/test_tree.png"
                )
        assert result is mock_tree

    def test_calls_get_topology_with_union(self) -> None:
        """NCBITaxa.get_topology must be called with the union of taxids."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                draw_circular_tree({1, 2}, {3}, output_file="/tmp/t.png")

        call_args = mock_ncbi_cls.return_value.get_topology.call_args
        taxid_arg = call_args[0][0]
        assert set(taxid_arg) == {1, 2, 3}

    def test_motif_taxids_default_to_empty_set(self) -> None:
        """Omitting motif_taxids should not raise an error."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                # Should not raise
                draw_circular_tree({1}, {2}, output_file="/tmp/t.png")

    def test_interactive_show_called_without_output_file(self) -> None:
        """When output_file is None, tree.explore() should be called."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                draw_circular_tree({1}, {2})

        mock_tree.explore.assert_called_once()

    def test_render_called_with_output_file(self) -> None:
        """When output_file is given, tree.render() should be called."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                draw_circular_tree({1}, {2}, output_file="/tmp/out.png")

        mock_tree.render.assert_called_once()
        render_path = mock_tree.render.call_args[0][0]
        assert render_path == "/tmp/out.png"

    def test_logs_category_summary(self) -> None:
        """draw_circular_tree should log a category summary via logger.info."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        import coevo.analysis.phylogeny as phylo_mod

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                with patch.object(phylo_mod.logger, "info") as mock_info:
                    # protein_taxids={1}, rna_taxids={2} -> 1 protein_only, 1 rna_only
                    draw_circular_tree({1}, {2}, output_file="/tmp/out.png")

        # At least one info call should contain the category summary
        logged_messages = " ".join(
            str(call_args) for call_args in mock_info.call_args_list
        )
        assert "Protein Only" in logged_messages
        assert "Rna Only" in logged_messages
        assert "Both" in logged_messages


class TestCollapseByCategory:
    def test_collapses_when_dominant_category(self) -> None:
        """Internal node with all-same-category children should be collapsed."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_PROTEIN_ONLY),
                ("2", CATEGORY_PROTEIN_ONLY),
                ("3", CATEGORY_PROTEIN_ONLY),
            ],
        )
        parent.up = MagicMock()  # make this an internal (non-root) node
        _collapse_by_category(parent, threshold=0.9)
        # All children detached -> parent is now a leaf
        assert parent.is_leaf
        assert parent.category == CATEGORY_PROTEIN_ONLY

    def test_no_collapse_when_mixed(self) -> None:
        """Node with evenly mixed categories should not be collapsed."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_PROTEIN_ONLY),
                ("2", CATEGORY_PROTEIN_ONLY),
                ("3", CATEGORY_RNA_ONLY),
                ("4", CATEGORY_RNA_ONLY),
                ("5", CATEGORY_RNA_ONLY),
            ],
        )
        _collapse_by_category(parent, threshold=0.9)
        # 3/5 = 60% RNA_ONLY, 2/5 = 40% PROTEIN_ONLY -> no collapse
        assert not parent.is_leaf

    def test_custom_threshold(self) -> None:
        """Threshold of 0.5 should collapse when one category has > 50%."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_RNA_ONLY),
                ("2", CATEGORY_RNA_ONLY),
                ("3", CATEGORY_PROTEIN_ONLY),
            ],
        )
        parent.up = MagicMock()  # make this an internal (non-root) node
        _collapse_by_category(parent, threshold=0.5)
        # 2/3 ≈ 0.67 >= 0.5 -> collapse
        assert parent.is_leaf
        assert parent.category == CATEGORY_RNA_ONLY
    def test_root_never_collapsed(self) -> None:
        """Root node must never be collapsed even when all children share a category.

        Collapsing the root would remove all edges and leave a single filled
        circle – the bug reported in the issue.
        """
        root = _make_node_with_children(
            parent_name="root",
            children=[
                ("1", CATEGORY_RNA_ONLY),
                ("2", CATEGORY_RNA_ONLY),
                ("3", CATEGORY_RNA_ONLY),
            ],
        )
        # root.up is None so it IS the root
        _collapse_by_category(root, threshold=0.9)
        # Root should keep its children; tree must not reduce to a single node
        assert not root.is_leaf
        assert len(root.children) == 3





class TestLimitVisibleNodes:
    def test_no_change_when_already_below_limit(self) -> None:
        tree = _make_simple_tree(n_leaves=5)
        _limit_visible_nodes(tree, max_nodes=10)
        assert len(list(tree.iter_leaves())) == 5

    def test_collapses_to_max_nodes(self) -> None:
        tree = _make_simple_tree(n_leaves=20)
        _limit_visible_nodes(tree, max_nodes=5)
        assert len(list(tree.iter_leaves())) <= 5

    def test_exact_limit_unchanged(self) -> None:
        tree = _make_simple_tree(n_leaves=10)
        _limit_visible_nodes(tree, max_nodes=10)
        assert len(list(tree.iter_leaves())) == 10

    def test_no_change_when_below_limit_get_leaves(self) -> None:
        """PhyloTree-like nodes expose get_leaves() instead of iter_leaves()."""
        tree = _make_get_leaves_tree(n_leaves=5)
        _limit_visible_nodes(tree, max_nodes=10)
        assert len(tree.get_leaves()) == 5

    def test_collapses_to_max_nodes_get_leaves(self) -> None:
        """PhyloTree-like nodes expose get_leaves() instead of iter_leaves()."""
        tree = _make_get_leaves_tree(n_leaves=20)
        _limit_visible_nodes(tree, max_nodes=5)
        assert len(tree.get_leaves()) <= 5


# ---------------------------------------------------------------------------
# Helpers used by tests
# ---------------------------------------------------------------------------


def _make_mock_treeview() -> MagicMock:
    """Return a mock ete4.treeview module."""
    mock_tv = MagicMock()
    mock_ts = MagicMock()
    mock_tv.TreeStyle.return_value = mock_ts
    mock_tv.NodeStyle.return_value = MagicMock()
    mock_tv.TextFace.return_value = MagicMock()
    mock_tv.RectFace.return_value = MagicMock()
    return mock_tv


def _make_mock_tree() -> MagicMock:
    """Return a mock ete4 Tree suitable for draw_circular_tree tests."""
    mock_tree = MagicMock()

    # traverse yields an empty list by default so the annotation loop is skipped
    mock_tree.traverse.return_value = []
    mock_tree.iter_leaves.return_value = []

    return mock_tree


def _make_node_with_children(
    parent_name: str,
    children: list[tuple[str, str]],
) -> MagicMock:
    """Build a minimal mock ete4 tree for collapse tests.

    Parameters
    ----------
    parent_name:
        Name attribute of the parent node.
    children:
        List of ``(name, category)`` pairs for child nodes.

    Returns
    -------
    MagicMock
        A mock parent node whose ``.children``, ``.is_leaf``,
        ``.traverse()``, and ``.add_feature()`` behave realistically.
    """

    class SimpleNode:
        def __init__(self, name: str, category: str | None = None) -> None:
            self.name = name
            self.category = category
            self.children: list[SimpleNode] = []
            self.up: SimpleNode | None = None

        @property
        def is_leaf(self) -> bool:
            return len(self.children) == 0

        def iter_leaves(self):  # type: ignore[override]
            if self.is_leaf:
                yield self
            else:
                for child in self.children:
                    yield from child.iter_leaves()

        def traverse(self, order: str = "levelorder"):
            if order == "postorder":
                for child in self.children:
                    yield from child.traverse("postorder")
                yield self
            else:
                yield self
                for child in self.children:
                    yield from child.traverse(order)

        def add_feature(self, name: str, value: object) -> None:
            setattr(self, name, value)

        def detach(self) -> None:
            if self.up is not None:
                self.up.children.remove(self)
                self.up = None

    parent = SimpleNode(parent_name)
    for child_name, child_cat in children:
        child = SimpleNode(child_name, category=child_cat)
        child.up = parent
        parent.children.append(child)

    return parent  # type: ignore[return-value]


def _make_simple_tree(n_leaves: int) -> MagicMock:
    """Build a two-level tree (root -> n_leaves leaves) for limit tests."""

    class SimpleNode:
        def __init__(self, name: str, category: str | None = None) -> None:
            self.name = name
            self.category = category
            self.children: list[SimpleNode] = []
            self.up: SimpleNode | None = None

        @property
        def is_leaf(self) -> bool:
            return len(self.children) == 0

        def iter_leaves(self):  # type: ignore[override]
            if self.is_leaf:
                yield self
            else:
                for child in self.children:
                    yield from child.iter_leaves()

        def traverse(self, order: str = "levelorder"):
            yield self
            for child in self.children:
                yield from child.traverse(order)

        def add_feature(self, name: str, value: object) -> None:
            setattr(self, name, value)

        def detach(self) -> None:
            if self.up is not None:
                self.up.children.remove(self)
                self.up = None

        def get_distance(self, other: object, topology_only: bool = False) -> float:
            # Simple two-level tree: leaves have depth 1
            return 1.0 if self.is_leaf else 0.0

    root = SimpleNode("root")
    for i in range(n_leaves):
        leaf = SimpleNode(str(i), category=CATEGORY_BOTH)
        leaf.up = root
        root.children.append(leaf)

    return root  # type: ignore[return-value]


def _make_get_leaves_tree(n_leaves: int) -> Any:
    """Build a two-level tree that exposes ``get_leaves()`` but NOT ``iter_leaves()``.

    This mimics the ``PhyloTree`` API used by ete4, where ``get_leaves()`` is
    available but ``iter_leaves()`` is not.
    """

    class PhyloNode:
        def __init__(self, name: str, category: str | None = None) -> None:
            self.name = name
            self.category = category
            self.children: list[PhyloNode] = []
            self.up: PhyloNode | None = None

        @property
        def is_leaf(self) -> bool:
            return len(self.children) == 0

        def get_leaves(self) -> list[PhyloNode]:
            if self.is_leaf:
                return [self]
            result: list[PhyloNode] = []
            for child in self.children:
                result.extend(child.get_leaves())
            return result

        def traverse(self, order: str = "levelorder"):
            yield self
            for child in self.children:
                yield from child.traverse(order)

        def add_feature(self, name: str, value: object) -> None:
            setattr(self, name, value)

        def detach(self) -> None:
            if self.up is not None:
                self.up.children.remove(self)
                self.up = None

        def get_distance(self, other: object, topology_only: bool = False) -> float:
            # Simple two-level tree: leaves have depth 1
            return 1.0 if self.is_leaf else 0.0

    root = PhyloNode("root")
    for i in range(n_leaves):
        leaf = PhyloNode(str(i), category=CATEGORY_BOTH)
        leaf.up = root
        root.children.append(leaf)

    return root  # type: ignore[return-value]
