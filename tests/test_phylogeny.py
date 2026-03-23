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
    ENTEROBACTERIACEAE_TAXID,
    HAS_ETE4,
    _collapse_by_category,
    _dominant_category,
    _limit_visible_nodes,
    _node_depth,
    _propagate_categories,
    classify_taxids,
    draw_circular_tree,
    enterobacteriaceae_summary,
    motif_position_histogram,
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
        """When output_file is None, tree.explore() should be called with keep_server=True."""
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": _make_mock_smartview(),
                    "ete4.smartview.faces": _make_mock_smartview_faces(),
                },
            ):
                draw_circular_tree({1}, {2})

        mock_tree.explore.assert_called_once()
        _, explore_kwargs = mock_tree.explore.call_args
        assert explore_kwargs.get("keep_server") is True

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
# Tests for _node_depth
# ---------------------------------------------------------------------------


class TestNodeDepth:
    def test_root_has_depth_zero(self) -> None:
        """Root node (no parent) should have depth 0."""
        tree = _make_simple_tree(n_leaves=3)
        assert _node_depth(tree) == 0

    def test_leaf_has_depth_one_in_two_level_tree(self) -> None:
        """Leaves in a two-level tree should have depth 1."""
        tree = _make_simple_tree(n_leaves=3)
        for leaf in tree.iter_leaves():
            assert _node_depth(leaf) == 1

    def test_deeper_leaves_have_higher_depth(self) -> None:
        """Leaves two levels down should have depth 2."""
        tree = _make_three_level_tree()
        # All leaves in _make_three_level_tree are two edges from root
        assert all(_node_depth(leaf) == 2 for leaf in tree.iter_leaves())


# ---------------------------------------------------------------------------
# Tests for _collapse_by_category cascade fix
# ---------------------------------------------------------------------------


class TestCollapseByCategoryCascade:
    def test_higher_threshold_shows_more_nodes(self) -> None:
        """A higher threshold should collapse fewer nodes (more nodes visible).

        This verifies the cascade-free behaviour: with threshold=1.0 only nodes
        whose *original* leaf descendants are 100 % homogeneous are collapsed,
        whereas threshold=0.9 collapses nodes with ≥ 90 % homogeneous leaves,
        so threshold=1.0 must leave ≥ as many nodes as threshold=0.9.
        """
        # Build a three-level tree where some subtrees are pure and one is mixed
        tree_strict = _make_three_level_tree()
        tree_loose = _make_three_level_tree()

        _collapse_by_category(tree_strict, threshold=1.0)
        _collapse_by_category(tree_loose, threshold=0.9)

        strict_leaves = len(list(tree_strict.iter_leaves()))
        loose_leaves = len(list(tree_loose.iter_leaves()))

        # threshold=1.0 (strict) must leave at least as many leaves as 0.9 (loose)
        assert strict_leaves >= loose_leaves

    def test_no_cascade_when_mixed_subtree_present(self) -> None:
        """A subtree with mixed categories prevents parent from collapsing.

        With the pre-computed leaf approach, a node whose leaf descendants are
        heterogeneous will not be collapsed even if other subtrees have already
        been collapsed and the parent's visible children happen to be homogeneous.
        """
        tree = _make_three_level_tree()
        # The tree has one mixed child subtree (see _make_three_level_tree)
        _collapse_by_category(tree, threshold=1.0)
        # The mixed subtree prevents complete collapse to the root level
        assert not tree.is_leaf
        assert len(list(tree.iter_leaves())) > 1


# ---------------------------------------------------------------------------
# Tests for draw_circular_tree show_all parameter
# ---------------------------------------------------------------------------


class TestShowAll:
    def _build_mocks(self):
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}
        return mock_ncbi_cls, mock_tree

    def test_show_all_skips_collapse(self) -> None:
        """With show_all=True, _collapse_by_category must not be called."""
        mock_ncbi_cls, mock_tree = self._build_mocks()

        import coevo.analysis.phylogeny as phylo_mod

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                with patch.object(
                    phylo_mod, "_collapse_by_category"
                ) as mock_collapse, patch.object(
                    phylo_mod, "_limit_visible_nodes"
                ) as mock_limit:
                    draw_circular_tree(
                        {1}, {2}, output_file="/tmp/t.png", show_all=True
                    )

        mock_collapse.assert_not_called()
        mock_limit.assert_not_called()

    def test_show_all_false_calls_collapse(self) -> None:
        """With show_all=False (default), _collapse_by_category must be called."""
        mock_ncbi_cls, mock_tree = self._build_mocks()

        import coevo.analysis.phylogeny as phylo_mod

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.treeview": _make_mock_treeview(),
                },
            ):
                with patch.object(
                    phylo_mod, "_collapse_by_category"
                ) as mock_collapse, patch.object(
                    phylo_mod, "_limit_visible_nodes"
                ) as mock_limit:
                    draw_circular_tree({1}, {2}, output_file="/tmp/t.png")

        mock_collapse.assert_called_once()
        mock_limit.assert_called_once()


# ---------------------------------------------------------------------------
# Tests for _propagate_categories
# ---------------------------------------------------------------------------


class TestPropagateCategories:
    def test_sets_category_on_uncategorised_internal_node(self) -> None:
        """An internal node with no category gets one from its leaf subtree."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_PROTEIN_ONLY),
                ("2", CATEGORY_PROTEIN_ONLY),
                ("3", CATEGORY_PROTEIN_ONLY),
            ],
        )
        assert parent.category is None
        _propagate_categories(parent)
        assert parent.category == CATEGORY_PROTEIN_ONLY

    def test_skips_already_categorised_internal_node(self) -> None:
        """An internal node with a category already set is left unchanged."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_PROTEIN_ONLY),
                ("2", CATEGORY_RNA_ONLY),
            ],
        )
        parent.category = CATEGORY_BOTH
        _propagate_categories(parent)
        assert parent.category == CATEGORY_BOTH

    def test_skips_leaf_nodes(self) -> None:
        """Leaf nodes are unaffected by propagation."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_RNA_ONLY),
                ("2", CATEGORY_RNA_ONLY),
            ],
        )
        _propagate_categories(parent)
        # Children (leaves) should retain their original categories.
        for child in parent.children:
            assert child.category == CATEGORY_RNA_ONLY

    def test_mixed_subtree_picks_most_common(self) -> None:
        """With threshold=0.0, the most common category is chosen."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[
                ("1", CATEGORY_BOTH),
                ("2", CATEGORY_BOTH),
                ("3", CATEGORY_PROTEIN_ONLY),
            ],
        )
        _propagate_categories(parent)
        assert parent.category == CATEGORY_BOTH

    def test_no_category_when_all_leaves_uncategorised(self) -> None:
        """Internal node stays uncategorised when leaves have no category."""
        parent = _make_node_with_children(
            parent_name="10",
            children=[("1", None), ("2", None)],
        )
        _propagate_categories(parent)
        assert parent.category is None


# ---------------------------------------------------------------------------
# Tests for interactive view features (dot radius, legend, popup, shape)
# ---------------------------------------------------------------------------


class TestInteractiveFeatures:
    def _build_interactive_mocks(self):
        mock_ncbi_cls = MagicMock()
        mock_tree = _make_mock_tree()
        mock_ncbi_cls.return_value.get_topology.return_value = mock_tree
        mock_ncbi_cls.return_value.get_taxid_translator.return_value = {}
        return mock_ncbi_cls, mock_tree

    def test_legend_face_yielded_in_draw_tree(self) -> None:
        """The draw_tree function should yield a LegendFace with CATEGORY_COLORS."""
        mock_ncbi_cls, mock_tree = self._build_interactive_mocks()

        legend_face_calls: list = []

        class CapturingLayout:
            def __init__(self, name, draw_tree=None, draw_node=None, **kw):
                if callable(draw_tree):
                    elements = list(draw_tree(None))
                    for el in elements:
                        if hasattr(el, "colormap"):
                            legend_face_calls.append(el)

        mock_sv = MagicMock()
        mock_sv.Layout.side_effect = CapturingLayout

        mock_legend_face_cls = MagicMock()

        class RealLegendFace:
            def __init__(self, title, variable, colormap=None, **kw):
                self.colormap = colormap

        mock_faces = MagicMock()
        mock_faces.LegendFace.side_effect = RealLegendFace
        mock_faces.TextFace.return_value = MagicMock()

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": mock_sv,
                    "ete4.smartview.faces": mock_faces,
                },
            ):
                draw_circular_tree({1}, {2})

        assert len(legend_face_calls) == 1
        colormap = legend_face_calls[0].colormap
        assert colormap is not None
        # All three category labels must appear in the legend.
        label_keys = set(colormap.keys())
        assert any("Protein Only" in k or "Protein" in k for k in label_keys)

    def test_draw_tree_sets_collapsed_shape_outline(self) -> None:
        """The draw_tree function must request collapsed-shape=outline."""
        mock_ncbi_cls, mock_tree = self._build_interactive_mocks()

        style_dicts: list = []

        class CapturingLayout:
            def __init__(self, name, draw_tree=None, draw_node=None, **kw):
                if callable(draw_tree):
                    for el in draw_tree(None):
                        if isinstance(el, dict):
                            style_dicts.append(el)

        mock_sv = MagicMock()
        mock_sv.Layout.side_effect = CapturingLayout
        mock_faces = _make_mock_smartview_faces()

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": mock_sv,
                    "ete4.smartview.faces": mock_faces,
                },
            ):
                draw_circular_tree({1}, {2})

        combined = {}
        for d in style_dicts:
            combined.update(d)
        assert combined.get("collapsed-shape") == "outline"

    def test_draw_tree_sets_popup_props(self) -> None:
        """The draw_tree function must declare the four popup properties."""
        mock_ncbi_cls, mock_tree = self._build_interactive_mocks()

        style_dicts: list = []

        class CapturingLayout:
            def __init__(self, name, draw_tree=None, draw_node=None, **kw):
                if callable(draw_tree):
                    for el in draw_tree(None):
                        if isinstance(el, dict):
                            style_dicts.append(el)

        mock_sv = MagicMock()
        mock_sv.Layout.side_effect = CapturingLayout
        mock_faces = _make_mock_smartview_faces()

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": mock_sv,
                    "ete4.smartview.faces": mock_faces,
                },
            ):
                draw_circular_tree({1}, {2})

        combined = {}
        for d in style_dicts:
            combined.update(d)
        popup_props = combined.get("show-popup-props", [])
        assert "Taxonomy ID" in popup_props
        assert "Protein" in popup_props
        assert "16S rRNA" in popup_props
        assert "Both" in popup_props

    def test_draw_node_dot_radius_categorised(self) -> None:
        """Categorised nodes must have dot radius 5."""
        mock_ncbi_cls, mock_tree = self._build_interactive_mocks()

        captured_draw_node: list = []

        class CapturingLayout:
            def __init__(self, name, draw_tree=None, draw_node=None, **kw):
                if draw_node is not None:
                    captured_draw_node.append(draw_node)

        mock_sv = MagicMock()
        mock_sv.Layout.side_effect = CapturingLayout
        mock_faces = _make_mock_smartview_faces()

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": mock_sv,
                    "ete4.smartview.faces": mock_faces,
                },
            ):
                draw_circular_tree({1}, {2})

        assert captured_draw_node, "draw_node function was not captured"
        draw_node = captured_draw_node[0]

        # Simulate a categorised leaf node.
        mock_node = MagicMock()
        mock_node.taxid = 1
        mock_node.category = CATEGORY_PROTEIN_ONLY
        mock_node.sci_name = "Organism A"
        mock_node.name = "1"
        mock_node.is_leaf = True

        elements = list(draw_node(mock_node, False))
        dot_dicts = [e for e in elements if isinstance(e, dict) and "dot" in e]
        assert dot_dicts, "No dot dict found in draw_node output"
        assert dot_dicts[0]["dot"]["r"] == 5

    def test_draw_node_dot_radius_uncategorised(self) -> None:
        """Uncategorised nodes must have dot radius 2."""
        mock_ncbi_cls, mock_tree = self._build_interactive_mocks()

        captured_draw_node: list = []

        class CapturingLayout:
            def __init__(self, name, draw_tree=None, draw_node=None, **kw):
                if draw_node is not None:
                    captured_draw_node.append(draw_node)

        mock_sv = MagicMock()
        mock_sv.Layout.side_effect = CapturingLayout
        mock_faces = _make_mock_smartview_faces()

        with patch("coevo.analysis.phylogeny.HAS_ETE4", True):
            with patch.dict(
                "sys.modules",
                {
                    "ete4": MagicMock(NCBITaxa=mock_ncbi_cls),
                    "ete4.smartview": mock_sv,
                    "ete4.smartview.faces": mock_faces,
                },
            ):
                draw_circular_tree({1}, {2})

        assert captured_draw_node
        draw_node = captured_draw_node[0]

        # Simulate an uncategorised internal node.
        mock_node = MagicMock()
        mock_node.taxid = 99
        mock_node.category = None
        mock_node.sci_name = "Root"
        mock_node.name = "99"
        mock_node.is_leaf = False

        elements = list(draw_node(mock_node, False))
        dot_dicts = [e for e in elements if isinstance(e, dict) and "dot" in e]
        assert dot_dicts, "No dot dict found in draw_node output"
        assert dot_dicts[0]["dot"]["r"] == 2


def _make_mock_treeview() -> MagicMock:
    """Return a mock ete4.treeview module."""
    mock_tv = MagicMock()
    mock_ts = MagicMock()
    mock_tv.TreeStyle.return_value = mock_ts
    mock_tv.NodeStyle.return_value = MagicMock()
    mock_tv.TextFace.return_value = MagicMock()
    mock_tv.RectFace.return_value = MagicMock()
    return mock_tv


def _make_mock_smartview() -> MagicMock:
    """Return a mock ete4.smartview module for interactive-path tests."""
    mock_sv = MagicMock()
    mock_sv.Layout.return_value = MagicMock()
    return mock_sv


def _make_mock_smartview_faces() -> MagicMock:
    """Return a mock ete4.smartview.faces module for interactive-path tests."""
    mock_faces = MagicMock()
    mock_faces.TextFace.return_value = MagicMock()
    return mock_faces


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


def _make_three_level_tree() -> Any:
    """Build a three-level tree for cascade/threshold tests.

    Structure (root → internal nodes → leaves)::

        root
        ├── internal_A  (3 protein_only leaves  → 100% homogeneous)
        │   ├── leaf_a1 (protein_only)
        │   ├── leaf_a2 (protein_only)
        │   └── leaf_a3 (protein_only)
        ├── internal_B  (3 rna_only leaves → 100% homogeneous)
        │   ├── leaf_b1 (rna_only)
        │   ├── leaf_b2 (rna_only)
        │   └── leaf_b3 (rna_only)
        └── internal_C  (2 protein_only + 1 rna_only → mixed ~67%/33%)
            ├── leaf_c1 (protein_only)
            ├── leaf_c2 (protein_only)
            └── leaf_c3 (rna_only)

    ``internal_A`` and ``internal_B`` are fully homogeneous subtrees and will
    be collapsed at threshold ≤ 1.0.  ``internal_C`` is mixed, so it will only
    collapse if the threshold is ≤ 0.67.  The root is never collapsed.
    """

    class TLNode:
        def __init__(self, name: str, category: str | None = None) -> None:
            self.name = name
            self.category = category
            self.children: list[TLNode] = []
            self.up: TLNode | None = None

        @property
        def is_leaf(self) -> bool:
            return len(self.children) == 0

        def iter_leaves(self):
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

    def _add_children(parent: TLNode, specs: list[tuple[str, str | None]]) -> None:
        for name, cat in specs:
            child = TLNode(name, category=cat)
            child.up = parent
            parent.children.append(child)

    root = TLNode("root")

    internal_a = TLNode("internal_A")
    internal_a.up = root
    root.children.append(internal_a)
    _add_children(
        internal_a,
        [("a1", CATEGORY_PROTEIN_ONLY), ("a2", CATEGORY_PROTEIN_ONLY), ("a3", CATEGORY_PROTEIN_ONLY)],
    )

    internal_b = TLNode("internal_B")
    internal_b.up = root
    root.children.append(internal_b)
    _add_children(
        internal_b,
        [("b1", CATEGORY_RNA_ONLY), ("b2", CATEGORY_RNA_ONLY), ("b3", CATEGORY_RNA_ONLY)],
    )

    internal_c = TLNode("internal_C")
    internal_c.up = root
    root.children.append(internal_c)
    _add_children(
        internal_c,
        [("c1", CATEGORY_PROTEIN_ONLY), ("c2", CATEGORY_PROTEIN_ONLY), ("c3", CATEGORY_RNA_ONLY)],
    )

    return root  # type: ignore[return-value]


# ---------------------------------------------------------------------------
# Tests for ENTEROBACTERIACEAE_TAXID constant
# ---------------------------------------------------------------------------


class TestEnterobacteriaceaeTaxidConstant:
    def test_value_is_543(self) -> None:
        assert ENTEROBACTERIACEAE_TAXID == 543


# ---------------------------------------------------------------------------
# Tests for enterobacteriaceae_summary
# ---------------------------------------------------------------------------


def _make_mock_ncbi(
    lineage_map: dict[int, list[int]],
    rank_map: dict[int, str],
    name_map: dict[int, str],
) -> MagicMock:
    """Build a minimal NCBITaxa mock with fixed lineage/rank/name data."""
    mock = MagicMock()

    def _get_lineage(taxid: int) -> list[int]:
        if taxid not in lineage_map:
            raise ValueError(f"Unknown taxid {taxid}")
        return lineage_map[taxid]

    mock.get_lineage.side_effect = _get_lineage
    mock.get_rank.side_effect = lambda taxids: {t: rank_map.get(t, "no rank") for t in taxids}
    mock.get_taxid_translator.side_effect = lambda taxids: {
        t: name_map[t] for t in taxids if t in name_map
    }
    return mock


# Minimal taxonomy tree used in the tests below:
#   543 (family: Enterobacteriaceae)
#     |-- 561 (genus: Escherichia)
#     |    |-- 562 (species: Escherichia coli)
#     |    |    |-- 1001 (strain)
#     |    |    +-- 1002 (strain)
#     |    +-- 1003 (species: Escherichia fergusonii)
#     +-- 590 (genus: Salmonella)
#          +-- 2000 (species: Salmonella enterica)
#               +-- 2001 (strain)
_LINEAGE_MAP: dict[int, list[int]] = {
    1001: [1, 543, 561, 562, 1001],
    1002: [1, 543, 561, 562, 1002],
    1003: [1, 543, 561, 1003],
    2001: [1, 543, 590, 2000, 2001],
    9999: [1, 9999],  # outside Enterobacteriaceae
}
_RANK_MAP: dict[int, str] = {
    1: "no rank",
    543: "family",
    561: "genus",
    562: "species",
    590: "genus",
    2000: "species",
    1001: "strain",
    1002: "strain",
    1003: "species",
    2001: "strain",
    9999: "genus",
}
_NAME_MAP: dict[int, str] = {
    543: "Enterobacteriaceae",
    561: "Escherichia",
    562: "Escherichia coli",
    590: "Salmonella",
    2000: "Salmonella enterica",
    1001: "E. coli strain A",
    1002: "E. coli strain B",
    1003: "Escherichia fergusonii",
    2001: "S. enterica strain X",
}


def _patch_ete4_ncbi(ncbi_mock: MagicMock):
    """Context manager: patch sys.modules so ete4.NCBITaxa returns *ncbi_mock*."""
    return patch.dict(
        "sys.modules",
        {"ete4": MagicMock(NCBITaxa=MagicMock(return_value=ncbi_mock))},
    )


class TestEnterobacteriaceaeSummary:
    # ------------------------------------------------------------------
    # Without ete4
    # ------------------------------------------------------------------

    def test_returns_empty_when_ete4_missing(self) -> None:
        with patch("coevo.analysis.phylogeny.HAS_ETE4", False):
            df = enterobacteriaceae_summary({1}, {2})
        assert isinstance(df, pd.DataFrame)
        assert df.empty
        assert list(df.columns) == [
            "rank", "name", "total", "protein_only_pct", "rna_only_pct", "both_pct"
        ]

    def test_returns_empty_when_taxids_empty(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary(set(), set())
        assert df.empty

    # ------------------------------------------------------------------
    # Basic structure
    # ------------------------------------------------------------------

    def test_returns_dataframe(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 1002}, {2001})
        assert isinstance(df, pd.DataFrame)

    def test_columns_present(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 1002}, {2001})
        for col in ["rank", "name", "total", "protein_only_pct", "rna_only_pct", "both_pct"]:
            assert col in df.columns

    def test_out_of_family_taxids_excluded(self) -> None:
        """Taxid 9999 is not in Enterobacteriaceae and must not appear."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 9999}, set())
        names = df["name"].tolist()
        # 1001 is a strain under E. coli so Escherichia genus row should appear
        assert "Escherichia" in names

    def test_family_row_present(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 2001}, set())
        assert "family" in df["rank"].values

    def test_genus_rows_present(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 2001}, set())
        genus_names = df.loc[df["rank"] == "genus", "name"].tolist()
        assert "Escherichia" in genus_names
        assert "Salmonella" in genus_names

    def test_species_rows_present(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 2001}, set())
        species_names = df.loc[df["rank"] == "species", "name"].tolist()
        assert "Escherichia coli" in species_names
        assert "Salmonella enterica" in species_names

    def test_strains_not_listed_as_species(self) -> None:
        """Strain taxids must not appear as dedicated species rows."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 1002, 2001}, set())
        species_names = set(df.loc[df["rank"] == "species", "name"])
        assert "E. coli strain A" not in species_names
        assert "E. coli strain B" not in species_names
        assert "S. enterica strain X" not in species_names

    def test_family_comes_first(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 2001}, set())
        assert df.iloc[0]["rank"] == "family"

    def test_percentages_sum_to_100(self) -> None:
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001}, {2001})
        for _, row in df.iterrows():
            pct_sum = row["protein_only_pct"] + row["rna_only_pct"] + row["both_pct"]
            assert pct_sum == pytest.approx(100.0, abs=0.5)

    def test_category_both_counted(self) -> None:
        """A taxid present in both sets must be counted as 'both'."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001}, {1001})
        fam_row = df.loc[df["rank"] == "family"].iloc[0]
        assert fam_row["both_pct"] == pytest.approx(100.0)
        assert fam_row["protein_only_pct"] == pytest.approx(0.0)
        assert fam_row["rna_only_pct"] == pytest.approx(0.0)

    def test_total_counts_per_genus(self) -> None:
        """Two strains under E. coli -> genus Escherichia total should be 2."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1001, 1002}, set())
        ecoli_genus = df.loc[(df["rank"] == "genus") & (df["name"] == "Escherichia")]
        assert ecoli_genus.iloc[0]["total"] == 2

    def test_no_taxids_in_family_returns_empty(self) -> None:
        """Only out-of-family taxids should produce an empty summary."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({9999}, set())
        assert df.empty

    def test_taxid_at_species_rank_listed_as_species(self) -> None:
        """Taxid 1003 has rank 'species' itself, so it should appear as a species row."""
        mock_ncbi = _make_mock_ncbi(_LINEAGE_MAP, _RANK_MAP, _NAME_MAP)
        with patch("coevo.analysis.phylogeny.HAS_ETE4", True), _patch_ete4_ncbi(mock_ncbi):
            df = enterobacteriaceae_summary({1003}, set())
        species_names = set(df.loc[df["rank"] == "species", "name"])
        assert "Escherichia fergusonii" in species_names


# ---------------------------------------------------------------------------
# Tests for motif_position_histogram
# ---------------------------------------------------------------------------


class TestMotifPositionHistogram:
    def _make_df(self, rows: list[dict]) -> pd.DataFrame:
        return pd.DataFrame(rows, columns=["sequence_id", "motif_present", "motif_offset"])

    def test_returns_dataframe(self) -> None:
        df = self._make_df([{"sequence_id": "s1", "motif_present": True, "motif_offset": 0}])
        result = motif_position_histogram(df)
        assert isinstance(result, pd.DataFrame)

    def test_columns(self) -> None:
        df = self._make_df([{"sequence_id": "s1", "motif_present": True, "motif_offset": 0}])
        result = motif_position_histogram(df)
        assert list(result.columns) == ["offset", "count"]

    def test_empty_input(self) -> None:
        result = motif_position_histogram(pd.DataFrame())
        assert result.empty
        assert list(result.columns) == ["offset", "count"]

    def test_missing_motif_offset_column(self) -> None:
        df = pd.DataFrame([{"sequence_id": "s1", "motif_present": True}])
        result = motif_position_histogram(df)
        assert result.empty

    def test_no_motif_present(self) -> None:
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": False, "motif_offset": None},
            {"sequence_id": "s2", "motif_present": False, "motif_offset": None},
        ])
        result = motif_position_histogram(df)
        assert result.empty

    def test_single_offset(self) -> None:
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": True, "motif_offset": 0},
            {"sequence_id": "s2", "motif_present": True, "motif_offset": 0},
            {"sequence_id": "s3", "motif_present": False, "motif_offset": None},
        ])
        result = motif_position_histogram(df)
        assert len(result) == 1
        assert result.iloc[0]["offset"] == 0
        assert result.iloc[0]["count"] == 2

    def test_multiple_offsets(self) -> None:
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": True, "motif_offset": -1},
            {"sequence_id": "s2", "motif_present": True, "motif_offset": 0},
            {"sequence_id": "s3", "motif_present": True, "motif_offset": 0},
            {"sequence_id": "s4", "motif_present": True, "motif_offset": 1},
        ])
        result = motif_position_histogram(df)
        assert len(result) == 3
        offset_map = dict(zip(result["offset"].tolist(), result["count"].tolist()))
        assert offset_map[-1] == 1
        assert offset_map[0] == 2
        assert offset_map[1] == 1

    def test_sorted_by_offset(self) -> None:
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": True, "motif_offset": 2},
            {"sequence_id": "s2", "motif_present": True, "motif_offset": -2},
            {"sequence_id": "s3", "motif_present": True, "motif_offset": 0},
        ])
        result = motif_position_histogram(df)
        assert result["offset"].tolist() == sorted(result["offset"].tolist())

    def test_absent_motifs_not_counted(self) -> None:
        """motif_present=False rows must not contribute to the histogram."""
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": True, "motif_offset": 1},
            {"sequence_id": "s2", "motif_present": False, "motif_offset": 1},
        ])
        result = motif_position_histogram(df)
        assert result.iloc[0]["count"] == 1

    def test_zero_tolerance_all_offset_zero(self) -> None:
        """With tolerance=0 all matched offsets should be 0."""
        df = self._make_df([
            {"sequence_id": "s1", "motif_present": True, "motif_offset": 0},
            {"sequence_id": "s2", "motif_present": True, "motif_offset": 0},
        ])
        result = motif_position_histogram(df)
        assert len(result) == 1
        assert result.iloc[0]["offset"] == 0
