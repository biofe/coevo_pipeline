"""Unit tests for coevo.analysis.cooccurrence and coevo.analysis.statistics."""

from __future__ import annotations

import pytest

from coevo.analysis.cooccurrence import compute_cooccurrence
from coevo.analysis.statistics import contingency_table, fisher_exact_test


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

PROTEIN_TAXA = {1, 2, 3, 4, 5}       # 5 organisms with protein
RNA_TAXA = {3, 4, 5, 6, 7}            # 5 organisms with motif
# Intersection: {3, 4, 5}  (3 organisms)
# Union: {1..7}  (7 organisms)
# Jaccard: 3/7 ≈ 0.4286


# ---------------------------------------------------------------------------
# Tests for compute_cooccurrence
# ---------------------------------------------------------------------------


class TestComputeCooccurrence:
    def test_returns_dict(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert isinstance(result, dict)

    def test_intersection(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["intersection"] == 3

    def test_union(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["union"] == 7

    def test_jaccard_index(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["jaccard_index"] == pytest.approx(3 / 7)

    def test_protein_total(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["protein_total"] == 5

    def test_rna_total(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["rna_total"] == 5

    def test_protein_only(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["protein_only"] == 2  # {1, 2}

    def test_rna_only(self) -> None:
        result = compute_cooccurrence(PROTEIN_TAXA, RNA_TAXA)
        assert result["rna_only"] == 2  # {6, 7}

    def test_empty_sets(self) -> None:
        result = compute_cooccurrence(set(), set())
        assert result["intersection"] == 0
        assert result["union"] == 0
        assert result["jaccard_index"] == 0.0

    def test_identical_sets(self) -> None:
        taxa = {1, 2, 3}
        result = compute_cooccurrence(taxa, taxa)
        assert result["jaccard_index"] == pytest.approx(1.0)
        assert result["intersection"] == 3
        assert result["union"] == 3

    def test_disjoint_sets(self) -> None:
        result = compute_cooccurrence({1, 2}, {3, 4})
        assert result["intersection"] == 0
        assert result["jaccard_index"] == pytest.approx(0.0)

    def test_subset(self) -> None:
        result = compute_cooccurrence({1, 2, 3}, {2, 3})
        assert result["intersection"] == 2
        assert result["union"] == 3
        assert result["jaccard_index"] == pytest.approx(2 / 3)


# ---------------------------------------------------------------------------
# Tests for contingency_table
# ---------------------------------------------------------------------------


class TestContingencyTable:
    def test_returns_dict(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA)
        assert isinstance(table, dict)

    def test_both(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA)
        assert table["both"] == 3

    def test_protein_only(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA)
        assert table["protein_only"] == 2

    def test_motif_only(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA)
        assert table["motif_only"] == 2

    def test_neither_default_universe(self) -> None:
        # Default universe = union, so neither == 0
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA)
        assert table["neither"] == 0

    def test_neither_with_all_taxa(self) -> None:
        all_taxa = set(range(1, 11))  # 10 organisms total
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA, all_taxa=all_taxa)
        # universe has 10; protein∪rna = 7; neither = 10 - 7 = 3
        assert table["neither"] == 3

    def test_empty_sets(self) -> None:
        table = contingency_table(set(), set())
        assert table["both"] == 0
        assert table["protein_only"] == 0
        assert table["motif_only"] == 0
        assert table["neither"] == 0


# ---------------------------------------------------------------------------
# Tests for fisher_exact_test
# ---------------------------------------------------------------------------


class TestFisherExactTest:
    def test_returns_dict(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA, all_taxa=set(range(1, 11)))
        result = fisher_exact_test(table)
        assert isinstance(result, dict)

    def test_has_odds_ratio(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA, all_taxa=set(range(1, 11)))
        result = fisher_exact_test(table)
        assert "odds_ratio" in result

    def test_has_p_value(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA, all_taxa=set(range(1, 11)))
        result = fisher_exact_test(table)
        assert "p_value" in result

    def test_p_value_range(self) -> None:
        table = contingency_table(PROTEIN_TAXA, RNA_TAXA, all_taxa=set(range(1, 11)))
        result = fisher_exact_test(table)
        assert 0.0 <= result["p_value"] <= 1.0

    def test_perfect_association(self) -> None:
        # All organisms with protein also have motif, and vice versa
        taxa = {1, 2, 3}
        all_taxa = {1, 2, 3, 4, 5}
        table = contingency_table(taxa, taxa, all_taxa=all_taxa)
        result = fisher_exact_test(table)
        # Perfect separation should give p-value <= 0.5 (small sample)
        assert result["p_value"] <= 0.5

    def test_no_association(self) -> None:
        # Balanced: 2 in both, 2 protein-only, 2 rna-only, 4 neither
        protein = {1, 2, 3, 4}
        rna = {3, 4, 5, 6}
        all_taxa = set(range(1, 11))
        table = contingency_table(protein, rna, all_taxa=all_taxa)
        result = fisher_exact_test(table)
        assert isinstance(result["odds_ratio"], float)
