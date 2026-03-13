"""Unit tests for coevo.blast.blast_parser."""

from __future__ import annotations

import io
import textwrap
from pathlib import Path

import pandas as pd
import pytest

from coevo.blast.blast_parser import (
    parse_blast_tabular,
    extract_taxids,
    filter_blast_hits,
    BLAST_COLUMNS,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_TSV = textwrap.dedent(
    """\
    query1\thit1\t98.5\t300\t1e-50\t500.0\t12345
    query1\thit2\t45.2\t200\t1e-10\t200.0\t67890;11111
    query2\thit3\t25.0\t100\t0.01\t80.0\tN/A
    query2\thit4\t80.0\t250\t1e-30\t350.0\t22222
    """
)


@pytest.fixture
def blast_tsv_file(tmp_path: Path) -> Path:
    """Write a small synthetic BLAST tabular file and return its path."""
    p = tmp_path / "blast_hits.tsv"
    p.write_text(SAMPLE_TSV)
    return p


@pytest.fixture
def blast_df(blast_tsv_file: Path) -> pd.DataFrame:
    """Parsed BLAST DataFrame from the sample file."""
    return parse_blast_tabular(blast_tsv_file)


# ---------------------------------------------------------------------------
# Tests for parse_blast_tabular
# ---------------------------------------------------------------------------


class TestParseBlastTabular:
    def test_returns_dataframe(self, blast_df: pd.DataFrame) -> None:
        assert isinstance(blast_df, pd.DataFrame)

    def test_column_names(self, blast_df: pd.DataFrame) -> None:
        assert list(blast_df.columns) == BLAST_COLUMNS

    def test_row_count(self, blast_df: pd.DataFrame) -> None:
        assert len(blast_df) == 4

    def test_pident_dtype(self, blast_df: pd.DataFrame) -> None:
        assert blast_df["pident"].dtype == float

    def test_evalue_dtype(self, blast_df: pd.DataFrame) -> None:
        assert blast_df["evalue"].dtype == float

    def test_length_dtype(self, blast_df: pd.DataFrame) -> None:
        assert blast_df["length"].dtype == int

    def test_values(self, blast_df: pd.DataFrame) -> None:
        first = blast_df.iloc[0]
        assert first["qseqid"] == "query1"
        assert first["sseqid"] == "hit1"
        assert first["pident"] == pytest.approx(98.5)
        assert first["staxids"] == "12345"

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(FileNotFoundError):
            parse_blast_tabular(tmp_path / "nonexistent.tsv")

    def test_skips_comment_lines(self, tmp_path: Path) -> None:
        content = "# comment\n" + SAMPLE_TSV
        p = tmp_path / "commented.tsv"
        p.write_text(content)
        df = parse_blast_tabular(p)
        assert len(df) == 4


# ---------------------------------------------------------------------------
# Tests for extract_taxids
# ---------------------------------------------------------------------------


class TestExtractTaxids:
    def test_returns_set_of_ints(self, blast_df: pd.DataFrame) -> None:
        taxids = extract_taxids(blast_df)
        assert isinstance(taxids, set)
        assert all(isinstance(t, int) for t in taxids)

    def test_single_taxid(self, blast_df: pd.DataFrame) -> None:
        taxids = extract_taxids(blast_df)
        assert 12345 in taxids

    def test_semicolon_separated_taxids(self, blast_df: pd.DataFrame) -> None:
        taxids = extract_taxids(blast_df)
        # hit2 has "67890;11111"
        assert 67890 in taxids
        assert 11111 in taxids

    def test_na_values_ignored(self, blast_df: pd.DataFrame) -> None:
        taxids = extract_taxids(blast_df)
        # "N/A" should not produce a taxid
        # All expected valid taxids: 12345, 67890, 11111, 22222
        assert len(taxids) == 4

    def test_unique_taxids(self, tmp_path: Path) -> None:
        content = "q\ts\t90.0\t100\t1e-5\t200.0\t999\nq\ts2\t88.0\t90\t1e-4\t150.0\t999\n"
        p = tmp_path / "dup.tsv"
        p.write_text(content)
        df = parse_blast_tabular(p)
        taxids = extract_taxids(df)
        assert taxids == {999}

    def test_empty_dataframe(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        taxids = extract_taxids(df)
        assert taxids == set()


# ---------------------------------------------------------------------------
# Tests for filter_blast_hits
# ---------------------------------------------------------------------------


class TestFilterBlastHits:
    def test_filter_by_identity(self, blast_df: pd.DataFrame) -> None:
        filtered = filter_blast_hits(blast_df, min_identity=50.0)
        # Only hits with pident >= 50: 98.5 and 80.0 => 2 rows
        assert len(filtered) == 2

    def test_filter_keeps_all_above_threshold(self, blast_df: pd.DataFrame) -> None:
        filtered = filter_blast_hits(blast_df, min_identity=0.0)
        assert len(filtered) == len(blast_df)

    def test_filter_with_coverage(self, blast_df: pd.DataFrame) -> None:
        # query_length=300, min_coverage=0.9 => length/300 >= 0.9 => length >= 270
        # Rows with length >= 270: 300 and 300 (first), 250 (no)... let's check:
        # hit1: length=300/300=1.0 >= 0.9 and pident=98.5 >= 0 -> yes
        # hit2: length=200/300=0.67 -> no
        # hit3: length=100/300=0.33 -> no
        # hit4: length=250/300=0.83 -> no
        filtered = filter_blast_hits(blast_df, min_identity=0.0, min_coverage=0.9, query_length=300)
        assert len(filtered) == 1
        assert filtered.iloc[0]["sseqid"] == "hit1"
