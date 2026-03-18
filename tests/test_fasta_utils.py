"""Unit tests for coevo.sequences.fasta_utils.build_alignment_fasta."""

from __future__ import annotations

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from coevo.blast.blast_parser import BLAST_COLUMNS
from coevo.sequences.fasta_utils import build_alignment_fasta

# ---------------------------------------------------------------------------
# Helpers & fixtures
# ---------------------------------------------------------------------------

REFERENCE_RECORD = SeqRecord(Seq("ACGTACGT"), id="ref1", description="reference seq")


def _make_blast_df(rows: list[tuple]) -> pd.DataFrame:
    """Build a minimal BLAST DataFrame from tuples matching BLAST_COLUMNS."""
    return pd.DataFrame(rows, columns=BLAST_COLUMNS)


# ---------------------------------------------------------------------------
# Tests for build_alignment_fasta
# ---------------------------------------------------------------------------


class TestBuildAlignmentFasta:
    def test_returns_tuple_of_records_and_dataframe(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        result = build_alignment_fasta([REFERENCE_RECORD], df)
        records, metadata = result
        assert isinstance(records, list)
        assert isinstance(metadata, pd.DataFrame)

    def test_reference_is_first_entry(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert records[0].id == "ref1"

    def test_reference_sequence_preserved(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert str(records[0].seq) == str(REFERENCE_RECORD.seq)

    def test_empty_blast_returns_only_reference(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        records, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 1
        assert metadata.empty

    def test_unique_blast_sequences_appended(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "67890", "CCCCAAAA"),
        ]
        df = _make_blast_df(rows)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 3  # ref + 2 blast sequences

    def test_duplicate_blast_sequences_deduplicated(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "67890", "TTTTGGGG"),  # same seq
        ]
        df = _make_blast_df(rows)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 2  # ref + 1 unique blast sequence

    def test_blast_seq_identical_to_reference_not_duplicated_in_fasta(self) -> None:
        """A BLAST sequence identical to the reference must not be appended as a new record."""
        rows = [
            ("q", "hit1", 100.0, 8, 1e-50, 500.0, "12345", "ACGTACGT"),  # same as ref
        ]
        df = _make_blast_df(rows)
        records, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 1  # only reference, no duplicate
        # Taxid is linked to index 0 (the reference)
        assert len(metadata) == 1
        assert metadata.iloc[0]["seq_index"] == 0
        assert metadata.iloc[0]["taxid"] == 12345

    def test_blast_seq_with_gaps_identical_to_reference_not_duplicated_in_fasta(self) -> None:
        """A gapped BLAST sequence matching the reference after gap removal is not duplicated."""
        rows = [
            ("q", "hit1", 100.0, 8, 1e-50, 500.0, "12345", "ACGT-ACGT"),  # same as ref with gap
        ]
        df = _make_blast_df(rows)
        records, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 1  # only reference
        # Taxid is still linked to index 0
        assert metadata.iloc[0]["seq_index"] == 0
        assert metadata.iloc[0]["taxid"] == 12345

    def test_blast_seq_case_insensitive_dedup(self) -> None:
        """Sequences differing only in case count as duplicates."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "ttttgggg"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "67890", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 2  # ref + 1 unique

    def test_gap_characters_stripped_for_dedup(self) -> None:
        """Gaps (dashes) in sseq are stripped before deduplication."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTT-GGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "67890", "TTTTGGGG"),  # same seq w/o gap
        ]
        df = _make_blast_df(rows)
        records, _ = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(records) == 2  # ref + 1 unique

    # ------------------------------------------------------------------
    # Metadata tests
    # ------------------------------------------------------------------

    def test_metadata_has_correct_columns(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert list(metadata.columns) == ["seq_index", "taxid"]

    def test_metadata_seq_index_starts_at_one(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert metadata.iloc[0]["seq_index"] == 1

    def test_metadata_taxid_correct(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert metadata.iloc[0]["taxid"] == 12345

    def test_metadata_no_ref_index_when_blast_seq_differs(self) -> None:
        """Index 0 does not appear in metadata when no BLAST hit matches the reference."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert 0 not in metadata["seq_index"].values

    def test_metadata_multi_taxid_produces_multiple_rows(self) -> None:
        """Semicolon-separated taxids produce one metadata row each."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "12345;67890", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(metadata) == 2
        assert set(metadata["taxid"]) == {12345, 67890}
        assert all(metadata["seq_index"] == 1)

    def test_metadata_na_taxid_skipped(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "N/A", "TTTTGGGG"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert metadata.empty

    def test_metadata_two_sequences_correct_indices(self) -> None:
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "11111", "TTTTGGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "22222", "CCCCAAAA"),
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert set(metadata["seq_index"]) == {1, 2}
        idx1_taxids = set(metadata[metadata["seq_index"] == 1]["taxid"])
        idx2_taxids = set(metadata[metadata["seq_index"] == 2]["taxid"])
        assert idx1_taxids == {11111}
        assert idx2_taxids == {22222}

    def test_metadata_ref_identical_seq_multi_taxid_linked_to_index_zero(self) -> None:
        """All taxids from a ref-identical BLAST hit are linked to seq_index=0."""
        rows = [
            ("q", "hit1", 100.0, 8, 1e-50, 500.0, "11111;22222", "ACGTACGT"),  # same as ref
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert set(metadata["seq_index"]) == {0}
        assert set(metadata["taxid"]) == {11111, 22222}

    def test_metadata_dedup_seq_collects_taxids_from_all_occurrences(self) -> None:
        """A duplicate sequence contributes its taxids to the original seq_index."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "11111", "TTTTGGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "22222", "TTTTGGGG"),  # dup seq, different taxid
        ]
        df = _make_blast_df(rows)
        records, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        # Only one unique sequence stored (ref + 1 blast)
        assert len(records) == 2
        # Both taxids are linked to the same seq_index=1
        assert len(metadata) == 2
        assert set(metadata["taxid"]) == {11111, 22222}
        assert all(metadata["seq_index"] == 1)

    def test_metadata_dedup_taxids_deduplicated(self) -> None:
        """Identical (seq_index, taxid) pairs from repeated hits are deduplicated."""
        rows = [
            ("q", "hit1", 95.0, 8, 1e-10, 200.0, "11111", "TTTTGGGG"),
            ("q", "hit2", 90.0, 8, 1e-9, 180.0, "11111", "TTTTGGGG"),  # same seq, same taxid
        ]
        df = _make_blast_df(rows)
        _, metadata = build_alignment_fasta([REFERENCE_RECORD], df)
        assert len(metadata) == 1
        assert metadata.iloc[0]["taxid"] == 11111
        assert metadata.iloc[0]["seq_index"] == 1

    # ------------------------------------------------------------------
    # Error handling
    # ------------------------------------------------------------------

    def test_empty_reference_records_raises(self) -> None:
        df = pd.DataFrame(columns=BLAST_COLUMNS)
        with pytest.raises(ValueError, match="empty"):
            build_alignment_fasta([], df)
