"""Unit tests for coevo.sequences.motif_detection."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from coevo.sequences.motif_detection import (
    detect_motif_in_alignment,
    _check_motif,
    _ungapped_to_aligned_index,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

ALIGNMENT_FASTA = textwrap.dedent(
    """\
    >seq1
    ACGUACGUACGU
    >seq2
    ACGUACGUAAAU
    >seq3
    GCAUACGUACGU
    >seq4
    NNNNNNNNNNNN
    """
)

# Alignment fixture with gap characters to test dash-shift correction.
GAPPED_FASTA = textwrap.dedent(
    """\
    >gseq1
    A-CGU
    >gseq2
    --CGU
    """
)


@pytest.fixture
def alignment_file(tmp_path: Path) -> Path:
    """Write a tiny synthetic alignment FASTA and return its path."""
    p = tmp_path / "alignment.fasta"
    p.write_text(ALIGNMENT_FASTA)
    return p


@pytest.fixture
def gapped_alignment_file(tmp_path: Path) -> Path:
    """Write a synthetic alignment with gap characters and return its path."""
    p = tmp_path / "gapped_alignment.fasta"
    p.write_text(GAPPED_FASTA)
    return p


# ---------------------------------------------------------------------------
# Tests for _ungapped_to_aligned_index
# ---------------------------------------------------------------------------


class TestUngappedToAlignedIndex:
    def test_no_gaps(self) -> None:
        # No gaps: aligned index equals ungapped index
        assert _ungapped_to_aligned_index("ACGU", 0) == 0
        assert _ungapped_to_aligned_index("ACGU", 3) == 3

    def test_single_leading_gap(self) -> None:
        # "-CGU": ungapped indices 0,1,2 map to aligned 1,2,3
        assert _ungapped_to_aligned_index("-CGU", 0) == 1
        assert _ungapped_to_aligned_index("-CGU", 1) == 2

    def test_internal_gap(self) -> None:
        # "A-GU": ungapped 0->0, 1->2, 2->3
        assert _ungapped_to_aligned_index("A-GU", 0) == 0
        assert _ungapped_to_aligned_index("A-GU", 1) == 2
        assert _ungapped_to_aligned_index("A-GU", 2) == 3

    def test_out_of_range_returns_none(self) -> None:
        assert _ungapped_to_aligned_index("ACG", 10) is None

    def test_all_gaps_returns_none(self) -> None:
        assert _ungapped_to_aligned_index("---", 0) is None


# ---------------------------------------------------------------------------
# Tests for _check_motif (unit-level)
# ---------------------------------------------------------------------------


class TestCheckMotif:
    def test_motif_present(self) -> None:
        # positions 1,2,3 (1-based) -> ungapped indices 0,1,2 -> 'A','C','G'
        result = _check_motif("seq1", "ACGU", [1, 2, 3], ["A", "C", "G"])
        assert result["motif_present"] is True

    def test_motif_absent_wrong_residue(self) -> None:
        result = _check_motif("seq2", "ACGU", [1, 2, 3], ["A", "C", "U"])
        assert result["motif_present"] is False

    def test_case_insensitive(self) -> None:
        result = _check_motif("seq1", "acgu", [1, 2, 3], ["A", "C", "G"])
        assert result["motif_present"] is True

    def test_position_out_of_range(self) -> None:
        result = _check_motif("seq1", "ACG", [10], ["A"])
        assert result["motif_present"] is False

    def test_single_position(self) -> None:
        # position 4 (1-based) -> ungapped index 3 -> 'U'
        result = _check_motif("seq1", "ACGU", [4], ["U"])
        assert result["motif_present"] is True

    def test_sequence_id_preserved(self) -> None:
        result = _check_motif("my_seq", "AAAA", [1], ["A"])
        assert result["sequence_id"] == "my_seq"

    def test_dash_shifts_index(self) -> None:
        # "-CGU": position 1 (1-based) -> ungapped index 0 -> aligned index 1 -> 'C'
        result = _check_motif("gseq", "-CGU", [1], ["C"])
        assert result["motif_present"] is True

    def test_dash_shifts_index_mismatch(self) -> None:
        # "-CGU": position 1 (1-based) -> ungapped index 0 -> aligned 'C', not 'A'
        result = _check_motif("gseq", "-CGU", [1], ["A"])
        assert result["motif_present"] is False

    def test_multi_char_residue_match(self) -> None:
        # residue "ACG" means accept 'A', 'C', or 'G'; seq char at pos 1 is 'A'
        result = _check_motif("seq1", "ACGU", [1], ["ACG"])
        assert result["motif_present"] is True

    def test_multi_char_residue_no_match(self) -> None:
        # residue "CG" does not include 'A'; seq char at pos 1 is 'A'
        result = _check_motif("seq1", "ACGU", [1], ["CG"])
        assert result["motif_present"] is False


# ---------------------------------------------------------------------------
# Tests for detect_motif_in_alignment
# ---------------------------------------------------------------------------


class TestDetectMotifInAlignment:
    def test_returns_dataframe(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            residues=["A", "C", "G"],
        )
        assert isinstance(df, pd.DataFrame)

    def test_column_names(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            residues=["A"],
        )
        assert "sequence_id" in df.columns
        assert "motif_present" in df.columns

    def test_row_count(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            residues=["A"],
        )
        # 4 sequences in the alignment
        assert len(df) == 4

    def test_motif_present_correct_sequences(self, alignment_file: Path) -> None:
        # positions [1,2,3] (1-based) -> ungapped indices [0,1,2] -> ["A","C","G"]:
        # seq1: ACG -> present
        # seq2: ACG -> present
        # seq3: GCA -> absent
        # seq4: NNN -> absent
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            residues=["A", "C", "G"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] is True or results["seq1"] == True  # noqa: E712
        assert results["seq2"] is True or results["seq2"] == True  # noqa: E712
        assert results["seq3"] is False or results["seq3"] == False  # noqa: E712
        assert results["seq4"] is False or results["seq4"] == False  # noqa: E712

    def test_mismatched_positions_residues(self, alignment_file: Path) -> None:
        with pytest.raises(ValueError, match="same length"):
            detect_motif_in_alignment(
                alignment_file=str(alignment_file),
                positions=[1, 2],
                residues=["A"],
            )

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(Exception):
            detect_motif_in_alignment(
                alignment_file=str(tmp_path / "missing.fasta"),
                positions=[1],
                residues=["A"],
            )

    def test_no_motif_anywhere(self, alignment_file: Path) -> None:
        # Use a residue that doesn't appear at position 1 in any sequence
        # seq1 ungapped[0]=A, seq2 ungapped[0]=A, seq3 ungapped[0]=G, seq4 ungapped[0]=N
        # => "Z" matches none
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            residues=["Z"],
        )
        assert df["motif_present"].sum() == 0

    def test_all_motif_present(self, alignment_file: Path) -> None:
        # Position 6 (1-based) -> ungapped index 5 -> 'C' for seq1/seq2/seq3
        # seq1: ACGUACGUACGU  ungapped[5]=C
        # seq2: ACGUACGUAAAU  ungapped[5]=C
        # seq3: GCAUACGUACGU  ungapped[5]=C
        # seq4: NNNNNNNNNNNN  ungapped[5]=N
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[6],
            residues=["C"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq2"] == True  # noqa: E712
        assert results["seq3"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712

    def test_gapped_alignment_dash_correction(self, gapped_alignment_file: Path) -> None:
        # gseq1: "A-CGU" -> ungapped: ACGU
        #   position 2 (1-based) -> ungapped index 1 -> aligned index 2 -> 'C'
        # gseq2: "--CGU" -> ungapped: CGU
        #   position 1 (1-based) -> ungapped index 0 -> aligned index 2 -> 'C'
        df = detect_motif_in_alignment(
            alignment_file=str(gapped_alignment_file),
            positions=[2],
            residues=["C"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["gseq1"] == True  # noqa: E712
        # gseq2 ungapped position 2 is 'G', not 'C'
        assert results["gseq2"] == False  # noqa: E712

    def test_multi_char_residue_in_alignment(self, alignment_file: Path) -> None:
        # position 1 (1-based): seq1='A', seq3='G' -> residue "AG" matches both
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            residues=["AG"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq3"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712

