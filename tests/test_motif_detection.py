"""Unit tests for coevo.sequences.motif_detection."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from coevo.sequences.motif_detection import detect_motif_in_alignment, _check_motif


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


@pytest.fixture
def alignment_file(tmp_path: Path) -> Path:
    """Write a tiny synthetic alignment FASTA and return its path."""
    p = tmp_path / "alignment.fasta"
    p.write_text(ALIGNMENT_FASTA)
    return p


# ---------------------------------------------------------------------------
# Tests for _check_motif (unit-level)
# ---------------------------------------------------------------------------


class TestCheckMotif:
    def test_motif_present(self) -> None:
        result = _check_motif("seq1", "ACGU", [0, 1, 2], ["A", "C", "G"])
        assert result["motif_present"] is True

    def test_motif_absent_wrong_residue(self) -> None:
        result = _check_motif("seq2", "ACGU", [0, 1, 2], ["A", "C", "U"])
        assert result["motif_present"] is False

    def test_case_insensitive(self) -> None:
        result = _check_motif("seq1", "acgu", [0, 1, 2], ["A", "C", "G"])
        assert result["motif_present"] is True

    def test_position_out_of_range(self) -> None:
        result = _check_motif("seq1", "ACG", [10], ["A"])
        assert result["motif_present"] is False

    def test_single_position(self) -> None:
        result = _check_motif("seq1", "ACGU", [3], ["U"])
        assert result["motif_present"] is True

    def test_sequence_id_preserved(self) -> None:
        result = _check_motif("my_seq", "AAAA", [0], ["A"])
        assert result["sequence_id"] == "my_seq"


# ---------------------------------------------------------------------------
# Tests for detect_motif_in_alignment
# ---------------------------------------------------------------------------


class TestDetectMotifInAlignment:
    def test_returns_dataframe(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[0, 1, 2],
            residues=["A", "C", "G"],
        )
        assert isinstance(df, pd.DataFrame)

    def test_column_names(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[0],
            residues=["A"],
        )
        assert "sequence_id" in df.columns
        assert "motif_present" in df.columns

    def test_row_count(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[0],
            residues=["A"],
        )
        # 4 sequences in the alignment
        assert len(df) == 4

    def test_motif_present_correct_sequences(self, alignment_file: Path) -> None:
        # positions [0,1,2] with residues ["A","C","G"]:
        # seq1: ACG -> present
        # seq2: ACG -> present
        # seq3: GCA -> absent
        # seq4: NNN -> absent
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[0, 1, 2],
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
                positions=[0, 1],
                residues=["A"],
            )

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(Exception):
            detect_motif_in_alignment(
                alignment_file=str(tmp_path / "missing.fasta"),
                positions=[0],
                residues=["A"],
            )

    def test_no_motif_anywhere(self, alignment_file: Path) -> None:
        # Use a residue that doesn't appear at position 0 in any sequence
        # seq1[0]=A, seq2[0]=A, seq3[0]=G, seq4[0]=N  => "Z" matches none
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[0],
            residues=["Z"],
        )
        assert df["motif_present"].sum() == 0

    def test_all_motif_present(self, alignment_file: Path) -> None:
        # Position 4 is 'A' in seq1=ACGUACGUACGU (pos4='A'), seq2 (pos4='A'),
        # seq3 (pos4='A'), seq4 (pos4='N') -- seq4 won't match 'A'
        # Use position that all share: let's pick pos 4 = 'A' for first 3
        # Actually let's just use position 5 = 'C' for seq1/seq2/seq3 to be sure
        # seq1: ACGUACGUACGU  pos5=C
        # seq2: ACGUACGUAAAU  pos5=C
        # seq3: GCAUACGUACGU  pos5=C
        # seq4: NNNNNNNNNNNN  pos5=N
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[5],
            residues=["C"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq2"] == True  # noqa: E712
        assert results["seq3"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712
