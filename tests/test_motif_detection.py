"""Unit tests for coevo.sequences.motif_detection."""

from __future__ import annotations

import textwrap
from pathlib import Path

import pandas as pd
import pytest

from coevo.sequences.motif_detection import (
    ALPHABETS,
    _get_alphabet,
    _parse_fragment,
    _match_fragment_at,
    _check_motif,
    _ungapped_to_aligned_index,
    detect_motif_in_alignment,
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
# Tests for ALPHABETS / _get_alphabet
# ---------------------------------------------------------------------------


class TestGetAlphabet:
    def test_rna_lowercase(self) -> None:
        assert _get_alphabet("rna") == frozenset("acgu")

    def test_rna_uppercase(self) -> None:
        assert _get_alphabet("RNA") == frozenset("acgu")

    def test_dna_lowercase(self) -> None:
        assert _get_alphabet("dna") == frozenset("acgt")

    def test_dna_uppercase(self) -> None:
        assert _get_alphabet("DNA") == frozenset("acgt")

    def test_protein_lowercase(self) -> None:
        assert _get_alphabet("protein") == frozenset("acdefghiklmnpqrstvwy")

    def test_protein_uppercase(self) -> None:
        assert _get_alphabet("Protein") == frozenset("acdefghiklmnpqrstvwy")

    def test_unknown_raises(self) -> None:
        with pytest.raises(ValueError, match="Unknown molecule type"):
            _get_alphabet("lipid")


# ---------------------------------------------------------------------------
# Tests for _parse_fragment
# ---------------------------------------------------------------------------


class TestParseFragment:
    """Tests for fragment pattern parsing."""

    RNA = ALPHABETS["rna"]
    PROTEIN = ALPHABETS["protein"]

    # -- plain characters --

    def test_single_plain_char(self) -> None:
        result = _parse_fragment("a", self.RNA)
        assert result == [frozenset({"a"})]

    def test_multiple_plain_chars(self) -> None:
        result = _parse_fragment("agc", self.RNA)
        assert result == [frozenset({"a"}), frozenset({"g"}), frozenset({"c"})]

    def test_plain_chars_case_normalised(self) -> None:
        result = _parse_fragment("A", self.RNA)
        assert result == [frozenset({"a"})]

    # -- character classes [ag] --

    def test_char_class(self) -> None:
        result = _parse_fragment("[ag]", self.RNA)
        assert result == [frozenset({"a", "g"})]

    def test_char_class_uppercase_normalised(self) -> None:
        result = _parse_fragment("[AG]", self.RNA)
        assert result == [frozenset({"a", "g"})]

    def test_char_class_mixed_with_plain(self) -> None:
        result = _parse_fragment("a[cg]u", self.RNA)
        assert result == [frozenset({"a"}), frozenset({"c", "g"}), frozenset({"u"})]

    def test_unclosed_bracket_raises(self) -> None:
        with pytest.raises(ValueError, match="Unclosed"):
            _parse_fragment("[ag", self.RNA)

    # -- negated character classes {a} --

    def test_negated_class(self) -> None:
        # {a} with RNA alphabet: all of acgu except a -> cgu
        result = _parse_fragment("{a}", self.RNA)
        assert result == [frozenset({"c", "g", "u"})]

    def test_negated_class_multi_char(self) -> None:
        # {ag} with RNA: acgu - {a,g} = {c,u}
        result = _parse_fragment("{ag}", self.RNA)
        assert result == [frozenset({"c", "u"})]

    def test_negated_class_protein(self) -> None:
        # {PH} with protein: protein_alphabet - {p, h}
        expected = ALPHABETS["protein"] - frozenset({"p", "h"})
        result = _parse_fragment("{PH}", self.PROTEIN)
        assert result == [expected]

    def test_unclosed_curly_raises(self) -> None:
        with pytest.raises(ValueError, match="Unclosed"):
            _parse_fragment("{ag", self.RNA)

    # -- wildcard x --

    def test_wildcard_x(self) -> None:
        result = _parse_fragment("x", self.RNA)
        assert result == [self.RNA]

    def test_wildcard_X_uppercase(self) -> None:
        result = _parse_fragment("X", self.RNA)
        assert result == [self.RNA]

    def test_wildcard_x_protein(self) -> None:
        result = _parse_fragment("x", self.PROTEIN)
        assert result == [self.PROTEIN]

    # -- combined --

    def test_combined_pattern(self) -> None:
        # "a[cg]{a}x" has 4 tokens
        result = _parse_fragment("a[cg]{a}x", self.RNA)
        assert len(result) == 4
        assert result[0] == frozenset({"a"})
        assert result[1] == frozenset({"c", "g"})
        assert result[2] == frozenset({"c", "g", "u"})  # rna - {a}
        assert result[3] == self.RNA

    def test_empty_fragment(self) -> None:
        result = _parse_fragment("", self.RNA)
        assert result == []


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
# Tests for _match_fragment_at
# ---------------------------------------------------------------------------


class TestMatchFragmentAt:
    RNA = ALPHABETS["rna"]

    def test_plain_match(self) -> None:
        char_sets = _parse_fragment("acg", self.RNA)
        assert _match_fragment_at("ACGU", 0, char_sets) is True

    def test_plain_no_match(self) -> None:
        char_sets = _parse_fragment("acg", self.RNA)
        assert _match_fragment_at("GCAU", 0, char_sets) is False

    def test_fragment_beyond_sequence(self) -> None:
        char_sets = _parse_fragment("acg", self.RNA)
        assert _match_fragment_at("AC", 0, char_sets) is False

    def test_gapped_sequence(self) -> None:
        # "A-CGU": ungapped 0=A, 1=C, 2=G, 3=U
        char_sets = _parse_fragment("cg", self.RNA)
        assert _match_fragment_at("A-CGU", 1, char_sets) is True

    def test_wildcard_matches_any(self) -> None:
        char_sets = _parse_fragment("x", self.RNA)
        for ch in "acgu":
            assert _match_fragment_at(ch, 0, char_sets) is True

    def test_char_class_match(self) -> None:
        char_sets = _parse_fragment("[ag]", self.RNA)
        assert _match_fragment_at("A", 0, char_sets) is True
        assert _match_fragment_at("G", 0, char_sets) is True
        assert _match_fragment_at("C", 0, char_sets) is False

    def test_negated_class_match(self) -> None:
        # {a} with RNA matches c, g, u but not a
        char_sets = _parse_fragment("{a}", self.RNA)
        assert _match_fragment_at("C", 0, char_sets) is True
        assert _match_fragment_at("A", 0, char_sets) is False


# ---------------------------------------------------------------------------
# Tests for _check_motif (unit-level)
# ---------------------------------------------------------------------------


class TestCheckMotif:
    RNA = ALPHABETS["rna"]

    def _frag(self, pattern: str) -> list[frozenset[str]]:
        return _parse_fragment(pattern, self.RNA)

    def test_motif_present(self) -> None:
        # positions 1,2,3 (1-based) -> ungapped indices 0,1,2 -> 'A','C','G'
        result = _check_motif(
            "seq1", "ACGU",
            [1, 2, 3],
            [self._frag("A"), self._frag("C"), self._frag("G")],
        )
        assert result["motif_present"] is True

    def test_motif_absent_wrong_residue(self) -> None:
        result = _check_motif(
            "seq2", "ACGU",
            [1, 2, 3],
            [self._frag("A"), self._frag("C"), self._frag("U")],
        )
        assert result["motif_present"] is False

    def test_case_insensitive(self) -> None:
        result = _check_motif(
            "seq1", "acgu",
            [1, 2, 3],
            [self._frag("A"), self._frag("C"), self._frag("G")],
        )
        assert result["motif_present"] is True

    def test_position_out_of_range(self) -> None:
        result = _check_motif("seq1", "ACG", [10], [self._frag("A")])
        assert result["motif_present"] is False

    def test_single_position(self) -> None:
        # position 4 (1-based) -> ungapped index 3 -> 'U'
        result = _check_motif("seq1", "ACGU", [4], [self._frag("U")])
        assert result["motif_present"] is True

    def test_sequence_id_preserved(self) -> None:
        result = _check_motif("my_seq", "AAAA", [1], [self._frag("A")])
        assert result["sequence_id"] == "my_seq"

    def test_dash_shifts_index(self) -> None:
        # "-CGU": position 1 (1-based) -> ungapped index 0 -> aligned index 1 -> 'C'
        result = _check_motif("gseq", "-CGU", [1], [self._frag("C")])
        assert result["motif_present"] is True

    def test_dash_shifts_index_mismatch(self) -> None:
        result = _check_motif("gseq", "-CGU", [1], [self._frag("A")])
        assert result["motif_present"] is False

    def test_multi_char_fragment_match(self) -> None:
        # "acg" fragment at position 1 checks pos 1,2,3 -> 'A','C','G'
        result = _check_motif("seq1", "ACGU", [1], [self._frag("acg")])
        assert result["motif_present"] is True

    def test_multi_char_fragment_no_match(self) -> None:
        result = _check_motif("seq1", "ACGU", [1], [self._frag("acu")])
        assert result["motif_present"] is False

    def test_char_class_match(self) -> None:
        # [AG] at position 1: 'A' is in {a,g}
        result = _check_motif("seq1", "ACGU", [1], [self._frag("[AG]")])
        assert result["motif_present"] is True

    def test_negated_class_match(self) -> None:
        # {C} at position 1: 'A' is in RNA-alphabet minus {c} = {a,g,u}
        result = _check_motif("seq1", "ACGU", [1], [self._frag("{C}")])
        assert result["motif_present"] is True

    def test_negated_class_no_match(self) -> None:
        # {A} at position 1: 'A' is NOT in RNA-alphabet minus {a} = {c,g,u}
        result = _check_motif("seq1", "ACGU", [1], [self._frag("{A}")])
        assert result["motif_present"] is False

    def test_wildcard_matches(self) -> None:
        result = _check_motif("seq1", "ACGU", [1], [self._frag("x")])
        assert result["motif_present"] is True

    # -- tolerance --

    def test_tolerance_zero_exact(self) -> None:
        # "g" at position 5 (1-based); sequence "AACTTG" has 'g' at pos 6
        result = _check_motif("s", "AACTTG", [5], [self._frag("g")], tolerance=0)
        assert result["motif_present"] is False

    def test_tolerance_one_allows_shift(self) -> None:
        # "g" at position 5 ± 1 = [4,6]; "AACTTG" has 'g' at ungapped pos 6
        result = _check_motif("s", "AACTTG", [5], [self._frag("g")], tolerance=1)
        assert result["motif_present"] is True

    def test_tolerance_one_left_shift(self) -> None:
        # "g" at position 5 ± 1 = [4,6]; "AACG" has 'g' at ungapped pos 4
        result = _check_motif("s", "AACG", [5], [self._frag("g")], tolerance=1)
        assert result["motif_present"] is True

    def test_tolerance_does_not_go_below_zero(self) -> None:
        # tolerance should not cause negative ungapped index
        result = _check_motif("s", "AG", [1], [self._frag("a")], tolerance=2)
        assert result["motif_present"] is True

    # -- motif_offset --

    def test_motif_offset_zero_on_exact_match(self) -> None:
        result = _check_motif("s", "ACGU", [1], [self._frag("A")], tolerance=0)
        assert result["motif_present"] is True
        assert result["motif_offset"] == 0

    def test_motif_offset_none_when_absent(self) -> None:
        result = _check_motif("s", "ACGU", [1], [self._frag("G")], tolerance=0)
        assert result["motif_present"] is False
        assert result["motif_offset"] is None

    def test_motif_offset_positive_shift(self) -> None:
        # "G" at position 5 (1-based) ± 1; "AACTTG" has 'G' at ungapped pos 6 -> delta=+1
        result = _check_motif("s", "AACTTG", [5], [self._frag("g")], tolerance=1)
        assert result["motif_present"] is True
        assert result["motif_offset"] == 1

    def test_motif_offset_negative_shift(self) -> None:
        # "G" at position 5 ± 1; "AACG" has 'G' at ungapped pos 4 -> delta=-1
        result = _check_motif("s", "AACG", [5], [self._frag("g")], tolerance=1)
        assert result["motif_present"] is True
        assert result["motif_offset"] == -1

    def test_motif_offset_key_present_in_result(self) -> None:
        result = _check_motif("s", "ACGU", [1], [self._frag("A")])
        assert "motif_offset" in result




class TestDetectMotifInAlignment:
    def test_returns_dataframe(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            fragments=["A", "C", "G"],
        )
        assert isinstance(df, pd.DataFrame)

    def test_column_names(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["A"],
        )
        assert "sequence_id" in df.columns
        assert "motif_present" in df.columns
        assert "motif_offset" in df.columns

    def test_row_count(self, alignment_file: Path) -> None:
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["A"],
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
            fragments=["A", "C", "G"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] is True or results["seq1"] == True  # noqa: E712
        assert results["seq2"] is True or results["seq2"] == True  # noqa: E712
        assert results["seq3"] is False or results["seq3"] == False  # noqa: E712
        assert results["seq4"] is False or results["seq4"] == False  # noqa: E712

    def test_mismatched_positions_fragments(self, alignment_file: Path) -> None:
        with pytest.raises(ValueError, match="same length"):
            detect_motif_in_alignment(
                alignment_file=str(alignment_file),
                positions=[1, 2],
                fragments=["A"],
            )

    def test_file_not_found(self, tmp_path: Path) -> None:
        with pytest.raises(Exception):
            detect_motif_in_alignment(
                alignment_file=str(tmp_path / "missing.fasta"),
                positions=[1],
                fragments=["A"],
            )

    def test_no_motif_anywhere(self, alignment_file: Path) -> None:
        # Use a residue that doesn't appear at position 1 in any sequence
        # seq1 ungapped[0]=A, seq2 ungapped[0]=A, seq3 ungapped[0]=G, seq4 ungapped[0]=N
        # => "Z" matches none
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["Z"],
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
            fragments=["C"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq2"] == True  # noqa: E712
        assert results["seq3"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712

    def test_gapped_alignment_dash_correction(self, gapped_alignment_file: Path) -> None:
        # gseq1: "A-CGU" -> ungapped: ACGU
        #   position 2 (1-based) -> ungapped index 1 -> aligned index 2 -> 'C' ✓
        # gseq2: "--CGU" -> ungapped: CGU
        #   position 2 (1-based) -> ungapped index 1 -> aligned index 3 -> 'G' ✗
        df = detect_motif_in_alignment(
            alignment_file=str(gapped_alignment_file),
            positions=[2],
            fragments=["C"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["gseq1"] == True  # noqa: E712
        # gseq2 ungapped position 2 is 'G', not 'C'
        assert results["gseq2"] == False  # noqa: E712

    def test_multi_char_fragment_in_alignment(self, alignment_file: Path) -> None:
        # position 1 (1-based): seq1='A', seq3='G' -> fragment "[AG]" matches both
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["[AG]"],
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq3"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712

    def test_fragment_equivalence(self, alignment_file: Path) -> None:
        """positions:[1,2,3] fragments:["a","c","g"] == positions:[1] fragments:["acg"]."""
        df_verbose = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            fragments=["A", "C", "G"],
        )
        df_compact = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["ACG"],
        )
        pd.testing.assert_frame_equal(
            df_verbose.set_index("sequence_id").sort_index(),
            df_compact.set_index("sequence_id").sort_index(),
        )

    def test_molecule_type_protein(self, alignment_file: Path) -> None:
        # wildcard 'x' with protein should match any protein alphabet char
        # seq1 pos 1 = 'A', which is in protein alphabet -> x matches
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["x"],
            molecule_type="protein",
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712

    def test_unknown_molecule_type_raises(self, alignment_file: Path) -> None:
        with pytest.raises(ValueError, match="Unknown molecule type"):
            detect_motif_in_alignment(
                alignment_file=str(alignment_file),
                positions=[1],
                fragments=["A"],
                molecule_type="lipid",
            )

    def test_tolerance_in_alignment(self, alignment_file: Path) -> None:
        """Fragment at a shifted position still matches when tolerance allows it.

        seq1: ACGUACGUACGU
        position 4 (1-based) is 'U' (ungapped index 3).
        With tolerance=1, position 5 (1-based, ungapped index 4) = 'A'.
        Check that 'U' at declared position 4 also matches with tolerance=0,
        and that 'U' at declared position 5 matches with tolerance=1.
        """
        # exact: 'U' is at position 4 -> match for seq1
        df_exact = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[4],
            fragments=["U"],
            tolerance=0,
        )
        # tolerance=1: 'U' declared at position 5, but 'U' is actually at pos 4
        df_tol = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[5],
            fragments=["U"],
            tolerance=1,
        )
        exact_results = df_exact.set_index("sequence_id")["motif_present"]
        tol_results = df_tol.set_index("sequence_id")["motif_present"]
        assert exact_results["seq1"] == True  # noqa: E712
        assert tol_results["seq1"] == True  # noqa: E712

    def test_negated_class_in_alignment(self, alignment_file: Path) -> None:
        # {A} with RNA at position 1: accepts everything EXCEPT 'a'
        # seq1 pos1='A' -> no match; seq3 pos1='G' -> match
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["{A}"],
            molecule_type="rna",
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == False  # noqa: E712
        assert results["seq3"] == True  # noqa: E712

    def test_wildcard_in_alignment(self, alignment_file: Path) -> None:
        # 'x' (wildcard) at pos 1 matches any RNA character.
        # seq4 has 'N' which is NOT in the RNA alphabet, so should NOT match.
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1],
            fragments=["x"],
            molecule_type="rna",
        )
        results = df.set_index("sequence_id")["motif_present"]
        assert results["seq1"] == True  # noqa: E712
        assert results["seq4"] == False  # noqa: E712


    def test_motif_offset_zero_exact_position(self, alignment_file: Path) -> None:
        """Motif found at exact position -> offset 0 for matched sequences."""
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            fragments=["A", "C", "G"],
            tolerance=0,
        )
        present = df[df["motif_present"] == True]  # noqa: E712
        assert (present["motif_offset"] == 0).all()

    def test_motif_offset_none_for_absent(self, alignment_file: Path) -> None:
        """Sequences without the motif have motif_offset=None."""
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[1, 2, 3],
            fragments=["A", "C", "G"],
            tolerance=0,
        )
        absent = df[df["motif_present"] == False]  # noqa: E712
        assert absent["motif_offset"].isna().all()

    def test_motif_offset_with_tolerance(self, alignment_file: Path) -> None:
        """With tolerance, offset can be non-zero for shifted matches."""
        # position 4 (1-based) ± 1 = [3,5]; seq1 has 'U' at ungapped pos 4 (exact)
        # so offset should be 0 for exact match
        df = detect_motif_in_alignment(
            alignment_file=str(alignment_file),
            positions=[4],
            fragments=["U"],
            molecule_type="rna",
            tolerance=1,
        )
        seq1_row = df[df["sequence_id"] == "seq1"]
        assert seq1_row.iloc[0]["motif_present"] is True or seq1_row.iloc[0]["motif_present"] == True  # noqa: E712
        assert seq1_row.iloc[0]["motif_offset"] == 0
