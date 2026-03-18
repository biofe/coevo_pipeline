"""Motif detection in multiple sequence alignments."""

from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from loguru import logger


# ---------------------------------------------------------------------------
# Molecule alphabets
# ---------------------------------------------------------------------------

#: Character alphabets for supported molecule types (lower-case).
ALPHABETS: dict[str, frozenset[str]] = {
    "protein": frozenset("acdefghiklmnpqrstvwy"),
    "dna": frozenset("acgt"),
    "rna": frozenset("acgu"),
}

_MOLECULE_ALIASES: dict[str, str] = {
    "protein": "protein",
    "Protein": "protein",
    "dna": "dna",
    "DNA": "dna",
    "rna": "rna",
    "RNA": "rna",
}


def _get_alphabet(molecule_type: str) -> frozenset[str]:
    """Return the lower-case character alphabet for *molecule_type*.

    Parameters
    ----------
    molecule_type:
        One of ``"protein"`` / ``"Protein"``, ``"dna"`` / ``"DNA"``,
        ``"rna"`` / ``"RNA"`` (case-insensitive).

    Raises
    ------
    ValueError
        If *molecule_type* is not recognised.
    """
    key = _MOLECULE_ALIASES.get(molecule_type)
    if key is None:
        valid = sorted({k.lower() for k in _MOLECULE_ALIASES})
        raise ValueError(
            f"Unknown molecule type {molecule_type!r}. "
            f"Choose from: {', '.join(valid)}"
        )
    return ALPHABETS[key]


# ---------------------------------------------------------------------------
# Fragment parsing
# ---------------------------------------------------------------------------


def _parse_fragment(fragment: str, alphabet: frozenset[str]) -> list[frozenset[str]]:
    """Parse a fragment pattern into a list of per-position character sets.

    Each *token* in the pattern string produces one entry in the returned
    list, representing the set of characters accepted at that residue
    position:

    - A plain character (e.g. ``'a'``) → ``frozenset({'a'})``
    - ``[ag]`` – character class → ``frozenset({'a', 'g'})``
    - ``{a}`` – negated character class → *alphabet* ``-`` ``frozenset({'a'})``
    - ``x`` – wildcard → *alphabet* (any residue)

    Matching is always case-insensitive; all characters are lowered
    internally.

    Parameters
    ----------
    fragment:
        Fragment pattern string (e.g. ``"a[cg]{p}x"``).
    alphabet:
        Set of valid lower-case characters for the molecule type.

    Returns
    -------
    list[frozenset[str]]
        One frozenset per character position described by the fragment.

    Raises
    ------
    ValueError
        If a ``[`` or ``{`` bracket is not properly closed.
    """
    token_sets: list[frozenset[str]] = []
    i = 0
    while i < len(fragment):
        ch = fragment[i]
        if ch == "[":
            try:
                j = fragment.index("]", i + 1)
            except ValueError:
                raise ValueError(f"Unclosed '[' in fragment {fragment!r} at position {i}")
            chars = frozenset(c.lower() for c in fragment[i + 1 : j])
            token_sets.append(chars)
            i = j + 1
        elif ch == "{":
            try:
                j = fragment.index("}", i + 1)
            except ValueError:
                raise ValueError(f"Unclosed '{{' in fragment {fragment!r} at position {i}")
            excluded = frozenset(c.lower() for c in fragment[i + 1 : j])
            token_sets.append(alphabet - excluded)
            i = j + 1
        elif ch.lower() == "x":
            token_sets.append(alphabet)
            i += 1
        else:
            token_sets.append(frozenset([ch.lower()]))
            i += 1
    return token_sets


# ---------------------------------------------------------------------------
# Low-level matching helpers
# ---------------------------------------------------------------------------


def _ungapped_to_aligned_index(sequence: str, ungapped_idx: int) -> int | None:
    """Return the aligned index for the n-th (0-based) non-gap character.

    Iterates through *sequence* counting non-``'-'`` characters.  When the
    *ungapped_idx*-th non-gap character is reached its position in the aligned
    string is returned.  Returns ``None`` if *ungapped_idx* is beyond the last
    non-gap character.

    Example::

        sequence    = "A-CG"
        ungapped[0] = 'A'  -> aligned index 0
        ungapped[1] = 'C'  -> aligned index 2   (gap at index 1 is skipped)
        ungapped[2] = 'G'  -> aligned index 3

    Parameters
    ----------
    sequence:
        Aligned sequence string, potentially containing ``-`` gap characters.
    ungapped_idx:
        Zero-based index in the ungapped (gap-stripped) sequence.

    Returns
    -------
    int or None
        Position in the aligned sequence, or ``None`` if *ungapped_idx* is
        beyond the last non-gap character.
    """
    count = 0
    for i, ch in enumerate(sequence):
        if ch != "-":
            if count == ungapped_idx:
                return i
            count += 1
    return None


def _match_fragment_at(
    sequence: str,
    ungapped_start: int,
    char_sets: list[frozenset[str]],
) -> bool:
    """Return ``True`` if *char_sets* matches *sequence* starting at *ungapped_start*.

    Each element of *char_sets* is checked against consecutive non-gap
    residues beginning at the *ungapped_start*-th non-gap character (0-based).

    Parameters
    ----------
    sequence:
        Aligned sequence string (may contain ``-`` gap characters).
    ungapped_start:
        Zero-based starting index in the ungapped sequence.
    char_sets:
        Per-position sets of accepted characters (from :func:`_parse_fragment`).

    Returns
    -------
    bool
    """
    for offset, char_set in enumerate(char_sets):
        aligned_idx = _ungapped_to_aligned_index(sequence, ungapped_start + offset)
        if aligned_idx is None or sequence[aligned_idx].lower() not in char_set:
            return False
    return True


def _check_motif(
    seq_id: str,
    sequence: str,
    positions: list[int],
    parsed_fragments: list[list[frozenset[str]]],
    tolerance: int = 0,
) -> dict[str, Any]:
    """Check whether a single sequence contains the motif.

    Parameters
    ----------
    seq_id:
        Identifier of the sequence.
    sequence:
        Full aligned sequence string (may contain ``-`` gap characters).
    positions:
        One-based residue positions in the *ungapped* sequence.  Gap
        characters are skipped so these positions refer to the original
        biological sequence numbering.
    parsed_fragments:
        Per-fragment list of per-position character sets as returned by
        :func:`_parse_fragment`.
    tolerance:
        Maximum allowed shift (in residues) for each fragment's start
        position.  A value of ``0`` requires exact positions.

    Returns
    -------
    dict
        ``{"sequence_id": str, "motif_present": bool}``
    """
    for pos, char_sets in zip(positions, parsed_fragments):
        ungapped_start_0 = pos - 1  # convert 1-based to 0-based
        matched = False
        for delta in range(-tolerance, tolerance + 1):
            candidate = ungapped_start_0 + delta
            if candidate < 0:
                continue
            if _match_fragment_at(sequence, candidate, char_sets):
                matched = True
                break
        if not matched:
            return {"sequence_id": seq_id, "motif_present": False}
    return {"sequence_id": seq_id, "motif_present": True}


def _check_motif_star(args: tuple) -> dict[str, Any]:
    """Unpack tuple and call :func:`_check_motif` (for multiprocessing)."""
    return _check_motif(*args)


def detect_motif_in_alignment(
    alignment_file: str | Path,
    positions: list[int],
    fragments: list[str],
    molecule_type: str = "rna",
    tolerance: int = 0,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Detect a nucleotide/amino-acid motif in each sequence of an alignment.

    The motif is defined by a list of *one-based* start positions
    (*positions*) and a corresponding list of *fragment* patterns
    (*fragments*).  Each fragment describes one or more consecutive residues
    starting at the given position.  A sequence is considered to have the
    motif when **all** fragments match (case-insensitive).

    Gap characters (``-``) introduced by the aligner are transparently skipped
    so that *positions* always refer to residue numbers in the original,
    ungapped sequence.

    **Fragment patterns** support the following syntax:

    - A plain character (``'a'``) matches exactly that residue.
    - ``[ag]`` matches any character in the bracket set (``'a'`` or ``'g'``).
    - ``{a}`` matches any character in the molecule alphabet *except* those
      listed inside the braces.
    - ``x`` matches any character in the molecule alphabet.
    - Multi-character patterns span consecutive ungapped positions, so
      ``"agc"`` at position 100 requires positions 100, 101, 102 to be
      ``'a'``, ``'g'``, ``'c'`` respectively.  This means that
      ``positions=[100, 101, 102], fragments=["a", "g", "c"]`` is
      equivalent to ``positions=[100], fragments=["agc"]``.

    Parameters
    ----------
    alignment_file:
        Path to the aligned FASTA file.
    positions:
        One-based start positions (in the *ungapped* sequence) for each
        fragment.
    fragments:
        Fragment pattern strings (same order as *positions*).
    molecule_type:
        Molecule type used to resolve wildcard alphabets.  One of
        ``"rna"`` (default), ``"dna"``, or ``"protein"`` (case-insensitive).
    tolerance:
        Maximum shift (in ungapped residues) allowed around each fragment's
        specified start position.  Defaults to ``0`` (exact match).
    n_jobs:
        Number of parallel worker processes.  Use ``1`` to disable
        multiprocessing (useful in tests).

    Returns
    -------
    pandas.DataFrame
        Table with columns ``sequence_id`` and ``motif_present`` (bool).

    Raises
    ------
    ValueError
        If *positions* and *fragments* have different lengths, or if
        *molecule_type* is not recognised.
    """
    if len(positions) != len(fragments):
        raise ValueError(
            f"positions and fragments must have the same length, "
            f"got {len(positions)} and {len(fragments)}"
        )

    alphabet = _get_alphabet(molecule_type)
    parsed_fragments = [_parse_fragment(frag, alphabet) for frag in fragments]

    alignment_file = Path(alignment_file)
    logger.info(f"Loading alignment from {alignment_file}")
    alignment: MultipleSeqAlignment = AlignIO.read(str(alignment_file), "fasta")
    logger.info(
        f"Alignment has {len(alignment)} sequences, "
        f"{alignment.get_alignment_length()} columns"
    )

    records = list(alignment)
    args_list = [
        (str(rec.id), str(rec.seq), positions, parsed_fragments, tolerance)
        for rec in records
    ]

    if n_jobs == 1:
        rows = [_check_motif(*args) for args in args_list]
    else:
        with ProcessPoolExecutor(max_workers=n_jobs) as pool:
            rows = list(pool.map(_check_motif_star, args_list))

    df = pd.DataFrame(rows, columns=["sequence_id", "motif_present"])
    n_present = df["motif_present"].sum()
    logger.info(
        f"Motif detected in {n_present}/{len(df)} sequences "
        f"(positions={positions}, fragments={fragments}, tolerance={tolerance})"
    )
    return df
