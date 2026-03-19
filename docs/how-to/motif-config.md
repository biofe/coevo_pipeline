# How to change the motif definition

The structural motif is defined in `config.yaml` under the `motif` key.  You can change positions, residue patterns, tolerance, and molecule type without modifying any code.

---

## Minimal example

```yaml
motif:
  positions: [1408, 1409, 1410]
  fragments: ["A", "C", "G"]
  molecule_type: rna
  tolerance: 0
```

---

## Change the positions

Set `motif.positions` to the one-based, ungapped positions you want to interrogate:

```yaml
motif:
  positions: [530, 531, 532]
  fragments: ["G", "A", "U"]
  molecule_type: rna
  tolerance: 0
```

---

## Use pattern syntax in fragments

Each entry in `motif.fragments` is a **pattern string** describing one or more consecutive residues:

| Syntax | Meaning |
|--------|---------|
| `A` | Exact match (residue A) |
| `[AG]` | A or G (character class) |
| `{A}` | Anything except A |
| `x` | Match any residue |
| `AGC` | Three consecutive residues A, G, C |

A single fragment can span multiple positions.  The following two configurations are **equivalent**:

```yaml
# Three separate single-character fragments
positions: [100, 101, 102]
fragments: ["A", "G", "C"]

# One three-character fragment starting at position 100
positions: [100]
fragments: ["AGC"]
```

---

## Allow positional tolerance

Set `motif.tolerance` to allow a ±N position shift when searching for each fragment.  This is useful when the alignment introduces small positional uncertainty:

```yaml
motif:
  positions: [1408]
  fragments: ["ACG"]
  molecule_type: rna
  tolerance: 2   # match anywhere within positions 1406–1410
```

---

## Change the molecule type

Set `motif.molecule_type` to `rna`, `dna`, or `protein` to control how sequences are normalised before matching:

```yaml
motif:
  molecule_type: dna   # T instead of U; no case normalisation to RNA alphabet
```

---

## Re-run motif detection after changing the config

```bash
coevo detect-motif --config config.yaml
```

Or, with Snakemake (force re-run the motif step):

```bash
snakemake --cores 8 --configfile config.yaml --forcerun detect_motif
```
