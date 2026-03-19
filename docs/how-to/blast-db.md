# How to use a different BLAST database

By default the pipeline searches against NCBI's `nr` (protein) and `core_nt` (nucleotide) databases.  You can point it at any local BLAST-formatted database.

---

## Change the database path in config.yaml

Edit the `blast` section:

```yaml
blast:
  protein_db: /data/blast/nr          # absolute or relative path, or a BLAST alias
  nucleotide_db: /data/blast/core_nt
  threads: 16
  evalue: 0.001
  max_target_seqs: 10000
```

`protein_db` and `nucleotide_db` accept any value that BLAST's `-db` flag accepts:

- A local path to the database files (without extension): `/data/blast/nr`
- A BLAST alias file: `nt` (if `$BLASTDB` environment variable is set)

---

## Build a custom local database

### Protein database

```bash
makeblastdb -in my_proteins.fasta -dbtype prot -out my_prot_db -parse_seqids
```

```yaml
blast:
  protein_db: my_prot_db
```

### Nucleotide database

```bash
makeblastdb -in my_16s_sequences.fasta -dbtype nucl -out my_nucl_db -parse_seqids
```

```yaml
blast:
  nucleotide_db: my_nucl_db
```

!!! warning "taxid support"
    The pipeline extracts taxids from the `staxids` BLAST output field.  Custom databases must be built with taxonomy information (e.g., using `-taxid_map`) for taxid extraction to work correctly.  Without taxids, `coevo extract-taxa` will produce empty taxa files.

---

## Set the BLASTDB environment variable

If your databases are in a shared directory, set `$BLASTDB` so BLAST can locate them by name:

```bash
export BLASTDB=/data/blast
```

Then use just the database name in `config.yaml`:

```yaml
blast:
  protein_db: nr
  nucleotide_db: core_nt
```

---

## Verify the database is usable

```bash
blastdbcheck -db /data/blast/nr -dbtype prot
blastdbcheck -db /data/blast/core_nt -dbtype nucl
```
