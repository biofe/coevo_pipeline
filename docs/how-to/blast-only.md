# How to run only the BLAST steps

Use this guide if you want to run just the database searches (BLASTP and/or BLASTN) without proceeding to alignment and motif detection.

---

## Run BLASTP only

```bash
coevo blast-protein examples/protein_query.fasta --config config.yaml
```

Output: `results/blast/protein_hits.tsv`

---

## Run BLASTN only

```bash
coevo blast-16s examples/small_16s.fasta --config config.yaml
```

Output: `results/blast/rna_hits.tsv`

---

## Run both BLAST steps and stop

```bash
coevo blast-protein examples/protein_query.fasta --config config.yaml
coevo blast-16s     examples/small_16s.fasta     --config config.yaml
```

---

## Using Snakemake targets

You can ask Snakemake to produce only the BLAST output files by specifying them as explicit targets:

```bash
snakemake --cores 8 --configfile config.yaml \
    results/blast/protein_hits.tsv \
    results/blast/rna_hits.tsv
```

Snakemake will only execute the rules needed to produce those files.

---

## What to do with the BLAST outputs

Once you have the hit tables you can:

1. **Extract taxids** — run `coevo extract-taxa` when ready.
2. **Inspect hits** — open `results/blast/protein_hits.tsv` in any TSV viewer (pandas, Excel, etc.).
3. **Adjust thresholds** — tune `filters.min_identity` in `config.yaml`, then re-run.

See [Config reference](../reference/config.md) for all filtering options.
