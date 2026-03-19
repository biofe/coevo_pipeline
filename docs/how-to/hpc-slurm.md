# How to run on HPC (Slurm) with Snakemake

For large BLAST searches, running on an HPC cluster with Slurm is recommended.

---

## Prerequisites

- Snakemake ≥ 8.0 (uses the plugin-based executor interface)
- `snakemake-executor-plugin-slurm` installed in your environment

```bash
pip install snakemake-executor-plugin-slurm
```

---

## Basic Slurm submission

```bash
snakemake \
    --cores 8 \
    --configfile config.yaml \
    --executor slurm \
    --default-resources \
        slurm_partition=main \
        mem_mb=16000 \
        runtime=120
```

This submits each Snakemake rule as a separate Slurm job using the `main` partition, 16 GB RAM, and a 2-hour wall-time limit.

---

## Per-rule resource overrides

Create a `profile/` directory with a `config.yaml` for rule-level resources:

```
profile/
└── config.yaml
```

```yaml
# profile/config.yaml
executor: slurm
default-resources:
  slurm_partition: main
  mem_mb: 8000
  runtime: 60

set-resources:
  blast_protein:
    slurm_partition: high-mem
    mem_mb: 64000
    runtime: 480
    cpus_per_task: 16
  blast_16s:
    slurm_partition: high-mem
    mem_mb: 64000
    runtime: 480
    cpus_per_task: 16
```

Then run with the profile:

```bash
snakemake --profile profile/ --configfile config.yaml
```

---

## Adjust BLAST threads to match the allocated CPUs

Set `blast.threads` in `config.yaml` to match the `cpus_per_task` you request from Slurm:

```yaml
blast:
  threads: 16
```

---

## Conda environments on HPC

If your cluster does not have BLAST+ and MAFFT in the default `$PATH`, use the conda environment:

```bash
snakemake --profile profile/ --configfile config.yaml --use-conda
```

Snakemake will create the conda environment from `environment.yml` on the first run.

---

## Monitor jobs

```bash
squeue -u $USER
snakemake --profile profile/ --configfile config.yaml --summary
```
