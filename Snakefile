"""
Snakemake workflow for the coevo_rnase_16s pipeline.

Steps:
    1. BLAST protein search
    2. BLAST 16S rRNA search
    3. Parse hits & extract taxids
    4. Retrieve 16S sequences (placeholder – relies on external fetch step)
    5. Align 16S sequences with MAFFT
    6. Motif detection
    7. Co-occurrence analysis
    8. Summary statistics

Usage:
    snakemake --cores 8 --configfile config.yaml
"""

configfile: "config.yaml"

RESULTS = config["output"]["results_dir"]
CONDA_ENV = config.get("conda_env", "environment.yml")
CONFIG_FILE = workflow.configfiles[0] if workflow.configfiles else "config.yaml"


include: "workflows/rules_blast.smk"
include: "workflows/rules_alignment.smk"
include: "workflows/rules_analysis.smk"


rule all:
    input:
        expand("{results}/analysis/cooccurrence.tsv", results=RESULTS),
        expand("{results}/analysis/phylum_summary.tsv", results=RESULTS),
        expand("{results}/analysis/motif_results.tsv", results=RESULTS),
