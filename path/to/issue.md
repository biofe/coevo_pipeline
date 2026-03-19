## Pipeline map (Snakemake + CLI)

- [Protein BLAST](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/protein_blast_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [16S BLAST](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/16s_blast_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [Taxa Extraction](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/taxa_extraction_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [Prepare Alignment FASTA/Metadata](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/prepare_alignment_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [MAFFT Alignment](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/mafft_alignment_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [Motif Detection](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/motif_detection_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [Co-occurrence Analysis](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/cooccurrence_analysis_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.
- [Summary Statistics](https://github.com/biofe/coevo_pipeline/blob/97d718943884499206b9f28c319e7d8bbebfcc0e/path/to/summary_statistics_rule.smk#L<line_number>) - Links to the Snakemake rule and Python module.

## Note
Please update the README.md Configuration section to reflect that `motif: residues` is now `motif: fragments` and there is a `motif.tolerance` option, pointing to `config.yaml` and `coevo/config.py` that normalize legacy keys.
