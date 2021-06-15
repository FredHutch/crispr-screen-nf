# crispr-screen-nf
Nextflow workflow for analyzing CRISPR screen NGS data

```
Usage:
nextflow run FredHutch/crispr-screen-nf

Required Arguments:
    --treatment_fastq   Path to FASTQ data for treatment samples
    --control_fastq     Path to FASTQ data for control samples
                        Multiple files can be specified with wildcards and commas, e.g.
                            /path/to/inputs/A/*.fastq.gz,/path/to/inputs/B/*.fq.gz
    --library           Text file describing sgRNA library
                            As described at https://sourceforge.net/p/mageck/wiki/input/
    --output            Path to output directory
    --output_prefix     Prefix for all output files

Optional Arguments:
    --ntc_list          Path to file describing negative controls
                            As described in https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list

```