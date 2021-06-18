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
    --mle_designmat     To use MAGeCK-mle to call gene essentiality, use this flag
                            to specify the path a design matrix file as described in
                            https://sourceforge.net/p/mageck/wiki/demo/#the-fourth-tutorial-using-mageck-mle-module
    --organism          Organism string provided for MAGeCK-Flute (default: hsa)
    --scale_cuttoff     Parameter 'scale_cutoff' for MAGeCK-Flute (default: 1)
    --skip_flute        MAGeCK-Flute is only compatible with human (hsa) or mouse (mmu) gene symbols.
                        If the guide library contains gene symbols which are not compatible, set this
                        flag to skip the MAGeCK-Flute portion of the analysis.

```