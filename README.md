# crispr-screen-nf
Nextflow workflow for analyzing CRISPR screen NGS data

## Background

The purpose of this workflow is to analyze NGS data generated from a CRISPR screen,
which consists of two groups of FASTQ files. One of the groups of FASTQ files is
considered a 'control', meaning that it was generated from a set of CRISPR sgRNA
guides which did not undergo an experimental selection process. The other group of
FASTQ files are the 'treatment', indicating that they did undergo some experimental
selection. The purpose of this analysis is to quantify and summarize the degree of
that experimental selection process on each of the sgRNA guides in the library. In
other words, to determine what the biological impact of knocking down each particular
gene was in this biological context.

## Library Definition

In addition to the two groups of FASTQ files, the user must provide a 'library'
file which describes the set of sgRNA guides used in the screen, along with the
genes which each guide targets. Documentation describing this library file can
be found [here](https://sourceforge.net/p/mageck/wiki/input/). The file format
expected by this particular workflow is CSV, with the three column headers named:
`guide,sgrna,gene`.

## Read Trimming

After generating NGS data from the sgRNA libraries used in this screen, the user
must specify the number of bases which need to be trimmed from the beginning (5')
and end (3') of each read in order to yield the exact region of the sgRNA sequence.
Depending on the read length and library design, these values may vary between
experiments. The default values are to trim 32bp from the beginning of each read,
and 8bp from the end. Note that the number of bases trimmed from the end of the
read is indicated with a negative number (as `-8`).

## Prerequisites

Before using this pipeline, the user must install Nextflow and configure Nextflow
to use the appropriate set of computational resources which are available. The
software dependencies for this pipeline are managed using Docker images (instead
of Conda environments), and so support must be added for Docker or Singularity
as appropriate. The usage documentation listed below assumes a funcational
Nextflow configuration present on the system used to run the workflow.

## Usage Documentation

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

    --trim_5_prime      Amount to trim from 5 prime end (default: 32)
    --trim_3_prime      Amount to trim from 3 prime end (default: -8)
    --organism          Organism string provided for MAGeCK-Flute (default: hsa)
    --scale_cutoff      Parameter 'scale_cutoff' for MAGeCK-Flute (default: 1)
    --gmt               Pathway GMT File
    --deepmap_effect    "https://ndownloader.figshare.com/files/20234073" - Effect File Can't Download From Different Region
    --deepmap_samples   "https://ndownloader.figshare.com/files/20274744" - Sample File Can't Download From Different Region
    --ntc_list          Path to file describing negative controls
                            As described in https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list
```