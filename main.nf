#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.scale_cutoff = 1
params.trim_3_prime = -8
params.trim_5_prime = 32

// List of extensions to be remove to determine sample name
params.suffix_list = "trimmed gz fq fastq fna fasta"

// Containers
params.container__pandas = "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"
params.container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
params.container__multiqc = "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
params.container__cutadapt = "quay.io/biocontainers/cutadapt:3.4--py37h73a75cf_1"
params.container__mageck = "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"
params.container__mageckflute = "quay.io/biocontainers/bioconductor-mageckflute:1.12.0--r41hdfd78af_0"
params.container__mageckvispr = "quay.io/biocontainers/mageck-vispr:0.5.6--py_0"
params.container__rmd = "rocker/r-rmd:latest"

// Import the modules
include {
    cutadapt_trim as cutadapt_trim_treatment;
    cutadapt_trim as cutadapt_trim_control;
} from './module.cutadapt'

include {
    mageck_count as mageck_count_treatment;
    mageck_count as mageck_count_control;
    mageck_rra;
    mageck_mle;
    mageck_pathway;
    mageck_merge;
} from './module.mageck'

include {
    mageckflute_rra;
    mageckflute_mle
} from './module.mageckflute'

include {
    rmd_pdf as Rmd_Pdf_Treatment;
    rmd_pdf as Rmd_Pdf_Control;
} from './module.rmd'

include {
    populate_ntc;
    fastqc;
    parse_fastqc;
    multiqc;
} from './module.general'

include {
    mageckvispr_export as mageckvispr_export_rra;
    mageckvispr_export as mageckvispr_export_mle;
} from './module.mageckvispr'

// Validate Input Parameters / Print Help
def validate(params) {

    if (params.treatment_fastq 
        && params.control_fastq 
        && params.library
        && params.organism 
        && params.scale_cutoff 
        && params.gmt 
        && params.output)
        { return; }
    
    log.info"""
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
            --depmap_effect    "https://ndownloader.figshare.com/files/20234073" - Effect File Can't Download From Different Region
            --depmap_samples   "https://ndownloader.figshare.com/files/20274744" - Sample File Can't Download From Different Region
            --ntc_list          Path to file describing negative controls
                                    As described in https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list
        """
    
    exit 1
}

workflow {

    main:

    // Validate Input + Print Help On Fail 
    validate(params)
    
    // Channel : Treatment FastQ
    Channel.fromPath(params.treatment_fastq.split(',').toList()).set{Channel_Fastq_Treatment}

    // Chanel : Control FastQ 
    Channel.fromPath(params.control_fastq.split(',').toList()).set{Channel_Fastq_Control}

    // Chanel : SGRNA Library
    Channel.fromPath(params.library).set{Channel_Library}

    // Process : FastQC
    // Compute quality control metrics for all input FASTQ data
    fastqc(Channel_Fastq_Control.mix(Channel_Fastq_Treatment))

    // Process : Parse FastQC Output
    parse_fastqc(fastqc.out.data)

    // Process : MultiQC
    // Combine quality control reporting across all inputs
    multiqc(fastqc.out.zip.toSortedList())

    // Process : Cutadapt Trim
    // Remove Extra Adaptor Sequences From Reads
    cutadapt_trim_treatment(Channel_Fastq_Treatment, 'cutadapt/trim/treatment')
    cutadapt_trim_control(Channel_Fastq_Control, 'cutadapt/trim/control')

    // Process : Mageck Count
    // Map the raw FASTQ data to reference library file and count the reads for each sgRNA
    mageck_count_treatment(cutadapt_trim_treatment.out.combine(Channel_Library), 'mageck/count/treatment')
    mageck_count_control(cutadapt_trim_control.out.combine(Channel_Library), 'mageck/count/control')

    // Process : Mageck Merge
    // Concat All Count Data
    mageck_merge(mageck_count_treatment.out.counts.toSortedList(), mageck_count_control.out.counts.toSortedList(), 'mageck/count/combined')

    // Add back counts for duplicated guides (if any)
    populate_ntc(mageck_merge.out.merged.combine(Channel_Library))
    
    // Process : Mageck Rra
    // MAGeCK RRA (identifying CRISPR screen hits by calculating the RRA enrichment score to indicate the essentiality of a gene)
    mageck_rra(populate_ntc.out, 'rra', 'mageck/rra')

    // Process : Mageck MLE
    // MAGeCK MLE (identifying CRISPR screen hits by calculating a ‘beta score’ for each targeted gene to measure the degree of selection after the target is perturbed)
    mageck_mle(populate_ntc.out, params.treatment_fastq, params.control_fastq, 'mle', 'mageck/mle')

    // Process : Mageck Pathway
    mageck_pathway(mageck_rra.out.geneSummary.combine(Channel.fromPath(params.gmt)), 'gene', 'mageck/rra/pathway')
    
    // Run Mageck Flute RRA
    mageckflute_rra(mageck_rra.out.geneSummary, mageck_rra.out.sgrnaSummary, params.scale_cutoff, params.organism, 'mageckflute/rra')

    // Process Mageck Flute MLE
    mageckflute_mle(mageck_mle.out.geneSummary, Channel.fromPath(params.depmap_effect), Channel.fromPath(params.depmap_samples), 'mageckflute/mle')

    // Export Mageck RRA as JSON via MAGeCKVispr
    mageckvispr_export_rra(
        mageck_rra.out.geneSummary,
        mageck_rra.out.normCounts,
        "mageck/rra/vispr"
    )

    // Export Mageck RRA as JSON via MAGeCKVispr
    mageckvispr_export_mle(
        mageck_mle.out.geneSummary,
        mageck_rra.out.normCounts,
        "mageck/mle/vispr"
    )
    
    // Process : Rmd To Pdf
    // Rmd_Pdf_Treatment(mageck_count_treatment.out.r, 'mageck/count/treatment/pdf')
    // Rmd_Pdf_Control(mageck_count_control.out.r, 'mageck/count/control/pdf')
   
}
