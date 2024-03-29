#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the modules
include {
    trim as trim_treatment;
    trim as trim_control;
} from './modules/cutadapt'

include {
    mageck_count as mageck_count_treatment;
    mageck_count as mageck_count_control;
    mageck_rra;
    mageck_rra_ntc;
    mageck_mle;
    mageck_mle_ntc;
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

            --insert_length     Length of sgRNA guides (default: 20)
            --trim_5_prime      Number of bases to trim from the 5' of each read (default: 0)
            --adapter_5_prime   (optional) Sequence of 5' adapter to be trimmed from each read
            --organism          Organism string provided for MAGeCK-Flute (default: hsa)
            --scale_cutoff      Parameter 'scale_cutoff' for MAGeCK-Flute (default: 1)
            --gmt               Pathway GMT File
            --depmap_effect    "https://ndownloader.figshare.com/files/20234073" - Effect File Can't Download From Different Region
            --depmap_samples   "https://ndownloader.figshare.com/files/20274744" - Sample File Can't Download From Different Region
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
    trim_treatment(Channel_Fastq_Treatment, 'cutadapt/trim/treatment')
    trim_control(Channel_Fastq_Control, 'cutadapt/trim/control')

    // Process : Mageck Count
    // Map the raw FASTQ data to reference library file and count the reads for each sgRNA
    mageck_count_treatment(trim_treatment.out.combine(Channel_Library), 'mageck/count/treatment')
    mageck_count_control(trim_control.out.combine(Channel_Library), 'mageck/count/control')

    // Process : Mageck Merge
    // Concat All Count Data
    mageck_merge(mageck_count_treatment.out.counts.toSortedList(), mageck_count_control.out.counts.toSortedList(), 'mageck/count/combined')

    // Add back counts for duplicated guides (if any)
    populate_ntc(mageck_merge.out.merged.combine(Channel_Library))

    // If the user selected the flag for NTC-based normalization
    if ( params.use_control_normalization ){

        // Chanel : Control sgRNA list (NTCs)
        Channel.fromPath(params.control_sgrna).set{Channel_Control_Sgrna}

        // Process : Mageck Rra
        // MAGeCK RRA (identifying CRISPR screen hits by calculating the RRA enrichment score to indicate the essentiality of a gene)
        mageck_rra_ntc(populate_ntc.out, 'rra', 'mageck/rra', Channel_Control_Sgrna)
        mageck_rra_out = mageck_rra_ntc.out

        // Process : Mageck MLE
        // MAGeCK MLE (identifying CRISPR screen hits by calculating a ‘beta score’ for each targeted gene to measure the degree of selection after the target is perturbed)
        mageck_mle_ntc(populate_ntc.out, params.treatment_fastq, params.control_fastq, 'mle', 'mageck/mle', Channel_Control_Sgrna)
        mageck_mle_out = mageck_mle_ntc.out

    } else {

        // Process : Mageck Rra
        // MAGeCK RRA (identifying CRISPR screen hits by calculating the RRA enrichment score to indicate the essentiality of a gene)
        mageck_rra(populate_ntc.out, 'rra', 'mageck/rra')
        mageck_rra_out = mageck_rra.out

        // Process : Mageck MLE
        // MAGeCK MLE (identifying CRISPR screen hits by calculating a ‘beta score’ for each targeted gene to measure the degree of selection after the target is perturbed)
        mageck_mle(populate_ntc.out, params.treatment_fastq, params.control_fastq, 'mle', 'mageck/mle')
        mageck_mle_out = mageck_mle.out

    }

    // Process : Mageck Pathway
    mageck_pathway(mageck_rra_out.geneSummary.combine(Channel.fromPath(params.gmt)), 'gene', 'mageck/rra/pathway')
    
    // Run Mageck Flute RRA
    mageckflute_rra(mageck_rra_out.geneSummary, mageck_rra_out.sgrnaSummary, params.scale_cutoff, params.organism, 'mageckflute/rra')

    // Process Mageck Flute MLE
    mageckflute_mle(mageck_mle_out.geneSummary, Channel.fromPath(params.depmap_effect), Channel.fromPath(params.depmap_samples), 'mageckflute/mle')

    // Export Mageck RRA as JSON via MAGeCKVispr
    mageckvispr_export_rra(
        mageck_rra_out.geneSummary,
        mageck_rra_out.normCounts,
        "mageck/rra/vispr"
    )

    // Export Mageck RRA as JSON via MAGeCKVispr
    mageckvispr_export_mle(
        mageck_mle_out.geneSummary,
        mageck_rra_out.normCounts,
        "mageck/mle/vispr"
    )
    
    // Process : Rmd To Pdf
    // Rmd_Pdf_Treatment(mageck_count_treatment.out.r, 'mageck/count/treatment/pdf')
    // Rmd_Pdf_Control(mageck_count_control.out.r, 'mageck/count/control/pdf')
   
}
