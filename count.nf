#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the modules
include { trim } from './modules/cutadapt'

include { mageck_count } from './module.mageck'

include {
    populate_ntc;
    fastqc;
    parse_fastqc;
    multiqc;
} from './module.general'

// Validate Input Parameters / Print Help
def validate(params) {

    if (params.fastq
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
            --fastq             Path to FASTQ data for treatment samples
            --library           Text file describing sgRNA library
                                    As described at https://sourceforge.net/p/mageck/wiki/input/
            --output            Path to output directory

            --insert_length     Length of sgRNA guides (default: 20)
            --trim_5_prime      Number of bases to trim from the 5' of each read (default: 0)
            --adapter_5_prime   (optional) Sequence of 5' adapter to be trimmed from each read
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
    Channel.fromPath(params.fastq.split(',').toList()).set{Channel_Fastq}

    // Chanel : SGRNA Library
    Channel.fromPath(params.library).set{Channel_Library}

    // Process : FastQC
    // Compute quality control metrics for all input FASTQ data
    fastqc(Channel_Fastq)

    // Process : Parse FastQC Output
    parse_fastqc(fastqc.out.data)

    // Process : MultiQC
    // Combine quality control reporting across all inputs
    multiqc(fastqc.out.zip.toSortedList())

    // Process : Cutadapt Trim
    // Remove Extra Adaptor Sequences From Reads
    trim(Channel_Fastq, 'cutadapt/trim/')

    // Process : Mageck Count
    // Map the raw FASTQ data to reference library file and count the reads for each sgRNA
    mageck_count(trim.out.combine(Channel_Library), 'mageck/count')

}
