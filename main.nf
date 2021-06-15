#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.fastq = false
params.library = false
params.ntc_list = false
params.output = false
params.output_prefix = false

// Space delimited list of file endings to be removed from
// FASTQ file names to yield the samples names that they
// correspond to
params.suffix_list = "gz fq fastq fna fasta"

// Import the modules
include {
    mageck as treatment_mageck;
    mageck as control_mageck;
    join_counts;
    mageck_test;
    mageck_test_ntc;
} from './modules' params(
    suffix_list: params.suffix_list,
    output: params.output,
    output_prefix: params.output_prefix
)

// Function which prints help message text
def helpMessage() {
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
    --output_prefix     Prefix for all output files

Optional Arguments:
    --ntc_list          Path to file describing negative controls
                            As described in https://sourceforge.net/p/mageck/wiki/input/#negative-control-sgrna-list

"""
}

workflow {

    main:
        
    // If the user used the --help flag
    if(params.help){

        // Display the help text
        helpMessage()

        // And exit with an error code
        exit 1
    }

    // If the user did not specify a FASTQ input path for treatment samples
    if(!params.treatment_fastq){

        // Inform them of the error
        log.info"""

        ERROR: The --treatment_fastq flag must be provided to specify input files for treatment samples
        Use the --help flag for more details

        """
        // And exit with an error code
        exit 1
    }

    // If the user did not specify a FASTQ input path for control samples
    if(!params.control_fastq){

        // Inform them of the error
        log.info"""

        ERROR: The --control_fastq flag must be provided to specify input files for control samples
        Use the --help flag for more details

        """
        // And exit with an error code
        exit 1
    }

    // If the user did not specify a library input path
    if(!params.library){

        // Inform them of the error
        log.info"""

        ERROR: The --library flag must be provided to specify input files
        Use the --help flag for more details

        """
        // And exit with an error code
        exit 1
    }

    // If the user did not specify an output path
    if(!params.output){
        // Inform them of the error
        log.info"""

        ERROR: The --output flag must be provided to specify an output directory
        Use the --help flag for more details
        
        """
        // And exit with an error code
        exit 1
    }

    // Make a channel with the input reads for treatment samples
    Channel
        .fromPath(params.treatment_fastq.split(",").toList())
        .set{treatment_reads_ch}

    // Make a channel with the input reads for control samples
    Channel
        .fromPath(params.control_fastq.split(",").toList())
        .set{control_reads_ch}

    // Read the sgRNA library file
    Channel
        .fromPath(params.library)
        .set{sgrna_library}

    // Run MAGeCK on the treatment FASTQ files
    treatment_mageck(
        treatment_reads_ch.combine(sgrna_library)
    )

    // Run MAGeCK on the control FASTQ files
    control_mageck(
        control_reads_ch.combine(sgrna_library)
    )

    // Join together the counts from all samples
    join_counts(
        treatment_mageck.out.toSortedList(),
        control_mageck.out.toSortedList(),
    )

    // If the user supplied a list of guides used as negative controls
    if(params.ntc_list){

        // Run mageck test with the control-sgrna option
        mageck_test_ntc(
            join_counts.out.combine(
                Channel
                    .fromPath(params.ntc_list)
            )
        )

    // Otherwise
    }else{

        // Run mageck test without the control-sgrna option
        mageck_test(
            join_counts.out
        )

    }
}