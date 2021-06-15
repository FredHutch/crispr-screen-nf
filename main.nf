#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.fastq = false
params.library = false
params.output = false

// Space delimited list of file endings to be removed from
// FASTQ file names to yield the samples names that they
// correspond to
params.suffix_list = "gz fq fastq fna fasta"

// Import the modules
include {
    mageck
} from './modules' params(
    suffix_list: params.suffix_list
)

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:
nextflow run FredHutch/crispr-screen-nf

Required Arguments:
    --fastq             Path to input FASTQ data
                        Multiple files can be specified with wildcards and commas, e.g.
                            /path/to/inputs/A/*.fastq.gz,/path/to/inputs/B/*.fq.gz
    --library           Text file describing sgRNA library
                            As described at https://sourceforge.net/p/mageck/wiki/input/
    --output            Path to output directory

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

    // If the user did not specify a FASTQ input path
    if(!params.fastq){

        // Inform them of the error
        log.info"""

        ERROR: The --fastq flag must be provided to specify input files
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

    // Make a channel with the input reads
    Channel
        .fromPath(params.fastq.split(","))
        .set{reads_ch}

    // Read the sgRNA library file
    Channel
        .fromPath(params.library)
        .set{sgrna_library}

    // Run MAGeCK on the FASTQ files
    mageck(
        reads_ch.combine(sgrna_library)
    )
}