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

// Process used to run MAGeCK
process mageck {
    container "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"
    label "mem_medium"

    input:
        file fastq
        file library

    output:
        file "*.count.txt"

    script:
"""/bin/bash

set -Eeuo pipefail

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

echo FASTQ file is ${fastq.name}
echo Sample name is \$SAMPLE_NAME

mageck count -l ${library} -n \$SAMPLE_NAME --sample-label \$SAMPLE_NAME  --fastq ${fastq}
"""

}

workflow {
        
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
    sgrna_library = file(params.library)

    // Run MAGeCK on the FASTQ files
    mageck(
        reads_ch,
        sgrna_library
    )
}