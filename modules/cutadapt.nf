// Process : Cutadapt Trim - Fixed Number of Bases
process trim_fixed {
    
    container "${params.container__cutadapt}"
    label "mem_medium"


    publishDir "${params.output}/${prefix}/fastq/", mode: "copy", overwrite: "true"

    input:
        file(fastq)
        val(prefix)
        
    output: 
        file "${fastq.name.replaceAll(".gz","")}"

    script:
    template "cutadapt_trim.sh"
}

// Process : Cutadapt Trim - Paired Adapter Sequences
process trim_adapters {
    
    container "${params.container__cutadapt}"
    label "mem_medium"


    publishDir "${params.output}/${prefix}/fastq/", mode: "copy", overwrite: "true"

    input:
        file(fastq)
        val(prefix)
        
    output: 
        file "${fastq.name.replaceAll(".gz","")}"

    script:
    template "cutadapt_adapters.sh"
}

workflow trim {

    take:
    fastq_ch
    prefix

    main:

    // If the user provided a set of paired adapter sequences
    if ( params.paired_adapters ){

        // Perform trimming using those paired adapters
        trim_adapters(fastq_ch, prefix)

        // Use the output as that process as the output of the workflow
        trimmed_fastqs = trim_adapters.out

    // If adapter sequences were not found
    } else {

        // Perform trimming using a fixed number of bases from each end
        trim_fixed(fastq_ch, prefix)

        // Use the output as that process as the output of the workflow
        trimmed_fastqs = trim_fixed.out

    }

    emit:
    trimmed_fastqs

}