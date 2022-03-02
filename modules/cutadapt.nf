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

    // If the user provided a 5' adapter sequence
    if ( params.adapter_5_prime ){

        // Perform trimming using those paired adapters
        // Also trim the read length down to params.insert_length
        trim_adapters(fastq_ch, prefix)

        // Use the output as that process as the output of the workflow
        trimmed_fastqs = trim_adapters.out

    // If adapter sequences were not found
    } else {

        // Perform trimming using a fixed number of bases (--trim_5_prime, defaults to 0)
        // Also trim the read length down to params.insert_length
        trim_fixed(fastq_ch, prefix)

        // Use the output as that process as the output of the workflow
        trimmed_fastqs = trim_fixed.out

    }

    emit:
    trimmed_fastqs

}