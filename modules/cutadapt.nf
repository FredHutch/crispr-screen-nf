// Process : Cutadapt Trim
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

workflow trim {

    take:
    fastq_ch
    prefix

    main:

    trim_fixed(fastq_ch, prefix)

    emit:
    trim_fixed.out

}