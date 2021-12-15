// Process : Cutadapt Trim
process cutadapt_trim {
    
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