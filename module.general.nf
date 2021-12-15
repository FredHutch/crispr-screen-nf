// Process : Populate duplicated NTCs in the count table
// Context : Sometimes the input library will have duplicated guides used to generate
// synthetic NTC genes. Those duplicated guides will be lost by mageck count (only
// keeping one of the unique values). Using this process, we will add back in the counts
// for any duplicated guides
process populate_ntc {

    container "${params.container__pandas}"
    label "io_limited"

    input:
        tuple file("input.counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt"), file(library)

    output:
        tuple file("counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt")

    script:
    template 'populate_ntc.py'
}

// Assess quality of reads after removing barcodes
process fastqc {
    container "${params.container__fastqc}"
    label "io_limited"
    publishDir "${params.output}/fastqc/", mode: "copy", overwrite: "true", pattern: "*/*.html"
    
    input:
    path fastq

    output:
    path "*/*.html", emit: html
    path "*/*.zip", emit: zip
    path "*/*.txt", emit: data

    script:
    template 'fastqc.sh'

}

// Format the FastQC output as JSON
process parse_fastqc {
    container "${params.container__pandas}"
    label "io_limited"
    publishDir "${params.output}/fastqc/${fastqc_txt.name.replaceAll(/.txt/, '')}/", mode: "copy", overwrite: "true"
    
    input:
    path fastqc_txt

    output:
    path "*.json"

    script:
    template 'parse_fastqc.sh'

}

// Combine all FASTQC data into a single report
process multiqc {
    container "${params.container__multiqc}"
    publishDir "${params.output}/fastqc/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    template 'multiqc.sh'

}