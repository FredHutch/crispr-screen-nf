// Process : Build PDFs Using RMD Files
process rmd_pdf { 

    container "${params.container__rmd}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.pdf"

    input:
        tuple file("*.R"), file("*")
        val(prefix)

    output:
        path "*"

    script:
    template "rmd_pdf.sh"
}