// Process : MAGeCKVispr Export JSON
process mageckvispr_export {

    container "${params.container__mageckvispr}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    

    input:
        file("gene_summary.txt")
        file("count_normalized.txt")
        val(prefix)

    output:
        path "*"

    script:
    template "mageckvispr_export.py"
}