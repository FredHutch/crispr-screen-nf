// Process MAGeCKFlute RRA
process mageckflute_rra {

    container "${params.container__mageckflute}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    
    input:
        file(gene_summary)
        file(sgrna_summary)
        val(scale_cutoff)
        val(organism)
        val(prefix)
        
    output:
        path "*"
        
    script:
    template "mageckflute_rra.sh"
}

// Process : MAGeCKFlute Mle
process mageckflute_mle {

    container "${params.container__mageckflute}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    

    input:
        file(gene_summary)
        file(depmap_effect)
        file(depmap_samples)
        val(prefix)

    output:
        path "*"

    script:
    template "mageckflute_mle.sh"
}