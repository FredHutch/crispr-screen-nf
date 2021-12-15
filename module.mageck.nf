// Process : MAGeCK Count
process mageck_count {

    container "${params.container__mageck}"
    label "mem_medium"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log/", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(fastq), file(library)
        val(prefix)

    output:
        tuple file("*.R"), file("*"), emit: r
        path '*.count.txt', emit: counts
        
    script:
    template "mageck_count.sh"

}

// Process : Mageck Merge
process mageck_merge {

    container "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"
    label "io_limited"

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"

    input:
        file "treatment/treatment_*.txt"
        file "control/control_*.txt"
        val(prefix)

    output:
        tuple file("counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt"), emit: merged

    script:
    template "mageck_merge.py"
}

// Process : MAGeCK RRA
process mageck_rra {

    container "${params.container__mageck}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log/", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(counts_tsv), file(treatment_samples), file(control_samples)
        val(output_prefix)
        val(prefix)

    output:
        tuple file("*.R"), file("*"), emit: r
        path "${output_prefix}.gene_summary.txt", emit: geneSummary
        path "${output_prefix}.sgrna_summary.txt", emit: sgrnaSummary
        path "${output_prefix}.normalized.txt", emit: normCounts
        
    script:
    template "mageck_rra.sh"

}

// Process : MAGeCK MLE
process mageck_mle {

    container "${params.container__mageck}"
    label "mem_medium" 

    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*.txt"
    publishDir "${params.output}/${prefix}/log", mode: "copy", overwrite: "true", pattern: "*.log"

    input:
        tuple file(counts_tsv), file(sample_names), file(control_names)
        val(treatment)
        val(control)
        val(output_prefix)
        val(prefix)

    output:
        path "${output_prefix}.gene_summary.txt", emit: geneSummary
        path "${output_prefix}.sgrna_summary.txt", emit: sgrnaSummary

    script:
    template "mageck_mle.sh"

}

// Process : MAGeCK Pathway
process mageck_pathway { 
    
    container "${params.container__mageck}"
    label "mem_medium"
    
    publishDir "${params.output}/${prefix}/", mode: "copy", overwrite: "true", pattern: "*"
    
    output:
        path "gene.pathway_summary.txt", emit: gene_summary

    input:
        tuple file(counts_tsv), file(gmt)
        val(output_prefix)
        val(prefix)

    script:
    template "mageck_pathway.sh"
    
}