// Container used to run mageck
mageck_container = "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"
mageckflute_container = "quay.io/biocontainers/bioconductor-mageckflute:1.12.0--r41hdfd78af_0"

// Process used to run MAGeCK count
process mageck {
    container "${mageck_container}"
    label "mem_medium"
    publishDir "${params.output}/count/", mode: "copy", overwrite: "true"

    input:
        tuple file(fastq), file(library)

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

// Process used to run MAGeCK test
process mageck_test {
    container "${mageck_container}"
    label "io_limited"
    publishDir "${params.output}/rra/", mode: "copy", overwrite: "true"

    input:
        tuple file(counts_tsv), file(treatment_samples), file(control_samples)

    output:
        file "${params.output_prefix}.*"

    script:
"""/bin/bash

set -Eeuo pipefail

mageck test \
    -k ${counts_tsv} \
    -t "\$(cat ${treatment_samples})" \
    -c "\$(cat ${control_samples})" \
    -n "${params.output_prefix}"

ls -lahtr
"""

}

// Process used to run MAGeCK test with the --control-sgrna option
process mageck_test_ntc {
    container "${mageck_container}"
    label "io_limited"
    publishDir "${params.output}/rra/", mode: "copy", overwrite: "true"

    input:
        tuple file(counts_tsv), file(treatment_samples), file(control_samples), file(ntc_list)

    output:
        file "${params.output_prefix}.*"

    script:
"""/bin/bash

set -Eeuo pipefail

mageck test \
    -k ${counts_tsv} \
    -t "\$(cat ${treatment_samples})" \
    -c "\$(cat ${control_samples})" \
    -n "${params.output_prefix}" \
    --control-sgrna ${ntc_list} \
    --norm-method control

ls -lahtr
"""

}

// Process used to run MAGeCK-mle test with the --mle_designmat option
process mageck_test_mle {
    container "${mageck_container}"
    label "io_limited"
    publishDir "${params.output}/mle/", mode: "copy", overwrite: "true"

    input:
        tuple file(counts_tsv), file(treatment_samples), file(control_samples), file(designmat)

    output:
        file "${params.output_prefix}.*"

    script:
"""/bin/bash

set -Eeuo pipefail

mageck mle \
    -k ${counts_tsv} \
    -d ${designmat} \
    -n "${params.output_prefix}"

ls -lahtr
"""

}

// Process used to run MAGeCK FluteRRA
process mageck_flute_rra {
    container "${mageckflute_container}"
    label "io_limited"
    publishDir "${params.output}/rra_flute/", mode: "copy", overwrite: "true"

    input:
        file "*"

    output:
        file "*"

    script:
"""/usr/bin/env R --no-save

# Set up the paths to the input files
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"), "${params.output_prefix}.gene_summary.txt")

file2 = file.path(system.file("extdata", package = "MAGeCKFlute"), "${params.output_prefix}.sgrna_summary.txt")

# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
FluteRRA(file1, file2, proj="${params.output_prefix}", organism="${params.organism}", scale_cutoff = ${params.scale_cutoff}, outdir = "./")

"""

}

// Process used to join the outputs from mageck / counts
process join_counts {
    container "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"
    label "io_limited"

    input:
        file "treatment/treatment_*.txt"
        file "control/control_*.txt"

    output:
        tuple file("${params.output_prefix}.counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt")

    script:
"""
set -Eeuo pipefail

join_counts.py "${params.output_prefix}"

"""

}
