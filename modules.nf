// Container used to run mageck
mageck_container = "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"

// Process used to run MAGeCK count
process mageck {
    container "${mageck_container}"
    label "mem_medium"

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
    publishDir "${params.output}"

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
    publishDir "${params.output}"

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

// Process used to join the outputs from mageck / counts
process join_counts {
    container "quay.io/fhcrc-microbiome/python-pandas:v1.2.1_latest"
    label "io_limited"

    input:
        file "treatment/*"
        file "control/*"

    output:
        tuple file("${params.output_prefix}.counts.txt"), file("treatment_sample_names.txt"), file("control_sample_names.txt")

    script:
"""
set -Eeuo pipefail

join_counts.py "${params.output_prefix}"

"""

}
