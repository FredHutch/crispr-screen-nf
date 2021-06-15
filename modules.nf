
// Process used to run MAGeCK
process mageck {
    container "quay.io/biocontainers/mageck:0.5.9.4--py38h8c62d01_1"
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
