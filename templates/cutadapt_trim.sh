#!/bin/bash

set -Eeuo pipefail

echo FASTQ file is ${fastq.name}
echo Prime5 is ${params.trim_5_prime}
echo Prime3 is ${params.trim_3_prime}

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

cutadapt \
    -u "${params.trim_5_prime}" \
    -u "${params.trim_3_prime}" \
    -o \$SAMPLE_NAME".fastq" "${fastq.name}" \
    
ls -lahtr
