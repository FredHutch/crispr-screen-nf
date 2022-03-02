#!/bin/bash

set -Eeuo pipefail

echo "FASTQ file is ${fastq.name}"
echo "5-prime adapter sequence is ${params.adapter_5_prime}"
echo "Trimming reads to ${params.insert_length}bp in total length"

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

cutadapt \
    -m ${params.insert_length} \
    -g ${params.adapter_5_prime} \
    --length ${params.insert_length} \
    -o \$SAMPLE_NAME".fastq" "${fastq.name}" \
    
ls -lahtr
