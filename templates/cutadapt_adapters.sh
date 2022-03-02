#!/bin/bash

set -Eeuo pipefail

echo FASTQ file is ${fastq.name}
echo Adapter sequences are ${params.paired_adapters}

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

cutadapt \
    -m ${params.min_readlen} \
    -M ${params.max_readlen} \
    -a ${params.paired_adapters} \
    -o \$SAMPLE_NAME".fastq" "${fastq.name}" \
    
ls -lahtr
