#!/bin/bash

set -Eeuo pipefail

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

echo FASTQ file is ${fastq.name}
echo Sample name is \$SAMPLE_NAME

mageck count \
    -l ${library} \
    -n \$SAMPLE_NAME \
    --sample-label \$SAMPLE_NAME \
    --fastq ${fastq} \
    --trim-5 0 \
    --pdf-report

ls -lahtr