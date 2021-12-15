#!/bin/bash

set -Eeuo pipefail

# Parse the name of the sample from the name of the FASTQ file
SAMPLE_NAME="${fastq.name}"
for suffix in ${params.suffix_list}; do
    SAMPLE_NAME=\$(echo \$SAMPLE_NAME | sed "s/.\$suffix\$//")
done

echo "Creating output folder"
mkdir \$SAMPLE_NAME

echo "Running FASTQC"
fastqc --threads ${task.cpus} --extract -o \$SAMPLE_NAME "$fastq"

# Rename the FASTQC data file to contain the sample name
mv \$SAMPLE_NAME/\${SAMPLE_NAME}_fastqc/fastqc_data.txt \$SAMPLE_NAME/\$SAMPLE_NAME.txt

echo "DONE"