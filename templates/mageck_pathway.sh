#!/bin/bash
set -Eeuo pipefail

mageck pathway \
    --gene-ranking ${counts_tsv} \
    --gmt-file ${gmt} \
    -n "${output_prefix}"

ls -lahtr