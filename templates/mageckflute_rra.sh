#!/bin/bash
set -Eeuo pipefail

mageckflute_rra.R \
    "${gene_summary}" \
    "${sgrna_summary}" \
    "${organism}" \
    "${scale_cutoff}"

ls -lahtr