#!/bin/bash
set -Eeuo pipefail

mageckflute_mle.R \
    "${gene_summary}" \
    "${params.organism}" \
    "treatment_control" \
    "depmap" \
    "${depmap_effect}" \
    "${depmap_samples}"
    
ls -lahtr