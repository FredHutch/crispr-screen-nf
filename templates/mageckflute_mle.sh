#!/bin/bash
set -Eeuo pipefail

mageckflute_mle.R \
    "${gene_summary}" \
    "${params.organism}" \
    "treatment_control" \
    "deepmap" \
    "${deepmap_effect}" \
    "${deepmap_samples}"
    
ls -lahtr
