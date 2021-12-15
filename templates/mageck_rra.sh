#!/bin/bash
set -Eeuo pipefail

echo COUNTS file is ${counts_tsv.name}
echo TREATMENT file is \$(cat ${treatment_samples})
echo CONTROLS file is \$(cat ${control_samples})
echo OUTPUT prefix is ${output_prefix}

mageck test \
    -k ${counts_tsv} \
    -t "\$(cat ${treatment_samples})" \
    -c "\$(cat ${control_samples})" \
    -n "${output_prefix}" \
    --normcounts-to-file \
    --keep-tmp \
    --pdf-report

ls -lahtr