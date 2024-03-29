#!/bin/bash
set -Eeuo pipefail

echo CONTROL SGRNA file is ${control_sgrna}

echo 'Samples	baseline	treatment_control' > design.mtx
for i in \$(echo ${treatment} | tr ',' '\n')
do
    echo \$(basename \$i | cut -d. -f1)'  1    1' >> design.mtx
done
for i in \$(echo ${control} | tr ',' '\n')
do
    echo \$(basename \$i | cut -d. -f1)'  1    0' >> design.mtx
done

mageck mle \
    -k ${counts_tsv} \
    -d 'design.mtx' \
    -n ${output_prefix} \
    --norm-method control \
    --control-sgrna "${control_sgrna}" \
    --threads ${task.cpus}
    
ls -lahtr