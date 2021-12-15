#!/bin/bash
set -Eeuo pipefail

for rmd in *.R; do
    if [[ -s "\${rmd}" ]]; then
        R < "\${rmd}\" --no-save
    fi
done

ls -lahtr