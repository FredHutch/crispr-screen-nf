#!/bin/bash

set -euo pipefail

parse_fastqc_data.py "${fastqc_txt}" "${fastqc_txt.name.replaceAll(/.txt/, '')}.json"
