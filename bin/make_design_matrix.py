#!/usr/bin/env python3

import os
import pandas as pd
import sys

sample_names = sys.argv[1]
assert os.path.exists(sample_names)
control_names = sys.argv[2]
assert os.path.exists(control_names)

def read_list(fp):
    return [
        l.rstrip()
        for l in open(fp, 'r')
    ]

sample_names = read_list(sample_names)
control_names = read_list(control_names)

pd.DataFrame([
    {
        "Samples": n,
        "baseline": 1,
        "treatment_control": n in control_names
    }
    for n in sample_names
]).reindex(
    columns=[
        "Samples", "baseline", "treatment_control"
    ]
).to_csv(
    "design_matrix.tsv",
    sep="\t",
    index=None
)