#!/usr/bin/env python3

import os
import pandas as pd

# Function to read in a set of TSV files and join by shared index columns
def read_tsv_list(fp_list, index_cols):

    return pd.concat(
        [
            pd.read_csv(
                fp,
                sep = "\t"
            ).set_index(
                index_cols
            )
            for fp in fp_list
        ],
        axis=1
    )

# Function to list filepaths in a folder
def list_files_in_folder(folder):
    return [
        os.path.join(folder, fp)
        for fp in os.listdir(folder)
    ]

# Shared index columns
index_cols = ['sgRNA', 'Gene']

# Read in all of the treatment data
count_data = read_tsv_list(
    list_files_in_folder("counts"),
    index_cols
)

# Write out a list of the specimens in the treatment list
with open("sample_names.txt", "w") as handle:
    handle.write(",".join(count_data.columns.values))

# Combine the counts and write out to a file
count_data.to_csv(
    f"counts.txt",
    sep="\t"
)

# Make a summary table
summary_data = pd.DataFrame([
    dict(
        sample=sample,
        n_guides=sample_counts.shape[0],
        n_detected=(sample_counts > 0).sum(),
        n_undetected=(sample_counts == 0).sum(),
        n_reads=sample_counts.sum(),
        mean_depth=sample_counts.mean(),
        median_depth=sample_counts.median(),
        n_ten_reads=(sample_counts >= 10).sum(),
        n_hundred_reads=(sample_counts >= 100).sum()
    )
    for sample, sample_counts in count_data.items()
])

# Write out the summary table to a file
summary_data.to_csv(
    "summary.txt",
    sep="\t",
    index=None
)