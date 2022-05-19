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
