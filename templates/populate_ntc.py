#!/usr/bin/env python3

import os
import pandas as pd


def read_csv(fp, sep=',', required_columns=['sgRNA']):
    """Read in a table and make sure that there is a column named sgRNA."""

    assert os.path.exists(fp)
    df = pd.read_csv(fp, sep=sep)

    # Account for the fact that sometimes the guide library
    # definition file has the sgRNA column in lowercase
    df = df.rename(
        columns=dict(
                sgrna='sgRNA',
                guide='Guide',
                gene='Gene'
            )
        )

    # Make sure that there is both a Guide and sgRNA column
    for cname in required_columns:
        msg = f"Could not find {cname} in columns ({', '.join(df.columns.values)})"
        assert cname in df.columns.values, msg

    return df

# Counts table
counts = read_csv("input.counts.txt", sep='\\t', required_columns=['sgRNA', 'Gene'])

# Library definition file
library = read_csv("${library}", sep=',', required_columns=['sgRNA', 'Guide'])

# Iterate over each of the guides in the library which share a sequence
for guide_seq, shared_guides in library.groupby('sgRNA'):

    # If there is more than one
    if shared_guides.shape[0] > 1:

        # Get the ID of the guide which does have valid counts
        guide_counts = counts.set_index('sgRNA').reindex(index=shared_guides.Guide.values).dropna()

        # There should be no more than one set of counts
        msg = f"Found multiple rows of counts from mageck for the identical guide {guide_seq}"
        assert guide_counts.shape[0] <= 1, shared_guides

        # If there are no counts
        if guide_counts.shape[0] == 0:

            # Skip it
            continue

        # Drop the guide and gene columns from the counts, convert to a dict
        guide_counts = guide_counts.drop(columns=['Gene']).iloc[0].to_dict()
        
        # Add rows to the counts table for each of the guides which were omitted by mageck count
        counts = pd.concat(
            [
                counts,
                pd.DataFrame([
                    {
                        'sgRNA': r.Guide,
                        'Gene': r.Gene,
                        **guide_counts
                    }
                    for _, r in shared_guides.iterrows()
                    if r.Guide not in counts.sgRNA.values
                ])
            ]
        )

# Make sure that all counts are numeric
counts = counts.astype(
    {
        cname: 'int'
        for cname in counts.columns.values
        # Skipping the columns with gene and guide names
        if cname not in ['Gene', 'sgRNA']
    }
)

# Sort by gene and guide
counts = counts.sort_values(by=['Gene', 'sgRNA'])

# Remove any guides which had 0 counts across all samples
zero_counts = counts.drop(
    columns=['Gene', 'sgRNA']
).sum(
    axis=1
).apply(
    lambda x: x == 0
)
if zero_counts.any():
    print(f"Removing {zero_counts.sum():,} guides with 0 counts")
    counts = counts.drop(
        index=counts.index.values[zero_counts]
    )

# Write out the counts table
counts.to_csv('counts.txt', sep='\t', index=None)
