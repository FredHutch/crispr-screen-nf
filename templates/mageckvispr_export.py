#!/usr/bin/env python3

import json
import logging
import os
from vispr import Screen
import yaml

# The input files which should be present in the working directory
# are count_normalized.txt and gene_summary.txt
assert os.path.exists("count_normalized.txt")
assert os.path.exists("gene_summary.txt")

# Make a YAML file informing VISPR where those files are
configpath = "config.yaml"
with open(configpath, "w") as handle:
    yaml.dump(
        dict(
            experiment="myexperiment",
            species="homo_sapiens",
            assembly="hg19",
            targets=dict(results="gene_summary.txt", genes=True),
            sgrnas=dict(counts="count_normalized.txt")
        ),
        handle
    )

# Use that config file to load the data
with open(configpath) as f:
    screen = Screen(yaml.safe_load(f), parentdir="./")

def write(dat, name):
    with open(name + ".json", "w") as out:
        json.dump(json.loads(dat)['data'], out, indent=2)
        logging.info(f"Error writing out data for {name}")

if screen.fastqc is not None:
    write(screen.fastqc.plot_gc_content(), "gc-content")
    write(screen.fastqc.plot_base_quality(), "base-quality")
    write(screen.fastqc.plot_seq_quality(), "seq-quality")
if screen.mapstats is not None:
    write(screen.mapstats.plot_mapstats(), "mapped-unmapped")
    write(screen.mapstats.plot_zerocounts(), "zerocounts")
    write(screen.mapstats.plot_gini_index(), "gini-index")
write(screen.rnas.plot_normalization(), "readcounts")
write(screen.rnas.plot_readcount_cdf(), "readcount-cdf")
write(screen.rnas.plot_correlation(), "correlation")
write(screen.rnas.plot_pca(1, 2), "pca-1-2")
write(screen.rnas.plot_pca(1, 3), "pca-1-3")
write(screen.rnas.plot_pca(2, 3), "pca-2-3")

for condition, results in screen.targets.items():
    for selection, results in results.items():
        pre = ".".join(([] if condition == "default" else [condition]) +
                        [selection.replace(" ", "-")])
        try:
            write(results.plot_pvals(), pre + ".p-values")
        except:
            logging.info("Error writing out p-values")
        try:
            write(results.plot_pval_hist(), pre + ".p-value-hist")
        except:
            logging.info("Error writing out log10-pvalues")
