#!/usr/bin/env Rscript

library(MAGeCKFlute)

args = commandArgs(trailingOnly=TRUE)

# Set up the paths to the input files
gene_summary_txt = args[1]
print(gene_summary_txt)
print(file.exists(gene_summary_txt))

sgrna_summary_txt = args[2]
print(sgrna_summary_txt)
print(file.exists(sgrna_summary_txt))

# Project name
proj = args[3]

# Organism
organism = args[4]

# Scale cutoff
scale_cutoff = args[5]

# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
FluteRRA(gene_summary_txt, sgrna_summary_txt, proj=proj, organism=organism, scale_cutoff = as.numeric(scale_cutoff), outdir = getwd())
