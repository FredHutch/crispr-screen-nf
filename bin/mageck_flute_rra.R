#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# Set up the paths to the input files
gene_summary_txt = args[1]
file1 = file.path(system.file("extdata", package = "MAGeCKFlute"), gene_summary_txt)
print(gene_summary_txt)
print(file.exists(gene_summary_txt))

sgrna_summary_txt = args[2]
file2 = file.path(system.file("extdata", package = "MAGeCKFlute"), sgrna_summary_txt)
print(sgrna_summary_txt)
print(file.exists(sgrna_summary_txt))

# Project name
proj = args[3]

# Organism
organism = args[4]

# Scale cutoff
scale_cutoff = args[5]

# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
FluteRRA(file1, file2, proj=proj, organism=organism, scale_cutoff = scale_cutoff, outdir = getwd())
