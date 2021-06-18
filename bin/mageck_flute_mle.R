#!/usr/bin/env Rscript

library(MAGeCKFlute)

args = commandArgs(trailingOnly=TRUE)

# Set up the paths to the input files
gene_summary_txt = args[1]
print(gene_summary_txt)
print(file.exists(gene_summary_txt))

# Project name
proj = args[2]

# Organism
organism = args[3]

# Name of treatment group
treatname = args[4]

# Name of control group
ctrlname = args[5]

# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
FluteMLE(gene_summary_txt, proj=proj, organism=organism, treatname=treatname, ctrlname=ctrlname)
