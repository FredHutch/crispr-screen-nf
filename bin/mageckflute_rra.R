#!/usr/bin/env Rscript

requireNamespace("ggplot2")
library(MAGeCKFlute)
library(ggplot2)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)

# Set up the paths to the input files
gene_summary_txt = args[1]
print(gene_summary_txt)
print(file.exists(gene_summary_txt))

sgrna_summary_txt = args[2]
print(sgrna_summary_txt)
print(file.exists(sgrna_summary_txt))

# Organism
organism = args[3]

# Scale cutoff
scale_cutoff = args[4]

# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
FluteRRA(
    gene_summary_txt,
    sgrna_summary_txt,
    proj="rra",
    organism=organism,
    scale_cutoff = as.numeric(scale_cutoff),
    outdir = getwd(),
    incorporateDepmap = FALSE
#     omitEssential = FALSE,
#     type = "KEGG+REACTOME+GOBP"
)
