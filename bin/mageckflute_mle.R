#!/usr/bin/env Rscript

library(MAGeCKFlute)
library(ggplot2)
library(stats) 
library(grDevices)
library(utils)
library(gridExtra)
library(grid)

args = commandArgs(trailingOnly=TRUE)

# Set up the paths to the input files
gene_summary_txt = args[1]
print(gene_summary_txt)
print(file.exists(gene_summary_txt))

# Organism
organism = args[2]

# Name of treatment group
treatname = args[3]

# Name of control group
ctrlname = args[4]

# Depmap Effect
depmapEffect = args[5]

# Depmap Samples
depmapSamples = args[6]

proj = ""

IncorporateDepmap <- function(dd, symbol = "id",
                              cell_lines = NA, lineages = "All",
                              na.rm = FALSE){
  ## Load Depmap data
  depmap_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_19Q3.rds")
  
  if(file.exists(depmap_rds)){
    Depmap_19Q3 = readRDS(depmap_rds)
  }else{  
    Depmap_19Q3 = t(read.csv(depmapEffect, header = TRUE,
                             row.names = 1, stringsAsFactors = FALSE, check.names = FALSE))
    rownames(Depmap_19Q3) = gsub(" .*", "", rownames(Depmap_19Q3))
    saveRDS(Depmap_19Q3, depmap_rds)
  }    
  meta_rds = file.path(system.file("extdata", package = "MAGeCKFlute"), "Depmap_sample_info.rds")
  if(file.exists(meta_rds)){   
    sampleinfo = readRDS(meta_rds)
  }else{
    sampleinfo = read.csv(depmapSamples,
                          row.names = 1, header = TRUE, stringsAsFactors = FALSE)                  
    saveRDS(sampleinfo, meta_rds)
  }
  sampleinfo = sampleinfo[colnames(Depmap_19Q3), ]
  idx1 = sampleinfo$lineage%in%lineages
  idx2 = sampleinfo$stripped_cell_line_name%in%cell_lines |
    sampleinfo$CCLE.Name%in%cell_lines |
    sampleinfo$alias%in%cell_lines
  genes = as.character(dd[, symbol])
  Depmap_19Q3 = as.data.frame(Depmap_19Q3, stringsAsFactors = FALSE)
  if(sum(idx2)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, idx2], na.rm = TRUE))
  }else if(sum(idx1)>0){
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, idx1], na.rm = TRUE))
  }else{
    dd = cbind(dd, Depmap = rowMeans(Depmap_19Q3[genes, ], na.rm = TRUE))
  }
  if(na.rm) dd = na.omit(dd)
  return(dd)
}

FluteMLE <- function(gene_summary, treatname, ctrlname = "Depmap",
                     keytype = "Symbol", organism = "hsa", # Input dataset
                     incorporateDepmap = TRUE,
                     cell_lines = NA, lineages = "All",
                     norm_method = "cell_cycle", posControl = NULL,
                     omitEssential = FALSE,
                     top = 10, toplabels = NA,
                     scale_cutoff = 2, limit = c(0,200),
                     pvalueCutoff=0.25, enrich_method = "ORT", proj = NA,
                     width = 10, height = 7, outdir = ".",
                     pathview.top = 4, verbose = TRUE){
	## Prepare the running environment ##
  {
    message(Sys.time(), " # Create output dir and pdf file...")
    outdir = file.path(outdir)
    dir.create(file.path(outdir), showWarnings = FALSE)
    output_pdf = paste0( norm_method, ".pdf")
    pdf(file.path(outdir, output_pdf), width = width, height = height)
  }

  ## Beta Score Preparation ##
  {
    beta = ReadBeta(gene_summary)
    if(tolower(keytype)!="entrez"){
      beta$EntrezID = TransGeneID(beta$Gene, keytype, "Entrez", organism = organism)
    }else{
      beta$EntrezID = beta$Gene
    }
    if(tolower(keytype)!="symbol"){
      beta$Symbol = TransGeneID(beta$Gene, keytype, "Symbol", organism = organism)
    }else{
      beta$Symbol = beta$Gene
    }
    if(organism != "hsa"){
      beta$HumanGene = TransGeneID(beta$Symbol, "Symbol", "Symbol",
                                   fromOrg = organism, toOrg = "hsa")
    }else{
      beta$HumanGene = beta$Symbol
    }

    message(Sys.time(), " # Transform id to official human gene name ...")
    idx1 = is.na(beta$EntrezID)
    idx2 = !is.na(beta$EntrezID) & duplicated(beta$EntrezID)
    idx = idx1|idx2
    if(sum(idx2)>0) message(sum(idx2), " genes have duplicate Entrez IDs: ",
                            paste0(beta$Gene[idx2], collapse = ", "))

    dd = beta[!idx, ]
    if(incorporateDepmap | "Depmap"%in%ctrlname)
      dd = IncorporateDepmap(dd, symbol = "HumanGene", cell_lines = cell_lines,
                             lineages = lineages)
    if(!all(c(ctrlname, treatname) %in% colnames(dd)))
      stop("Sample name doesn't match !!!")
    ## Normalization
    if(tolower(norm_method)=="cell_cycle")
      dd = NormalizeBeta(dd, id = "HumanGene", samples = c(ctrlname, treatname),
                         method = "cell_cycle", posControl = posControl)
    if(tolower(norm_method)=="loess")
      dd = NormalizeBeta(dd, id = "HumanGene", samples = c(ctrlname, treatname), method = "loess")
  }
  message("137")
	## Distribution of beta scores ##
	{
	  ## All genes ##
	  outputDir1 = file.path(outdir, "QC/")
	  dir.create(outputDir1, showWarnings = FALSE)
	  idx_distr = c(ctrlname, treatname)
	  P1 = ViolinView(dd[, idx_distr], ylab = "Beta score", main = "All genes",
	                  filename = paste0(outputDir1, "ViolinView_all_", norm_method, ".png"))
	  P2 = DensityView(dd[, idx_distr], xlab = "Beta score", main = "All genes",
	                   filename = paste0(outputDir1, "DensityView_all_", norm_method, ".png"))		   
	  P3 = ConsistencyView(dd, ctrlname, treatname, main = "All genes",
	                       filename = paste0(outputDir1, "Consistency_all_", norm_method, ".png"))
	  P4 = MAView(dd, ctrlname, treatname, main = "All genes",
	              filename = paste0(outputDir1, "MAView_all_", norm_method, ".png"))
	write.table(dd, paste0(outputDir1, "beta.all.genes.dist.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
	  ## Essential genes ##
	  Zuber_Essential = readRDS("Zuber_Essential.rds")

	  if(is.null(posControl))
	    idx = toupper(dd$HumanGene) %in% toupper(Zuber_Essential$GeneSymbol)
	  else
	    idx = toupper(dd$Gene) %in% toupper(posControl)

	  if(sum(idx) > 10){
	    P1 = ViolinView(dd[idx, idx_distr], ylab = "Essential.B.S.", main = "Essential genes",
	                    filename = paste0(outputDir1, "ViolinView_posctrl_", norm_method, ".png"))
	    P2 = DensityView(dd[idx, idx_distr], xlab = "Essential.B.S.", main = "Essential genes",
	                     filename = paste0(outputDir1, "DensityView_posctrl_", norm_method, ".png"))
	    P3 = ConsistencyView(dd[idx, ], ctrlname, treatname, main = "Essential genes",
	                         filename = paste0(outputDir1, "Consistency_posctrl_", norm_method, ".png"))
	    P4 = MAView(dd[idx, ], ctrlname, treatname, main = "Essential genes",
	                filename = paste0(outputDir1, "MAView_posctrl_", norm_method, ".png"))
	    grid.arrange(P1, P2, P3, P4, ncol = 2)
		write.table(dd, paste0(outputDir1, "beta.essential.genes.dist.txt"),
                sep = "\t", row.names = FALSE, quote = FALSE)
	  }
	}
	message("174")
  ## Combine replicates ##
  {
    dd$Control = Biobase::rowMedians(as.matrix(dd[, ctrlname, drop = FALSE]))
    dd$Treatment = Biobase::rowMedians(as.matrix(dd[, treatname, drop = FALSE]))
    dd$Diff = dd$Treatment - dd$Control
    x_cut = c(-CutoffCalling(dd$Control, scale_cutoff),
              CutoffCalling(dd$Control, scale_cutoff))
    y_cut = c(-CutoffCalling(dd$Treatment, scale_cutoff),
              CutoffCalling(dd$Treatment, scale_cutoff))
    intercept = c(-CutoffCalling(dd$Diff, scale_cutoff),
                  CutoffCalling(dd$Diff, scale_cutoff))
    if(omitEssential){
      dd = OmitCommonEssential(dd, symbol = "HumanGene", dependency = -0.4)
    }
    dd$Rank = rank(dd$Diff)
    dd$RandomIndex = sample(1:nrow(dd), nrow(dd))
    dd$Gene = dd$Symbol
	write.table(dd, paste0(outputDir1, "omitessential_data.txt"),
			sep = "\t", row.names = FALSE, quote = FALSE)
  }
message("197")
	## Drug-targeted genes ##
	{
	  outputDir2 = file.path(outdir, "Selection/")
	  dir.create(outputDir2, showWarnings = FALSE)

	  p1 = ScatterView(dd, "Control", "Treatment", groups = c("top", "bottom"),
	                   groupnames = c("GroupA", "GroupB"), intercept = intercept)
	  ggsave(paste0(outputDir2, "ScatterView_TreatvsCtrl_", norm_method, ".png"),
	         p1, width = 4, height = 3)
	  write.table(p1$data, paste0(outputDir2, "Data_ScatterView_TreatvsCtrl.txt"),
	              sep = "\t", row.names = FALSE, quote = FALSE)
	  p2 = ScatterView(dd, x = "Rank", y = "Diff", label = "Symbol",
	                   groups = c("top", "bottom"), groupnames = c("GroupA", "GroupB"),
	                   top = top, y_cut = y_cut)
	  ggsave(paste0(outputDir2, "RankView_Treat-Ctrl_", norm_method, ".png"),
	         p2, width = 3, height = 5)
	  p3 = ScatterView(dd[dd$Diff>0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
	                   y_cut = y_cut, groups = "top", groupnames = c("GroupA"), top = top)
	  ggsave(paste0(outputDir2, "ScatterView_Treat-Ctrl_Positive_", norm_method, ".png"),
	         p3, width = 4, height = 3)
	  p4 = ScatterView(dd[dd$Diff<0, ], x = "RandomIndex", y = "Diff", label = "Symbol",
	                   y_cut = y_cut, groups = "bottom", groupnames = c("GroupB"), top = top)
	  ggsave(paste0(outputDir2, "ScatterView_Treat-Ctrl_Negative_", norm_method, ".png"),
	         p4, width = 4, height = 3)

	  grid.arrange(p1, p2, p3, p4, ncol = 2)
	}
	
	## Enrichment analysis of negative and positive selected genes ##
	{
	  outputDir3 = file.path(outdir, "Enrichment/")
	  outputDir4 = file.path(outdir, "PathwayView/")
	  dir.create(outputDir3, showWarnings=FALSE)
	  dir.create(outputDir4, showWarnings=FALSE)

	  E1 = EnrichAB(p1$data, pvalue = pvalueCutoff, enrich_method = enrich_method,
	                organism = organism, limit = limit,
	                filename = norm_method, out.dir = outputDir3)
	  # EnrichedView
	  grid.arrange(E1$keggA$gridPlot, E1$reactomeA$gridPlot, E1$gobpA$gridPlot, E1$complexA$gridPlot, ncol = 2)
	  grid.arrange(E1$keggB$gridPlot, E1$reactomeB$gridPlot, E1$gobpB$gridPlot, E1$complexB$gridPlot, ncol = 2)

	  # Pathway view for top 4 pathway
	  if(!is.null(E1$keggA$enrichRes) && nrow(E1$keggA$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$keggA$enrichRes$ID),
	                    top = pathview.top, ncol = 2, title="Group A",
	                    organism=organism, output=outputDir4)
	  if(!is.null(E1$keggB$enrichRes) && nrow(E1$keggB$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$keggB$enrichRes$ID),
	                    top = pathview.top, ncol = 2,
	                    title="Group B", organism=organism,
	                    output=outputDir4)
		
  }

	## Nine-squares ##
	{
	  p1 = ScatterView(dd, x = "Control", y = "Treatment", label = "Symbol",
	                   groups = c("midleft", "topcenter", "midright", "bottomcenter"),
	                   groupnames = c("Group1", "Group2", "Group3", "Group4"),
	                   top = top, display_cut = TRUE,
	                   x_cut = x_cut, y_cut = y_cut, intercept = intercept)
	  ggsave(paste0(outputDir2, "SquareView_", norm_method, ".png"), p1, width = 5, height = 4)
	  write.table(p1$data, paste0(outputDir2, proj, "squareview_data.txt"),
	              sep = "\t", row.names = FALSE, quote = FALSE)
	  grid.arrange(p1, ncol = 1)
	}

	## Nine-Square grouped gene enrichment ##
	{
	  E1 = EnrichSquare(p1$data, id = "EntrezID", keytype = "entrez",
	                    x = "Control", y = "Treatment", organism=organism,
	                    pvalue = pvalueCutoff, enrich_method = enrich_method,
	                    filename=norm_method, limit = limit, out.dir=outputDir3)
    # EnrichView
	  grid.arrange(E1$kegg1$gridPlot, E1$reactome1$gridPlot, E1$gobp1$gridPlot, E1$complex1$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg2$gridPlot, E1$reactome2$gridPlot, E1$gobp2$gridPlot, E1$complex2$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg3$gridPlot, E1$reactome3$gridPlot, E1$gobp3$gridPlot, E1$complex3$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg4$gridPlot, E1$reactome4$gridPlot, E1$gobp4$gridPlot, E1$complex4$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg12$gridPlot, E1$reactome12$gridPlot, E1$gobp12$gridPlot, E1$complex12$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg13$gridPlot, E1$reactome13$gridPlot, E1$gobp13$gridPlot, E1$complex13$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg24$gridPlot, E1$reactome24$gridPlot, E1$gobp24$gridPlot, E1$complex24$gridPlot, ncol = 2)
	  grid.arrange(E1$kegg34$gridPlot, E1$reactome34$gridPlot, E1$gobp34$gridPlot, E1$complex34$gridPlot, ncol = 2)

	  # PathwayView
	  if(!is.null(E1$kegg1$enrichRes) && nrow(E1$kegg1$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg1$enrichRes$ID), ncol = 2,
	                    top = pathview.top, title = "Group 1", organism=organism,
	                    output=outputDir4)
	  if(!is.null(E1$kegg2$enrichRes) && nrow(E1$kegg2$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg2$enrichRes$ID), ncol = 2,
	                    top = pathview.top, title = "Group 2",
	                    organism=organism,output=outputDir4)
	  if(!is.null(E1$kegg3$enrichRes) && nrow(E1$kegg3$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg3$enrichRes$ID), ncol = 2,
	                    top = pathview.top, title = "Group 3",
	                    organism=organism, output=outputDir4)
	  if(!is.null(E1$kegg4$enrichRes) && nrow(E1$kegg4$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg4$enrichRes$ID), ncol = 2,
	                    title = "Group 4", organism = organism,
	                    top = pathview.top, output=outputDir4)
	  if(!is.null(E1$kegg12$enrichRes) && nrow(E1$kegg12$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg12$enrichRes$ID), ncol = 2,
	                    title = "Group 1 & Group 2", organism=organism,
	                    top = pathview.top, output=outputDir4)
	  if(!is.null(E1$kegg13$enrichRes) && nrow(E1$kegg13$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg13$enrichRes$ID), ncol = 2,
	                    title = "Group 1 & Group 3", organism=organism,
	                    top = pathview.top, output=outputDir4)
	  if(!is.null(E1$kegg24$enrichRes) && nrow(E1$kegg24$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg24$enrichRes$ID), ncol = 2,
	                    title = "Group 2 & Group 4", organism=organism,
	                    top = pathview.top, output=outputDir4)
	  if(!is.null(E1$kegg34$enrichRes) && nrow(E1$kegg34$enrichRes)>0)
	    arrangePathview(dd, gsub("KEGG_", "", E1$kegg34$enrichRes$ID), ncol = 2,
	                    title = "Group 3 & Group 4", organism=organism,
	                    top = pathview.top, output=outputDir4)
  }
	dev.off()
}


# Run the MAGeCK Flute pipeline for the output of MAGeCK test RRA
try(FluteMLE(gene_summary_txt, proj=proj, organism=organism, treatname=treatname, keytype="Symbol"))
