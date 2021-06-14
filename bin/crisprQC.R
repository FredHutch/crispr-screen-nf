args =  commandArgs(TRUE)
countsDir=args[1]
outFile=args[2]
group<-args[3]
if (file.exists(countsDir) == FALSE || length(args) != 3){
  writeLines ("Usage:\nRscript crisprQC.R countsDir outFileBaseName groups<colon-sep>\n");
  writeLines ("outFileBaseName should be the same as the basename for the counts file used as input\n");
  writeLines ("groups is colon separated list of group corresponding to each sample\n");
  quit()
}

#install.packages("ggplot2")
#install.packages("ggrepel")

library(edgeR)
library(ggplot2)
library(ggrepel)
library(reshape2)
options(stringsAsFactors = F)

setwd(countsDir)
d<-read.delim(file=paste0(countsDir,outFile,"_counts.txt"), stringsAsFactors=F)
d<-data.frame(d[2:ncol(d)], row.names=d$guideID)
d[is.na(d)]<-0

targets <- data.frame(Sample = names(d),
                      Group = names(d))
#targets$Group <- sub("[0-3]$", "", targets$Group)
#targets<-targets[order(targets$Group, targets$Sample),]
#group<-targets$Group
group <- strsplit(group,":")[[1]]
d<-d[,as.character(targets$Sample)]
d<-DGEList(counts=d,group=group)

d <- calcNormFactors(d)
nd <- data.frame(cpm(d, normalized.lib.sizes = T, log = T), check.names = F)

# nd <- data.frame(cpm(d, normalized.lib.sizes = F, log = T), check.names = F)

###boxplots
png(file=paste0(outFile, "_boxPlots.png"))
par(mar = c(12,5,3,3))
boxplot(nd, las = 3, main = outFile, ylab = "normalized log2(CPM)")
dev.off()


##correlation matrix
nd.cor <- data.frame(cpm(d, normalized.lib.sizes = T, log = T), check.names = F, row.names=NULL)
get_upper_tri <- function(cormat) {
   cormat[lower.tri(cormat)] <- NA
   return(cormat)
}
cormat=round(cor(nd.cor,method="spearman"),2)
upper_tri <- get_upper_tri(cormat)
#lower_tri
#head(cormat)
melted_cormat=melt(upper_tri, na.rm = TRUE)
#melted_cormat

heat <- ggplot(data=melted_cormat, aes(Var2, Var1, fill=value))+ geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white" , midpoint = 0, limit=c(-1,1), space="Lab",name="Spearman Correlation") + coord_fixed() + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background=element_blank(), panel.border=element_blank())
heat+geom_text(aes(Var2, Var1, label=value), color = "black", size=1.5) ##to print correlation vals on plot
ggsave(paste0(outFile, "_corr.pdf"),width=10, height=10, units="in", dpi=800)

##PCA plot
#sds <- apply(nd, 1, sd)
#sds <- sds[order(sds, decreasing = T)]
#sds.top <- sds[1:10000]
#nd.top <- nd[row.names(nd) %in% names(sds.top),]

#nd <- data.frame(cpm(d, normalized.lib.sizes = T, log = F), check.names = F)
pca_data <- prcomp(t(nd))
pca_data_perc <- round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)
df_pca_data <- data.frame(PC1=pca_data$x[,1], PC2=pca_data$x[,2], sample=colnames(nd), condition=group)
head(df_pca_data)

ggplot(df_pca_data,aes(PC1,PC2, color=condition)) + geom_point(size=4) + labs(x=paste0("PC1(",pca_data_perc[1],")"), y=paste0("PC2(",pca_data_perc[2],")")) + geom_text_repel(data=df_pca_data, aes(PC1, PC2, label=row.names(df_pca_data)))
#+ geom_text_repel(aes(label=rownames(data)))
ggsave(paste0(outFile, "_PCA.pdf"),width=10, height=10, units="in", dpi=800)
