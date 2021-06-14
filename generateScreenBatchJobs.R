options(stringsAsFactors=F)
sampleSheet<-read.csv("/fh/fast/_SR/Genomics/ngs/illumina/solexa/SampleSheets/190812_D00300_0802_BH3FV2BCX3_lcarter.csv")

names(sampleSheet)[3] <- "Sample"

#TO DO: add some sanity checks
#fastx toolkit requires that the sample names be alphanumeric
sampleSheet$Sample<-gsub("-","_",sampleSheet$Sample)
sampleSheet$Sample<-gsub(" ","_",sampleSheet$Sample)
#moabAcct<-"paddison_p"

sampleSheet <- data.frame(Sample = sampleSheet$Sample)
sampleSheet$User <- "lcarter"
sampleSheet$FlowCell <- "190812_D00300_0802_BH3FV2BCX3"
sampleSheet <- unique(sampleSheet); nrow(sampleSheet) #31
#associate each sample with the correct reference
sampleSheet$Reference <- "Kinetochore"
sampleSheet$Reference[grep("CPD", sampleSheet$Sample)] <- "CPD"
sampleSheet$Reference[grep("Brunello", sampleSheet$Sample)] <- "Brunello"
sampleSheet$Reference[grep("Spliceosome", sampleSheet$Sample)] <- "Spliceosome"


for(i in 1:nrow(sampleSheet)){
  cat("#!/bin/bash", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n")
  cat("#SBATCH -N1 -n1 -t 0-4 -p campus --mail-type=END --mail-user=pchanana@fhcrc.org ", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat("PATH=/home/solexa/apps/fastx_toolkit_0.0.13:/app/cutadapt/1.1/bin:/home/solexa/apps/bowtie/bowtie-1.0.0:$PATH", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  #demultiplex
  cat("runDir=/fh/fast/_SR/Genomics/user/pchanana/2019.08.15.lcarter/align", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat("fastqDir=$runDir/fastq", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat("cd $fastqDir", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat(paste("sampleName=",sampleSheet$Sample[i],sep=""), file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  #align
  if(sampleSheet$Reference[i] == "Brunello"){
    cat("bowtieGenome=/shared/solexa/solexa/Genomes/genomes/sgRNA/Brunello/Brunello.fa", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  }
  if(sampleSheet$Reference[i] == "CPD"){
    cat("bowtieGenome=/shared/solexa/solexa/Genomes/genomes/sgRNA/CPD/CPD.fa", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  }
  if(sampleSheet$Reference[i] == "Kinetochore"){
    cat("bowtieGenome=/shared/solexa/solexa/Genomes/genomes/sgRNA/Kinetochore/Kinetochore.fa", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  }
  if(sampleSheet$Reference[i] == "Spliceosome"){
    cat("bowtieGenome=/shared/solexa/solexa/Genomes/genomes/sgRNA/Spliceosome/Spliceosome.fa", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  }
  cat("bowtieDir=$runDir/bowtie", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat("mkdir -p $bowtieDir", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat(paste("zcat $sampleName.fastq.gz | bowtie -p 1 --trim5 30 --trim3 0 -n 0 $bowtieGenome - $bowtieDir/$sampleName.bt",sep=""), file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
  cat("exit 0", file=paste(sampleSheet$Sample[i],".sbatch",sep=""),sep="\n", append=TRUE)
}

#system("mv *.sbatch ~/Desktop/ngs/ngs/illumina/lcarter/190812_D00300_0802_BH3FV2BCX3/qsub_files")


