#!/usr/bin/env Rscript
#################################################################################################
# R script to make diagnostic plots for DGE analysis
# #
# RNA-seq data set with 33 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesion with normal tissue, and 16 adenoma lesion with normal tissue
# 
# Stephany Orjuela, January 2018
################################################################################################

library(edgeR)
library(Rsubread)

files <- as.matrix(read.table("new_exp_setup.csv"))

fc_gene <- featureCounts(files[,1], annot.ext="Homo_sapiens.GRCh37.gtf", 
                                   nthreads=1, isPairedEnd=T, useMetaFeatures=T, isGTFAnnotationFile=T)

#Create DGElist object
x <- DGEList(counts=fc_gene$counts, samples = files[,2], group = files[,3])
x <- calcNormFactors(x)
keep <- rowSums(cpm(x) > 1) >= 2
x <- x[keep,]
colnames(x$counts) = x$samples$samples

save(x, file="DGEobj_fc.RData")

cm <- x$counts
nf <-  x$samples$norm.factors
ave.libsize <- exp(mean(log(colSums(cm)))) 
nf <- (ave.libsize/(nf * colSums(cm)))
names(nf) <- paste0(x$samples$samples,".",x$samples$group)
files <- read.csv("bam_files.csv", header = F)$V1

for (i in 1:length(nf)) {
  cmd1 <- paste0("bedtools genomecov -split -scale ", nf[i], " -ibam ", files[i], " -bga > bedGraphs/", names(nf)[i], ".bedgraph")
  cmd2 <- paste0("bedGraphToBigWig bedGraphs/", names(nf)[i], ".bedgraph", " temp.hg19.chrom.sizes bigwigs/", names(nf)[i], ".bw")
  cat(cmd1, "\n")
  system(cmd1)
  cat(cmd2, "\n")
  system(cmd2)
}

