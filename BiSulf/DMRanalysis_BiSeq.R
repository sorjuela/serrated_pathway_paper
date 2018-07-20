#!/usr/bin/env Rscript
# run as [R < scriptName.R --no-save]

#####################################################################################################
# R script to run a DMR analysis workflow using BiSeq starting from .cov files obtained from Bismark
#
# BS-seq data set with 32 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
# 
# Most of code taken from BiSeq vignette
#
# Stephany Orjuela, January 2018
####################################################################################################

library(BiSeq)

### Upload metadata ###-----------------------------------------------------------------------------

files <- read.table("design_matrix.csv", stringsAsFactors = F) #prelesions
#files <- read.table("CRCmatrix.csv", stringsAsFactors = F) #lesions

### Read cov files into R ###-----------------------------------------------------------------------
# These files contain all the covered CpG sites from all the samples

BSr <- readBismark(files$V2, colData= DataFrame(group = files$V3, row.names = paste0(files$V1,".",files$V3))) 
#save(BSr, file="BSraw.RData") 

#Since BiSeq takes a long time to run, we speed this up by looping through each chromosome, 
#excluding X, Y and MT. The raw and smoothed count objects are saved for further plotting

for(i in 1:22){ 

load(file="BSraw.RData")
w <- which(seqnames(BSr) == paste0("chr",i))
BSrchr <- BSr[w, ]
rm(BSr)
save(BSrchr, file = paste0("chr",i,"BSraw.RData"))

### Define clusters ###----------------------------------------------------------------------------- 

BSr.clust.unlim <- clusterSites(BSrchr,
                                groups = factor(colData(BSrchr)$group),
                                perc.samples = 0.8,
                                min.sites = 3,
                                max.dist =  500,
                                mc.cores = 30)

### smooth the methylation values of CpG sites within the clusters ###----------------------------- 
predictedMeth <- predictMeth(object = BSr.clust.unlim, h=1000, mc.cores = 30)
#save(predictedMeth, file=paste0("predictedMethchr",i,".RData"))

### Model methylation within a beta regression for each comparison ###-----------------------------

#subset for groups, change depending on group of interest
#SSA
s <- which(colData(predictedMeth)$group == "SSA")
n <- s-1

#Adenoma
#a <- which(colData(predictedMeth)$group == "ADENOMA")
#n <- a + 1

SNpMeth <- predictedMeth[,c(s,n)]
colData(SNpMeth)$group <- factor(colData(SNpMeth)$group, levels = c("ADENOMA", "NORMAL"))

betaResultsSN <- betaRegression(formula = ~group, link = "probit", object = SNpMeth, type = "BR", mc.cores = 30) 


### Test CpG clusters for differential methylation ###----------------------------------------------

#1st: Transform pvals into Zscores, which are Normal under the Ho:no group effect
#model beta regression again for resampled data
predictedMethNullSN <- SNpMeth[,c(1:28)] #Subset original table if desired 
colData(predictedMethNullSN)$group.null <- rep(c(1,2), 14) #generate null group
betaResultsNullSN <- betaRegression(formula = ~group.null, link = "probit", object = predictedMethNullSN, type="BR", mc.cores = 30)

#2nd: estimate the variogram for the Z scores obtained from the resampled data
varioSN <- makeVariogram(betaResultsNullSN)
vario.sm <- smoothVariogram(varioSN, sill = 0.9)

#3rd: replace the pValsList object (which consists of the test results of the
#resampled data) by the test results of interest (for group effect)
vario.aux <- makeVariogram(betaResultsSN, make.variogram=FALSE)
vario.sm$pValsList <- vario.aux$pValsList

#4th: The correlation of the Z scores between two locations in a cluster can now be estimated: 
locCor <- estLocCor(vario.sm)

#5th: test each CpG cluster for the presence of at least one differentially methylated 
#location at q what can be interpreted as the size-weighted FDR on clusters:
clusters.rej <- testClusters(locCor, FDR.cluster = 0.05) 

#6th: Trim the rejected CpG clusters, that is to remove the not differentially methylated CpG sites at 
#q 1, what can be interpreted as the location-wise FDR:
#This is a table with all the diff sites
clusters.trimmedSN <- trimClusters(clusters.rej, FDR.loc = 0.01)
save(clusters.trimmedSN, file=paste0("BiSeqrunchr",i,"cimp.CRC.RData"))

#build regions from DMCs
DMRsSN <- findDMRs(clusters.trimmedSN, max.dist = 50, diff.dir = TRUE, alpha = 0.05)
save(DMRsSN, file=paste0("DMRschr",i,"cimp.CRC.RData"))
rm(list = ls())
}

### Do extra filters ####---------------------------------------------------------------------

#Load all DMR objects from all chromosomes (for each comparison) and make a single DMR table

load_filter_DMRS <- function(path, CpGpath, comparison){
  f <- list.files(path, paste0("DMRschr[0-9]+", comparison))
  numbers <- as.numeric(regmatches(f, regexpr("[0-9]+", f)))
  f <- paste0(path, f[order(numbers)])
  
  fc <- list.files(path, paste0("BiSeqrunchr[0-9]+", comparison))
  numbers <- as.numeric(regmatches(fc, regexpr("[0-9]+", fc)))
  fc <- paste0(path, fc[order(numbers)])
  
  
  load(file = f[1])
  load(file = fc[1])
  allDMRs <- DMRsSN
  allClusters.trimmed <- clusters.trimmedSN
  
  for(i in 2:22){
    load(file = f[i])
    load(file = fc[i])
    allDMRs <- c(allDMRs, DMRsSN)
    sorted <- clusters.trimmedSN[order(clusters.trimmedSN$pos),]
    allClusters.trimmed <- rbind(allClusters.trimmed, sorted)
  }
  
  biseqDMRs <- allDMRs[width(allDMRs) > 1]
  
  ###Check which DMRs overlap with a probe and how much
  
  probes <- read.table(file="probes/130912_HG19_CpGiant_4M_EPI.bed")
  probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3))
  
  hits <- findOverlaps(biseqDMRs, probes)
  gr.over <- pintersect(biseqDMRs[queryHits(hits)], probes[subjectHits(hits)])
  gr.counts <- tapply(gr.over, queryHits(hits), function(x) sum(width(x)))
  biseqDMRs$overlap <- 0
  biseqDMRs$overlap[as.numeric(names(gr.counts))] <- unname(gr.counts)
  biseqDMRs$percentageOverlap <- round(100*biseqDMRs$overlap / width(biseqDMRs))
  
  ###Get number of GpGs from region that are differential, and the total number of covered CpGs
  
  load(paste0(pathCpG, "GRangesAllCpGs.RData")) #contains all CpG sites covered in dataset
  
  #Total number of CpGs covered
  hits <- findOverlaps(biseqDMRs, GRallCpGs)
  biseqDMRs$totalCpGs <- 0
  biseqDMRs$totalCpGs[rle(queryHits(hits))$values] <- rle(queryHits(hits))$lengths
  
  # Number of diff CpGs
  diffCpGs <- GRanges(allClusters.trimmed$chr, IRanges(start = allClusters.trimmed$pos, width = 1))
  hits <- findOverlaps(biseqDMRs, diffCpGs)
  biseqDMRs$numDiffCpGs <- 0
  biseqDMRs$numDiffCpGs[rle(queryHits(hits))$values] <- rle(queryHits(hits))$lengths
  
  return(biseqDMRs)
}

#EXAMPLE for a comparison
biseqDMRs <- load_filter_DMRS(".", ".", "Aden" ) 

save(biseqDMRs, file = "AdenVsNorm.DMRs.RData")
save(allClusters.trimmed, file = "AdenVsNorm.DMCs.RData")

### Make master table

biseqTable <- (as(biseqDMRs, "data.frame"))

colnames(biseqTable) <- c("CHR", "START", "END", "WIDTH", "STRAND", "MEDIAN.P", "MEDIAN.METHYLATION.SSA", "MEDIAN.METHYLATION.NORM",
                          "MEDIA.METH.DIFFERENCE", "bpOVERLAP.WITH.PROBE", "PERCENTAGE.OVERLAP", "TOTAL.NUM.CpGs",
                          "NUM.DIFFERENTIAL.CpGs")

write.table(biseqTable, "serrated_table_DMRs_BiSeq_AdenVsNorm.csv", row.names=F, quote=FALSE, sep="\t")


