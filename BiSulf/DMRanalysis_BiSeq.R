#!/usr/bin/env Rscript
# run as [R < scriptName.R --no-save]
# command Rscript doesn't work on sherborne

#########################################################################################
# R script to run a DMR analysis workflow using BiSeq starting from .cov files
#
# BS-seq data set with 31 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 16 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
# 
# most of code taken from BiSeq vignette
#
# Stephany Orjuela, January 2018
#########################################################################################

setwd("/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/sorjuela/serrated_pathway_paper/BiSulf")

library(BiSeq)

### Upload design matrix###-----------------------------------------------------------------------------
##patient ID, file name, lesion, sex, and other stuff

#files <- read.table("design_matrix.csv", stringsAsFactors = F)
files <- read.table("CRCmatrix.csv", stringsAsFactors = F)

### Read cov files into R ###--------------------------------------------------------------------------
#Here are all the CpG sites from all samples

#This object is reused in other script to run other DMR analysis
BSr <- readBismark(files$V2, colData= DataFrame(group = files$V3, row.names = paste0(files$V1,".",files$V3))) 
#BSr <- readBismark(files$V3, colData= DataFrame(group = files$V2, row.names = paste0(files$V1,".",files$V2))) 
save(BSr, file="BSrawCRC.RData") 

for(i in 1:22){

load(file="BSrawCRC.RData")
w <- which(seqnames(BSr) == paste0("chr",i))
BSrchr <- BSr[w, ]
rm(BSr)
save(BSrchr, file = paste0("chr",i,"BSraw.CRC.RData"))

### Define clusters ###------------------------------------------------------------ 

BSr.clust.unlim <- clusterSites(BSrchr,
                                groups = factor(colData(BSrchr)$group),
                                perc.samples = 0.8,
                                min.sites = 3,
                                max.dist =  500,
                                mc.cores = 30) #Here are only sites within clusters
#save(BSr.clust.unlim, file = "ClusterSiteschr1.RData")

### smooth the methylation values of CpG sites within the clusters ###----------- 
predictedMeth <- predictMeth(object = BSr.clust.unlim, h=1000, mc.cores = 30)
save(predictedMeth, file=paste0("predictedMethchr",i,".CRC.RData"))

load(paste0("predictedMethchr",i,".CRC.RData"))
#predictedMeth <- predictedMeth[,-47] #just for SSA

#subset for groups
#s <- which(colData(predictedMeth)$group == "SSA")
#n <- s-1

#a <- which(colData(predictedMeth)$group == "ADENOMA")
#n <- a + 1

a <- which(colData(predictedMeth)$group == "cimp")
n <- a + 1

SNpMeth <- predictedMeth[,c(a,n)]
#colData(SNpMeth)$group <- factor(colData(SNpMeth)$group, levels = c("ADENOMA", "NORMAL"))
colData(SNpMeth)$group <- factor(colData(SNpMeth)$group, levels = c("cimp", "NORMAL.cimp"))

betaResultsSN <- betaRegression(formula = ~group, link = "probit", object = SNpMeth, type = "BR", mc.cores = 30) 


### Test CpG clusters for differential methylation ###---------------------------------
#For details go to other non-paired script

predictedMethNullSN <- SNpMeth[,c(1:3,4:6)] 
colData(predictedMethNullSN)$group.null <- rep(c(1,2), 3) # number of patiens in this comparison
betaResultsNullSN <- betaRegression(formula = ~group.null, link = "probit", object = predictedMethNullSN, type="BR", mc.cores = 30)

varioSN <- makeVariogram(betaResultsNullSN)
vario.sm <- smoothVariogram(varioSN, sill = 0.9)
vario.aux <- makeVariogram(betaResultsSN, make.variogram=FALSE)
vario.sm$pValsList <- vario.aux$pValsList
locCor <- estLocCor(vario.sm)
clusters.rej <- testClusters(locCor, FDR.cluster = 0.05) 
clusters.trimmedSN <- trimClusters(clusters.rej, FDR.loc = 0.01)
save(clusters.trimmedSN, file=paste0("BiSeqrunchr",i,"cimp.CRC.RData"))

#build regions
DMRsSN <- findDMRs(clusters.trimmedSN, max.dist = 50, diff.dir = TRUE, alpha = 0.05)
save(DMRsSN, file=paste0("DMRschr",i,"cimp.CRC.RData"))
rm(list = ls())
}

### Do extra filters ####-------------------------------------------------------------------

#Load all DMR objects from all chroms (for each comparison) and make a single DMR table

load_filter_DMRS <- function(path, CpGpath, comparison){
  #f <- list.files(".", 'DMRschr[0-9]+Aden')
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
  
  #set.seed(12345)
  hits <- findOverlaps(biseqDMRs, probes)
  gr.over <- pintersect(biseqDMRs[queryHits(hits)], probes[subjectHits(hits)])
  gr.counts <- tapply(gr.over, queryHits(hits), function(x) sum(width(x)))
  biseqDMRs$overlap <- 0
  biseqDMRs$overlap[as.numeric(names(gr.counts))] <- unname(gr.counts)
  biseqDMRs$percentageOverlap <- round(100*biseqDMRs$overlap / width(biseqDMRs))
  
  ###Get number of GpGs from region that are differential, and the total number of covered CpGs
  
  load(paste0(pathCpG, "GRangesAllCpGs.RData")) 
  
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

#Example for a comparison
biseqDMRs <- load_filter_DMRS("CRCdata/cimp/", "CRCdata/", "cimp" ) 

#save(biseqDMRs, file = "CIMPvsNorm.DMRs.RData")
#save(allClusters.trimmed, file = "nonCIMP.allclusters.trimmed.RData")

#Make master table

biseqTable <- (as(biseqDMRs, "data.frame"))

colnames(biseqTable) <- c("CHR", "START", "END", "WIDTH", "STRAND", "MEDIAN.P", "MEDIAN.METHYLATION.SSA", "MEDIAN.METHYLATION.NORM",
                          "MEDIA.METH.DIFFERENCE", "bpOVERLAP.WITH.PROBE", "PERCENTAGE.OVERLAP", "TOTAL.NUM.CpGs",
                          "NUM.DIFFERENTIAL.CpGs")

write.table(biseqTable, "serrated_table_DMRs_BiSeqAll_AdenVsNormFix.csv", row.names=F, quote=FALSE, sep="\t")


