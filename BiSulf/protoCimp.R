#!/usr/bin/env Rscript

#########################################################################################
# R script to plot CpGsites overlapping CpG islands and/or genes (proto-CIMP)
#
# BS-seq data set with 31 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 16 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
#
# Stephany Orjuela, February 2018
#########################################################################################

library(ggplot2)
library(gridExtra)
library(BiSeq)

###barplot ####---------------------------------------------------------------------------------

#CpG islands
islands = read.csv("probes/CpGIslandTrack.bed", header =F, sep="\t")
islandR = GRanges(seqnames = islands$V1, ranges = IRanges(islands$V2, islands$V3)) #28691

#probes
probes <- read.table(file="probes/130912_HG19_CpGiant_4M_EPI.bed")
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3)) #240131

#Get islands that overlap a probe
hits <- findOverlaps(islandR, probes)
red.islandR <- islandR[unique(queryHits(hits))] #26003
red.islands <- islands[unique(queryHits(hits)),]

#load SSA DMRs and select hypermethylated sites
load("SSA.allclusters.trimmed.RData")
SSAsites <- allClusters.trimmed
load("Adenoma.allclusters.trimmed.RData")
Adensites <- allClusters.trimmed

filterSites <- function(sites, red.islandR, tissue){
  allClust <- sites[sites$meth.diff > 0,]
  allclust.GR <- GRanges(allClust$chr, IRanges(start = allClust$pos, width = 1))
  
  #Overlap with islands and count
  hits <- findOverlaps(red.islandR, allclust.GR)
  red.islandR$numDMCs <- 0
  red.islandR$numDMCs[rle(queryHits(hits))$values] <- rle(queryHits(hits))$lengths
  
  DMCs <- red.islandR$numDMCs[order(red.islandR$numDMCs)]
  toPlot <- data.frame(counts = rle(DMCs)$lengths, values = rle(DMCs)$values, tissue = tissue)
  return(toPlot)
}

SSAtab <- filterSites(SSAsites, red.islandR, "SSA/Ps")
Adentab <- filterSites(Adensites, red.islandR, "Adenomas")
tab <- rbind(SSAtab, Adentab)

pdf("protoCimpPlot.pdf")
ggplot(tab, aes(x=values, y=log(counts), fill = tissue)) + 
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_light(base_size = 15) +
  scale_fill_manual(values=c('#CD3278', '#38678f')) +
  coord_cartesian(xlim = c(0, 15)) +
  labs(x="Number of hypermethylated CpGs per canonical CpG Island", 
       y="log - Number of CpG Islands")

dev.off()

##### Another alternative to showing protocimp from the extension perspective ####-----------------

#Grab TSSs from ensembl
library(biomaRt)
grch37 <- useEnsembl(biomart="ensembl",GRCh=37, dataset = "hsapiens_gene_ensembl")
tss <- getBM(attributes = c("start_position", "transcription_start_site", "chromosome_name",
                            "transcript_start", "transcript_end",
                            "strand",  "ensembl_gene_id",
                            "ensembl_transcript_id", "external_gene_name"),
             mart = grch37)

GRtss <- GRanges(tss$chromosome_name, 
                 IRanges(start = tss$transcription_start_site, width = 1),
                 strand = ifelse(tss$strand == 1, "+", "-"))

#slop them
GRproms <- GenomicRanges::promoters(GRtss, upstream=2000, downstream=2000)

#Turn clus.trimm df to GR
DFtoGR <- function(sites){
  allClust <- sites[sites$meth.diff >= 0.10 & sites$p.li <= 0.01,] #for hypermeth
  #chroms <- gsub("chr", "", allClust$chr) #Change if seqnames without "chr"
  chroms <- allClust$chr
  allclust.GR <- GRanges(chroms, IRanges(start = allClust$pos, width = 1))
  return(allclust.GR)
}
  
GRsitesSSA <- DFtoGR(SSAsites) #227145
GRsitesAden <- DFtoGR(Adensites) #156908

#Fill in columns on main window table for each treatment

tss$sitesSSA <- 0
tss$sitesAdenoma <- 0

makeDFtoPlot <- function(GRsites1, GRsites2, window, table){
  hits <- findOverlaps(window, GRsites1)
  q <- sort(queryHits(hits))
  table$sitesSSA[rle(q)$values] <- rle(q)$lengths
  hits <- findOverlaps(window, GRsites2)
  q <- sort(queryHits(hits))
  table$sitesAdenoma[rle(q)$values] <- rle(q)$lengths
  return(table)
}

tss <- makeDFtoPlot(GRsitesSSA, GRsitesAden, GRproms, tss) #215170

#Filter table
x <- grep("[A-Za-z]", tss$chromosome_name)
tss_filt <- tss[-x,]
tss_filt <- tss_filt[rowSums(tss_filt[,9:10]) > 0,] #17154

save(tss_filt, file = "NumCpGsitesperProm.2kbwindow.RData")
#load("NumCpGsitesperProm.2kbwindow.RData")

pdf("protoCimp_scatter.pdf")
ggplot(tss_filt, aes(sitesAdenoma, sitesSSA)) + 
  geom_point() + 
  #geom_bin2d(bins = 60) +
  theme_light(base_size = 15) + 
  labs(x="Numper of hyper-methylated CpG sites Adenoma", y="Numper of hyper-methylated CpG sites SSA", 
       title = "Window: -/+ 2kb from TSS") +
  geom_smooth(method = "loess")

### For CpG islands

red.islands$sitesSSA <- 0
red.islands$sitesAdenoma <- 0

red.islands <- makeDFtoPlot(GRsitesSSA, GRsitesAden, red.islandR, red.islands)

#Filter table
red.islands <- red.islands[rowSums(red.islands[,5:6]) > 0,] #7578
save(red.islands, file = "NumCpGsitesperIsland.RData")
load("NumCpGsitesperIsland.RData")

ggplot(red.islands, aes(sitesAdenoma, sitesSSA)) + 
  geom_point() + 
  #geom_bin2d(bins = 60) + 
  theme_light(base_size = 15) + 
  labs(x="Numper of hyper-methylated CpG sites Adenoma", y="Numper of hyper-methylated CpG sites SSA", 
       title = "Window: Full CpG Island") +
  geom_smooth(method = "lm")

dev.off()
