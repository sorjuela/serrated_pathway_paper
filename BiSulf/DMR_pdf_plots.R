#########################################################################################
# R script to make plots for DMRs
#
# BS-seq data set with 31 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
# And 6 paired samples with diagnosed colorectal cancer, divided in CIMP (3) and nonCIMP (3) 
# depending on the methylation state of MLH1
# 
# Stephany Orjuela, January 2018
#########################################################################################

library(ggplot2)
library(UpSetR)
library(BiSeq)
library(RColorBrewer)
library(ChIPpeakAnno)
library(genomation)
library(methylKit)
library(GenomicAlignments)
library(EnsDb.Hsapiens.v75)
library(Gviz)
library(cowplot)

#### DMRs plots ###-----------------------------------------------------------
load(file = "SSAVsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")
SSADMRsTable <- as(SSADMRs, "data.frame")

load(file = "AdenVsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo")
AdenDMRsTable <- as(AdenDMRs, "data.frame")


ann_text <- data.frame(lab = c(paste0("n = ",sum(SSADMRsTable$state == "hyper")), paste0("n = ",sum(SSADMRsTable$state == "hypo"))),
                       state = c("hyper", "hypo"))

pdf("Lengths.Hist.boxes.pdf")
hists <- ggplot(SSADMRsTable, aes(x=width)) + 
  scale_y_log10() +
  facet_grid(. ~ state) +
  geom_histogram(bins = 200, fill="#CD3278") +
  labs(x="", y="No. of DMRs") + 
  theme_light(base_size = 10) + 
  coord_cartesian(xlim = c(0, 1500)) +
  geom_text(data = ann_text, mapping = aes(x = 500, y = 2000, label = lab)) 

boxplots <- ggplot(SSADMRsTable, aes(1, width)) + 
  facet_grid(. ~ state) +
  coord_flip() +
  geom_boxplot(fill = "white", colour = "#CD3278") +
  theme_light(base_size = 10) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) +
  labs(y="DMR length", x="") + 
  scale_y_continuous(limit = c(0, 1500))

plot_grid(hists, boxplots, align = "hv", ncol = 1,
          rel_heights = c(1,0.3))
             #widths=c(4, 1), heights=c(1, 4))


ann_text <- data.frame(lab = c(paste0("n = ",sum(AdenDMRsTable$state == "hyper")), paste0("n = ",sum(AdenDMRsTable$state == "hypo"))),
                       state = c("hyper", "hypo"))

hists <- ggplot(AdenDMRsTable, aes(x=width)) + 
  scale_y_log10() +
  facet_grid(. ~ state) +
  geom_histogram(bins = 200, fill="#38678f") + 
  labs(x="", y="No. of DMRs") + 
  theme_light(base_size = 10) + 
  coord_cartesian(xlim = c(0, 1500)) +
  geom_text(data = ann_text, mapping = aes(x = 500, y = 2000, label = lab))

boxplots <- ggplot(AdenDMRsTable, aes(1, width)) + 
  facet_grid(. ~ state) +
  coord_flip() +
  geom_boxplot(fill = "white", colour = "#38678f") +
  theme_light(base_size = 10) +
  theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) +
  labs(y="DMR length", x="") + 
  scale_y_continuous(limit = c(0, 1500))

plot_grid(hists, boxplots, align = "hv", ncol = 1,
          rel_heights = c(1,0.3))

dev.off()

#### Volcano plots DMRs ####--------------------------------------------------------

#need all points even the ones without the filter
load(file = "SSAVsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[biseqDMRs$percentageOverlap >= 25]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")
SSADMRsTable <- as(SSADMRs, "data.frame")

load(file = "AdenVsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[biseqDMRs$percentageOverlap >= 25]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo")
AdenDMRsTable <- as(AdenDMRs, "data.frame")

volcano <- function(DMSites, name){
  level <- ifelse(DMSites$median.meth.diff >= 0.10, "Hyper", "--")
  level[DMSites$median.meth.diff <= -0.10]  = "Hypo"
  level[DMSites$median.p > 0.01] = "--" 
  
  #Add numbers to plot
  numbers <- c()
  numbers[1] <- paste0("n = ",dim(DMSites[level == "Hypo",])[1])
  numbers[2] <- paste0("n = ",dim(DMSites[DMSites$median.meth.diff > -0.10 & DMSites$median.meth.diff < 0.10 & DMSites$median.p <= 0.01,])[1])
  numbers[3] <- paste0("n = ",dim(DMSites[level == "Hyper",])[1])
  
  ggplot(DMSites, aes(x=median.meth.diff, y=-log10(median.p), colour = level)) + 
    geom_point() + 
    theme_bw(base_size = 15) +
    labs(x="Methylation difference per DMR", y = "-log10(P-value)") +
    geom_vline(xintercept=c(-0.1, 0.1), col="black", lty=2, lwd=0.5) +
    geom_hline(yintercept=2, col="black", lty=2, lwd=0.5) +
    scale_color_manual(values=c("black", "yellow", "blue")) +
    coord_cartesian(xlim = c(-0.6, 0.6), ylim =c(0,85)) +
    annotate("text", x = c(-0.3,0,0.3), y = c(75,75,75), label = numbers, size = 4)
  ggsave(paste0("Volcano_plot_",name,"DMRs.jpeg"), device = "jpeg", dpi = 300) #Too heavy for a .pdf
}

volcano(SSADMRsTable, "SN")
volcano(AdenDMRsTable, "AN")

#### Upset plots DMRs ####----------------------------------------------------------

#These have the full filters

load(file = "SSAVsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25 & biseqDMRs$median.p <= 0.01]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")

load(file = "AdenVsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25 & biseqDMRs$median.p <= 0.01]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo")


makeUpsetTable <- function(orientationA, orientationS, SSADMRs, AdenDMRs){
  x <- AdenDMRs[mcols(AdenDMRs)$state == orientationA]
  y <- SSADMRs[mcols(SSADMRs)$state == orientationS]

  hits <- findOverlaps(x, y)
  comm <- x[queryHits(hits)]
  
  uncomS <- y[-subjectHits(hits)]
  uncomA <- x[-queryHits(hits)]
  
  hyperm <- matrix(0, nrow = sum(length(comm), length(uncomA), length(uncomS)), ncol = 2)
  hyperm[0:length(comm),] <- 1
  hyperm[(length(comm)+1):(length(comm)+length(uncomA)),1] <- 1
  hyperm[(length(comm)+length(uncomA)+1):(length(comm)+length(uncomA)+length(uncomS)),2] <- 1
  
  
  colnames(hyperm)=c("Adenomas", "SSA/Ps")#, "NA-val", "NS-val")#,"ADENOMA-SSA")
  return(hyperm)
  #l <- list(uncomA, uncomS)
  #return(uncomS)
}

hypoM <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
hyperM <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)


#Make upset

pdf("upset.plots.DMRs.pdf")
upset(as.data.frame(hypoM), 
      sets = c("Adenomas","SSA/Ps"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMRs (hypomethylated)",
      main.bar.color = "#3b56d8",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 16000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMRs")

upset(as.data.frame(hyperM), 
      sets = c("Adenomas","SSA/Ps"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMRs (hypermethylated)",
      main.bar.color = "#e1cf22",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 16000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMRs")
dev.off()

pdf("upset.plot.DMRs.hyperSSA.hypoAden.pdf")
upset(as.data.frame(perSpoA), 
      sets = c("Adenomas","SSA/Ps"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMRs",
      main.bar.color = "#198C19",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      #mainbar.y.max = 16000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMRs")
dev.off()

#### Annotation cakes ####------------------------------------------------------

#Choose promoter window or grab genes to annotate to
annoData <- toGRanges(EnsDb.Hsapiens.v75)


quickanno <- function(){
annoM <- ChIPpeakAnno::annotatePeakInBatch(SSADMRs,
                                           FeatureLocForDistance="start",
                                           AnnotationData = annoData,
                                           output = "nearestLocation",
                                           PeakLocForDistance= "middle")




#Function to choose closest DMR per gene
choose1perGene <- function(annotatedDMRs, factor){
  redanno <- annotatedDMRs
  x <- unique(annotatedDMRs$feature)
  
  for(i in x){
    w = which(redanno$feature == i)
    if(length(w) > 1){
      column = grep(factor, colnames(mcols(redanno)))
      closest = which.min(mcols(redanno)[,column][w])
      redanno = redanno[-w[-closest],]
    }
  } 
  return(redanno)
}

annoSN <- choose1perGene(annoM, "shortestDistance")
annoAN <- choose1perGene(annoM, "shortestDistance")
}

#pdf("pieCharts.DMRfeatureOverlap.pdf", 10,7.5)
piecolor <- brewer.pal(11, "Spectral")[c(3:5,9:11)]
par(mfrow=c(1,2))

pie1(table(as.data.frame(annoSN)$insideFeature), 
     col = piecolor, 
     border = NA, 
     main = "SSA/Ps", 
     legend = T, 
     legendpos = "bottom", 
     radius = 1)

pie1(table(as.data.frame(annoAN)$insideFeature),
     col = piecolor, 
     border = NA, 
     main = "Adenomas", 
     legend = T, 
     legendpos = "bottom",
     radius = 1)
dev.off()

####Use refseq Genes to make cakes ####--------------------------------------------------
gene.obj <- readTranscriptFeatures("genesRefSeq.bed", up.flank=1000,down.flank=500)

piecolor <- brewer.pal(11, "Spectral")[c(3,5,9:11)]

pieTables <- list(SSADMRs[SSADMRs$state == "hyper"],
              SSADMRs[SSADMRs$state == "hypo"],
              AdenDMRs[AdenDMRs$state == "hyper"],
              AdenDMRs[AdenDMRs$state == "hypo"])

pieNames <- c("SSA.hyper","SSA.hypo","Adenoma.hyper","Adenoma.hypo")

pdf("Cakes_DMRs_GeneFeatures.pdf")
for(i in 1:4){
  diffGeneann <- annotateWithGeneParts(pieTables[[i]],gene.obj)
  genomation::plotTargetAnnotation(diffGeneann, 
                                   col = piecolor, 
                                   border = NA, 
                                   main = pieNames[i], 
                                   radius = 1)
  }
dev.off()

#cgis
cpg.obj <- readFeatureFlank("probes/CpGIslandTrack.bed", feature.flank.name=c("CpGi","shores"))

pdf("Cakes_DMRs_CpGiFeatures.pdf")
for(i in 1:4){
  diffCpGann <- annotateWithFeatureFlank(pieTables[[i]], 
                                         cpg.obj$CpGi,
                                         cpg.obj$shores, 
                                         feature.name="CpGi",
                                         flank.name="shores")
  genomation::plotTargetAnnotation(diffCpGann,
                                   col = piecolor, 
                                   border = NA, 
                                   main = pieNames[i], 
                                   radius = 1)
  }
dev.off()

####Extras ####--------------------------------
#Fix to get these numbers:
#uniqHypoAs <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
#uniqHyperSs <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)
#findOverlaps(uniqHyperSs, uniqHypoAs) #276

### Get number of DMRs per 4kb window and count for each comparison. Exclude ones overlapping
#annoData <- toGRanges(EnsDb.Hsapiens.v75)
#proms <- promoters(annoData, upstream = 2000, downstream = 2000) #64102

#Get all unique ranges for hypo
#hypoM <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
#hyperM <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)

# allM <- c(hypoM, hyperM)
# 
# mat <- matrix(0,length(proms), 4)
# 
# for(i in 1:4){
#   hits <- findOverlaps(proms, allM[[i]])
#   mat[rle(queryHits(hits))$values,i] <- rle(queryHits(hits))$lengths
# }
# 
# mcols(proms) <- mat
# colnames(mcols(proms)) <- c("unHypoA", "unHypoS", "unHyperA", "unHyperS")
# proms <- proms[rowSums(mat) > 0]
# gene_names <- annoData$gene_name[rowSums(mat) > 0]
# proms$gene_name <- gene_names
# 
# proms[mcols(proms)[,1] > 0 & mcols(proms)[,2] > 0] # 6 close hypos
# proms[mcols(proms)[,3] > 0 & mcols(proms)[,4] > 0] # 29 close hyper
# proms[mcols(proms)[,1] > 0 & mcols(proms)[,4] > 0] #122 close hypoA and hyperS
# 
# promdf <- as(proms, "data.frame")
# write.table(promdf, "closeDMRCount.4kbwindowFromTSS.csv", row.names=F, quote=FALSE, sep="\t")

