#!/usr/bin/env Rscript
# run as [R < scriptName.R --no-save]

#########################################################################################
# R script to annotate DMRs obtaibned from BiSeq, and make a scatter plot of logFC of expression
# vs DMR mean methylation, and add the number of points for each quadrant
# 
# Stephany Orjuela, February 2018
#########################################################################################

library(BiSeq)
library(ChIPpeakAnno)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(cowplot)
library(edgeR)
library(ensembldb)

####Make table to plot ####----------------------------------------------------------------------------
#to make full plot I use all DMRs without filtering
load(file = "SSAvsNorm.DMRs.RData")

SSADMRs <- biseqDMRs[biseqDMRs$percentageOverlap >= 25]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")
SSADMRsTable <- as(SSADMRs, "data.frame")

load(file = "AdenvsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[biseqDMRs$percentageOverlap >= 25]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo")
AdenDMRsTable <- as(AdenDMRs, "data.frame")

#Choose promoter window or grab genes to annotate to
annoData <- toGRanges(EnsDb.Hsapiens.v75) #option 1 
proms <- promoters(annoData, upstream = 2000, downstream = 2000) #option 2

#load expression
load(file="../RNAseq/Salmon/DGE_testsBlock_filtGenes_withVM_moreGenes_novEdit_cdnaFix_2fit.Rdata")

#function to filter out repeated DMRs
choose1perGene <- function(annotatedDMRs, factor){
  redanno <- annotatedDMRs
  x <- unique(annotatedDMRs$feature)
  
  redanno <- redanno[abs(redanno$shortestDistance) < 2000,]

  for(i in x){
    w <- which(redanno$feature == i)
    if(length(w) > 1){
      column <- grep(factor, colnames(mcols(redanno)))
      closest <- which.min(mcols(redanno)[,column][w])
      redanno <- redanno[-w[-closest],]
    }
  } 
  return(redanno)
}

#Annotation with proximity to TSS: one DMR per gene, one gene per DMR
plotCorr <- function(DMRs, annoData, de, lrt, comparison){
  annoM <- ChIPpeakAnno::annotatePeakInBatch(DMRs,
                                             FeatureLocForDistance="TSS",
                                             AnnotationData = annoData,
                                             output = "nearestLocation",
                                             PeakLocForDistance= "middle")
                                             
  #function from above
  redanno <- choose1perGene(annoM, "shortestDistance")
  redannoDF <- data.frame(ID = redanno$feature, 
                          median.p = redanno$median.p, 
                          median.diff = redanno$median.meth.diff)
  #x <- which(de != 0)
  df <- data.frame(ID = lrt$genes$ENSEMBLID,
                   Gene = lrt$genes$GENESYMBOL,
                   chr = lrt$genes$CHROM,
                   start = lrt$genes$START,
                   end = lrt$genes$END,
                   pval = lrt$table$PValue, 
                   logFC = lrt$table$logFC)
  #df <- df[x,]
  jointBDF <- merge(df, redannoDF, by="ID")
  return(jointBDF)
}

jointSN <- plotCorr(SSADMRs, annoData, de[,2], lrt[[2]], "SN") 
jointAN <- plotCorr(AdenDMRs, annoData, de[,1], lrt[[1]], "AN")

save(jointAN, jointSN, file ="correlationDataFull.RData")
#load("correlationDataFull.RData")

#Find interesting genes based on other authors
delkerGenes <- c("ENSG00000139800", "ENSG00000124233", "ENSG00000253368", "ENSG00000184956",
                 "ENSG00000163499", "ENSG00000075618", "ENSG00000043355", "ENSG00000076242",
                 "ENSG00000189280", "ENSG00000169903", "ENSG00000109511", "ENSG00000163347",
                 "ENSG00000175592", "ENSG00000144908", "ENSG00000129451", "ENSG00000036565",
                 "ENSG00000112299", "ENSG00000169876")

#LEF1-AS1 is a non-coding RNA
SSAmarkers <- c("ENSG00000113083", "ENSG00000106483", "ENSG00000232021")

#LINC is also a ncRNA
panmarkers <- c("ENSG00000170743", "ENSG00000188227")

####Setup to plot ####-----------------------------------------------------------

geneA <- jointSN[jointSN$Gene %in% c("ZIC2","HUNK"),]

#Get numbers for plot
getQ <- function(jointSN, meth.lim, exp.lim){
  quadrants <- list(
                jointSN[jointSN$median.diff >= meth.lim & jointSN$logFC >= exp.lim,], #pos.pos
                jointSN[jointSN$median.diff >= meth.lim & jointSN$logFC <= -exp.lim,], #pos.neg
                jointSN[jointSN$median.diff <= -meth.lim & jointSN$logFC <= -exp.lim,], #neg.neg
                jointSN[jointSN$median.diff <= -meth.lim & jointSN$logFC >= exp.lim,]) #neg.pos
  return(quadrants)
}

quadrants <- getQ(jointSN, 0.10,1)
dark.datos <- data.frame(  
       numb = sapply(quadrants, function(x){
                                   l <- length(x$ID)
                                   paste("n =", l)}))
quadrants <- getQ(jointSN,0,0)
all.datos <- data.frame(  
  numb = sapply(quadrants, function(x){
    l <- length(x$ID)
    paste("Total =", l)}))

p1 <- ggplot(jointSN, aes(x=median.diff, y=logFC, color = abs(logFC) >= 1 & abs(median.diff) >= 0.1)) +

  geom_point(size = 0.5) + 
  theme_light(base_size = 15) +
  guides(color=FALSE) + 
  scale_color_manual(values=c("#F5D6E4", "#CD3278")) + 
  geom_vline(xintercept=0, col="black", lty=1, lwd=0.5) +
  geom_hline(yintercept=0, col="black", lty=1, lwd=0.5) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-8,8)) +
  #labs(x="", y="logFC of expression", title = "SSA - Norm") +
  geom_label(data=geneA, aes(x=median.diff, y=logFC, label=Gene), size=5, hjust = 0, nudge_x = 0.02,colour = "black",label.size = 1) +
  annotate("text", x = c(0.3,0.3,-0.3,-0.3), y = c(7.5,-6.5,-6.5,7.5), label = dark.datos$numb, size = 4.5) + #colored
  annotate("text", x = c(0.3,0.3,-0.3,-0.3), y = c(7,-7,-7,7), label = all.datos$numb, size = 3.5) #total

geneA <- jointAN[jointAN$Gene %in% c("ZIC2","HUNK"),]

#Get numbers
quadrants <- getQ(jointAN, 0.10,1)
dark.datos <- data.frame(  
  numb = sapply(quadrants, function(x){
    l <- length(x$ID)
    paste("n =", l)}))

quadrants <- getQ(jointAN,0,0)

all.datos <- data.frame(  
  numb = sapply(quadrants, function(x){
    l <- length(x$ID)
    paste("Total =", l)}))


p2 <- ggplot(jointAN, aes(x=median.diff, y=logFC, color = abs(logFC) > 1 & abs(median.diff) > 0.1)) +
  geom_point(size = 0.5) + 
  theme_light(base_size = 15) +
  guides(color=FALSE) + 
  scale_color_manual(values=c( "#B9C9D7","#38678f")) + 
  #scale_color_manual(values="#38678f") + 
  geom_vline(xintercept=0, col="black", lty=1, lwd=0.5) +
  geom_hline(yintercept=0, col="black", lty=1, lwd=0.5) +
  coord_cartesian(xlim = c(-0.5, 0.5), ylim = c(-8,8)) +
  #labs(x="", y="", title = "Adenoma - Norm") +
  geom_label(data=geneA, aes(x=median.diff, y=logFC, label=Gene), size=5, hjust = 0, nudge_x = 0.02,colour = "black", label.size = 1) +
  annotate("text", x = c(0.3,0.3,-0.4,-0.4), y = c(7.5,-6.5,-6.5,7.5), label = dark.datos$numb, size = 4.5) +
  annotate("text", x = c(0.3,0.3,-0.4,-0.4), y = c(7,-7,-7,7), label = all.datos$numb, size = 3.5)

plot_grid(p1, p2)
ggsave("corrPlots.pdf", width = 15, height = 10)

#### Venn diagrams using genes ####-------------------------------

#also make master table with this filtering

load(file = "SSAvsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25 & biseqDMRs$median.p <= 0.01]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo") #26435

load(file = "AdenvsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25 & biseqDMRs$median.p <= 0.01]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo") #34382

library(VennDiagram)
annotateDMRs <- function(DMRs, annoData){
  annoM <- ChIPpeakAnno::annotatePeakInBatch(SSADMRs,
                                             FeatureLocForDistance="TSS",
                                             AnnotationData = annoData,
                                             output = "nearestLocation",
                                             PeakLocForDistance= "middle")
  
  #function from above
  redanno <- choose1perGene(annoM, "shortestDistance")
  return(redanno)
}

annoSN <- annotateDMRs(SSADMRs, annoData)
annoAN <- annotateDMRs(AdenDMRs, annoData)
save(annoSN, annoAN, file ="annotatedDMRsFiltered.RData")
#load("annotatedDMRsFiltered.RData")

pairVenn <- function(anno1, anno2, comp, colors, orientation){
  grid.newpage()
  
  As <- anno1$feature[anno1$state == orientation]
  Ss <- anno2$feature[anno2$state == orientation]
  hits <- intersect(As, Ss)
  draw.pairwise.venn(area1 = length(Ss),
                       area2 = length(As),
                       cross.area = length(hits),
                       #inverted = T,
                       category = comp,
                       lwd = c(5,5),
                       col = colors,
                       fill = colors,
                       alpha = rep(0.3, 2),
                       cex = rep(2,3),
                       label.col = rep("gray30", 3),
                       fontfamily = rep("sans", 3),
                       #cat.pos = c(-5,5),
                       cat.dist = rep(0.05,1 ),
                       cat.col = colors,
                       cat.fontfamily = rep("sans", 2))
  }

pdf("vennDiagramDMRsperGene.pdf")
pairVenn(annoAN, annoSN, c("SSA/Ps","Adenomas"), c("#CD3278","#38678f"), "hyper")
pairVenn(annoAN, annoSN, c("SSA/Ps","Adenomas"), c("#CD3278","#38678f"), "hypo")
dev.off()