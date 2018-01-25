#!/usr/bin/env Rscript
#########################################################################################
# R script to make diagnostic plots for DGE analysis
# #
# RNA-seq data set with 33 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesion with normal tissue, and 16 adenoma lesion with normal tissue
# 
# Stephany Orjuela, January 2018
#########################################################################################

setwd("/run/user/1000/gvfs/sftp:host=imlssherborne/home/sorjuela/RNAseq/Salmon")

#load libraries and data from DGEanalysis.R
library(edgeR)
library(rgl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(VennDiagram)
load("DGE_testsBlock_filtGenes_withVM.Rdata") #change
load("DGE_obj_Block_withVM.RData") #change

##Venn diagram ####---------------------------------------------------------------------------------------

pdf("VennDiagrams.pdf")

vennDiagram(des_b, include=c("up", "down" ), counts.col = c("gray30", "gray50"), 
            cex=1.2, circle.col = c("#38678f","#CD3278","#ffa500"), mar = rep(0,4), lwd = 5)

vennDiagram(des_b, counts.col = "gray30", 
            cex=1.2, circle.col = c("#38678f","#CD3278","#ffa500"), mar = rep(0,4), lwd = 5)

vennDiagram(des_b[,1:2], include=c("up", "down" ), counts.col = c("gray30", "gray50"), 
            cex=1.2, circle.col = c("#38678f","#CD3278"), mar = rep(0,4), lwd = 5)
dev.off()

#Another option
# pairVenn = function(n1, n2, counts, name, comp, colors){
#   #pdf(paste0("Venn_",name,".pdf"))
#   grid.newpage()
#   draw.pairwise.venn(area1 = n1, 
#                      area2 = n2, 
#                      cross.area = counts[4,3], 
#                      category = comp, 
#                      lty = rep("blank",2), 
#                      fill = colors, 
#                      alpha = rep(0.5, 2),
#                      cex = rep(1.5,3),
#                      label.col = rep("gray30", 3),
#                      fontfamily = rep("sans", 3),
#                      cat.pos = c(-5,5), 
#                      cat.dist = rep(0.025, 2),
#                      cat.col = rep("gray30", 2),
#                      cat.fontfamily = rep("sans", 2))
#   #dev.off()
# }
# 
# par(mfrow=c(2,2))
# #Total
# countsANSN = vennCounts(des_b[,1:2])
# countsANAS = vennCounts(des_b[,c(1,3)])
# countsSNAS = vennCounts(des_b[,2:3])
# AN = length(which(des_b$`ADENOMA-NORMAL` != 0))
# SN = length(which(des_b$`SSA-NORMAL` != 0))
# AS = length(which(des_b$`SSA-ADENOMA` != 0))
# pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))
# 
# #down
# countsANSN = vennCounts(des_b[,1:2]) #choose only down genes
# AN = length(which(des_b$`ADENOMA-NORMAL` != 0))
# SN = length(which(des_b$`SSA-NORMAL` != 0))
# pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))
# 
# #up
# pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))
# #pairVenn(AN, AS, countsANAS, "ANAS", c("Adenoma Vs Normal", "SSA/P Vs Adenoma"), c("#38678f","#ffa500"))
# #pairVenn(SN, AS, countsSNAS, "SNAS", c("SSA/P Vs Normal", "SSA/P Vs Adenoma"), c("#CD3278","#ffa500"))
# 
# counts = vennCounts(des_b)
# grid.newpage()
# draw.triple.venn(area1 =AN, 
#                  area2 = SN, 
#                  area3 = AS, 
#                  n12 = countsANSN[4,3], 
#                  n23 = countsSNAS[4,3], 
#                  n13 = countsANAS[4,3], 
#                  n123 = counts[8,4], 
#                  category = c("Adenoma Vs Normal", "SSA/P Vs Normal", "SSA/P Vs Adenoma"), 
#                  lty = "blank", 
#                  fill = c("#38678f","#CD3278","#ffa500"),
#                  alpha = 0.5,
#                  cex = 1.5,
#                  label.col = "gray30",
#                  fontfamily = "sans", 
#                  cat.dist = 0.05,
#                  cat.col = "gray30",
#                  cat.fontfamily = "sans",
#                  cat.cex = 0.8)
# dev.off()

###Smear plot ###----------------------------------------------------------------------------
smear = function(lrt, de, name){
  isde = which(de != 0)
  de.genes <- rownames(lrt$table[isde,])
  pdf(paste0("Smear_plot_",name,".pdf"))
  plotSmear(lrt, smooth.scatter = T)
  plotSmear(lrt, de.tags = de.genes, smooth.scatter = T)
  dev.off()
}

smear(lrt[[1]], de[,1], "AN")
smear(lrt[[2]], de[,2], "SN")
smear(lrt[[3]], de[,3], "SA")

###Volcano plots###-----------------------------------------------------------------------

#ggplot version
volcano <- function(lrt, name){
  level <- ifelse(lrt$table$logFC >= 1, "Up", "--")
  level[which(lrt$table$logFC <= -1)]  <- "Down"
  adjustedp <- p.adjust(lrt1$table$PValue,method="BH")
  level[which(adjustedp >= 0.05)] <- "--" 
  volcanoData = data.frame(logFC = lrt$table$logFC, negLogPval = -log10(adjustedp), level = level)
  ggplot(volcanoData, aes(x=logFC, y=negLogPval, colour = level)) + geom_point() + 
    theme_bw() +
    geom_vline(xintercept=c(-1, 1), col="black", lty=2, lwd=0.5) +
    geom_hline(yintercept=1.30, col="black", lty=2, lwd=0.5) +
    scale_color_manual(values=c("black", "red", "green")) + 
    coord_cartesian(xlim = c(-10, 10), ylim =c(0,130)) 
  ggsave(paste0("Volcano_plot_",name,".pdf"))
}

volcano(lrt[[1]], "AN")
volcano(lrt[[2]], "SN")
volcano(lrt[[3]], "SA")

#Smear version of volcano
volcano.smooth <- function(lrt, name){
  level <- ifelse(lrt1$table$logFC >= 1, "Up", "--")
  level[which(lrt1$table$logFC <= -1)]  = "Down"
  adjustedp <- p.adjust(lrt1$table$PValue,method="BH")
  level[which(adjustedp >= 0.05)] <- "--" 
  volcanoData = data.frame(logFC <- lrt1$table$logFC, negLogPval = -log10(adjustedp), level = level)
  
  pdf(paste0("SmoothVolcano_plot_",name,".pdf"))
  Lab.palette <- colorRampPalette(c("white", "orange", "red"), space = "Lab")
  smoothScatter(volcanoData$logFC, volcanoData$negLogPval, nbin = 200,
                colramp = Lab.palette)
  abline(v = c(-1,1), h = 1.30, lty = 2)
  dev.off()
}

volcano.smooth(lrt[[1]], "AN")
volcano.smooth(lrt[[2]], "SN")
volcano.smooth(lrt[[3]], "SA")

##Gene clustering ####-------------------------------------------------------------------------------------------------

#biased - DE genes
isde <- des_b != 0
table(isde)
heatmap_genes <- apply(as.data.frame(isde), 1, any)
table(heatmap_genes)

heatmap_cpm <- edgeR::cpm(y, log=TRUE, prior.count=10)
colnames(heatmap_cpm) <-paste0(y$samples$samples,"_",cond_forplots_fix)
heatmap_data <- heatmap_cpm[heatmap_genes, ]

annot_col <- data.frame(Tissue = cond_forplots_fix, 
                        Gender = gender,
                        Mutation = factor(c("KRAS","NA","NA","BRAF","NA", "BRAF","WT",
                                            "NA","NA","BRAF","WT","NA","NA","BRAF",
                                            "NA","BRAF","WT","NA","WT","NA","NA",
                                            "BRAF","NA","BRAF","NA","BRAF","NA","BRAF",
                                            "KRAS","NA","NA","BRAF","NA","BRAF","NA",
                                            "KRAS","NA","KRAS", "WT","WT","NA","KRAS",
                                            "NA","KRAS","NA","KRAS","NA", "WT","NA",
                                            "NA","BRAF","WT","NA","NA","BRAF","KRAS",
                                            "NA","NA","BRAF","NA","BRAF","NA","WT", "NA")),
                        Age = factor(agegroup))

rownames(annot_col) <- paste0(y$samples$samples,"_",cond_forplots_fix)
heatmap.colors <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")))(50)
ann_colors <- list(Tissue = c('NORMAL.ADENOMA'="#AADEDB",'ADENOMA'="#38678f",'NORMAL.SSA/P'="#EAABC8",'SSA/P'="#CD3278"), 
                  Gender = c('F'="gray56", 'M'="gray12"),
                  Mutation = c('NA'="white", 'WT'="#4DAF4A", 'KRAS'="#377EB8",'BRAF'="#E41A1C"),
                  Age = c('<50'="#CCEBC5", '50-70'="#FBB4AE", '>70'="#B3CDE3"))

pdf("biased.heatmap.DEGenes.pdf")
pheatmap(heatmap_data, fontsize = 5, color = heatmap.colors, annotation_col = annot_col,
         annotation_colors=ann_colors, show_rownames =F, border_color = NA, clustering_method = "complete", 
         scale="row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")
dev.off()

#unbiased - most variable genes
#The heatmap with the amount by which each gene deviates in a specific sample from the gene’s average across all samples.
#Hence, we center each genes’ values across samples

cps <- edgeR::cpm(y)
zz <- as.matrix(log2(1+cps))
set.seed(1)
rv <- rowVars(zz)
idx <- order(-rv)[1:50] #change number
mat  <- zz[idx,]
mat  <- mat - matrixStats::rowMeans2(mat)
rownames(mat) <- mapIds(org.Hs.eg.db, keys=row.names(mat), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
colnames(mat) <- paste0(y$samples$samples,"_",cond_forplots_fix)

pdf("unbiased.heatmap.50mostVarGenes.pdf")
pheatmap(mat, fontsize = 5, color = heatmap.colors, annotation_col = annot_col,
         annotation_colors=ann_colors, show_rownames =F, border_color = NA, clustering_method = "complete", 
         scale="row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")
dev.off()