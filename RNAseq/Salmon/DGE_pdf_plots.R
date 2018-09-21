#!/usr/bin/env Rscript
#########################################################################################
# R script to make diagnostic plots for DGE analysis
# #
# RNA-seq data set with 33 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesion with normal tissue, and 16 adenoma lesion with normal tissue
# 
# Stephany Orjuela, January 2018
#########################################################################################


#load libraries and data generated from DGEanalysis.R
library(edgeR)
library(rgl)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(ggplot2)
library(VennDiagram)

load(file = "RNAseq/Data/DGE_obj.RData")
load(file= "RNAseq/Data/DGE_tests.Rdata")

##Venn diagram ####---------------------------------------------------------------------------------------

#To check
vennDiagram(des_b[,1:2], include=c("up", "down" ), counts.col = c("gray30", "gray50"), 
           cex=1.2, circle.col = c("#38678f","#CD3278"), mar = rep(0,4), lwd = 5)


#For actual plot in paper

pairVenn <- function(n1, n2, counts, comp, colors){
  
  grid.newpage()
  draw.pairwise.venn(area1 = n1,
                     area2 = n2,
                     cross.area = counts,
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

pdf("DGEvennDiagrams.pdf")

#up
countsANSN <- sum(rowSums(des_b[,1:2] > 0 ) == 2)
AN = sum(des_b$`ADENOMA-NORMAL` == 1)
SN = sum(des_b$`SSA-NORMAL`== 1)
pairVenn(SN, AN, countsANSN, c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))

#down
countsANSN <- sum(rowSums(des_b[,1:2] < 0 ) == 2)
AN = sum(des_b$`ADENOMA-NORMAL` == -1)
SN = sum(des_b$`SSA-NORMAL`== -1)
pairVenn(SN, AN, countsANSN, c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))

dev.off()

###Volcano plots ###-------------------------------------------------------------------------

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
heatmap_genes <- apply(as.data.frame(isde), 1, any)

heatmap_cpm <- cpm(y, log=TRUE, prior.count=10)
colnames(heatmap_cpm) <-paste0(y$samples$samples,"_",cond_forplots_fix)
heatmap_data <- heatmap_cpm[heatmap_genes, ]

#Make table for annotation (not sure why the mutation was not in the metadata tables)
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

#Make asimetrical color scheme
heatmap.colors <- colorRampPalette(rev(brewer.pal(11, "RdYlGn")))(15)
w <- heatmap.colors[-c(6:8)]
heatmap.colors <- c(rep(heatmap.colors[1], 8), w, rep(heatmap.colors[15], 8))

#Make table for annotation colors
ann_colors <- list(Tissue = c('NORMAL.ADENOMA'="#AADEDB",'ADENOMA'="#38678f",'NORMAL.SSA/P'="#EAABC8",'SSA/P'="#CD3278"), 
                  Gender = c('F'="gray56", 'M'="gray12"),
                  Mutation = c('NA'="white", 'WT'="#4DAF4A", 'KRAS'="#377EB8",'BRAF'="#E41A1C"),
                  Age = c('<50'="#CCEBC5", '50-70'="#FBB4AE", '>70'="#B3CDE3"))

#plot
pdf("biased.heatmap.DEGenes.pdf")
pheatmap(heatmap_data, fontsize = 5, color = heatmap.colors, annotation_col = annot_col, kmeans_k = 20,
         annotation_colors=ann_colors, show_rownames =F, border_color = NA, clustering_method = "complete", 
         scale="row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")
dev.off()

#unbiased - most variable genes
#The heatmap with the amount by which each gene deviates in a specific sample from the gene’s average across all samples.
#So, I center each genes values across samples

zz <- as.matrix(log2(1+cpm(y)))
zz <- as.matrix(log2(y$counts+1))
rv <- rowVars(zz)

#top variable genes
idx <- order(-rv)[1:100] #change number
mat  <- zz[idx,]
mat  <- mat - matrixStats::rowMeans2(mat)
colnames(mat) <- paste0(y$samples$samples,"_",cond_forplots_fix)

pdf("unbiased.heatmap_top10000varGenes.pdf")
pheatmap(mat, fontsize = 10, color = heatmap.colors, annotation_col = annot_col,
         annotation_colors=ann_colors, show_rownames =F, border_color = NA, clustering_method = "complete", 
         scale="row", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",
         fontsize_col = 6)
dev.off()
