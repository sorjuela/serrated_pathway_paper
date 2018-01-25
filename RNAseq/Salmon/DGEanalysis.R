#########################################################################################
# R script to run a normal DGE analysis workflow starting from Salmon counts import
# ending with export of gene tables
#
# RNA-seq data set with 33 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesion with normal tissue, and 16 adenoma lesion with normal tissue
# 
# part of code taken from Charlotte Soneson and Mark Robinson
#
# Stephany Orjuela, January 2018
#########################################################################################


library(tximport)
library(ensembldb)
library(readr)
library(statmod)
library(edgeR)
library(EnsDb.Hsapiens.v75)
library(rgl)
library(ggplot2)
library(VennDiagram)
library(Biostrings)
library(pheatmap)
library(RColorBrewer)
library(genefilter)
library(matrixStats)

## load metadata ####-----------------------------------------------------------------------------------------------------
##patient ID, file name, lesion, sex, and other stuff

files = as.matrix(read.table("GeneTables/new_exp_setup.csv")) 
files2 = as.matrix(read.table("GeneTables/exp_setup_newFa.csv")) #without normal_SSA
bf = as.vector(files2[,1])

##Make tx2gene table ####-------------------------------------------------------------------------------------------------

## Read the fasta file in R using a function from the Biostrings package
cdna.ncrna <- readDNAStringSet("Homo_sapiens.GRCh37.cdna.ncrna.fa")

## Go through the sequence names and extract the required information
tx2gene <- data.frame(t(sapply(names(cdna.ncrna), function(nm) {
  tmp <- strsplit(nm, " ")[[1]]
  tx <- tmp[1]
  gene <- gsub("gene:", "", tmp[grep("^gene:", tmp)])
  c(tx = tx, gene = gene)
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL

##Process to input edgeR, normalize ####---------------------------------------------------------------------------------

txi_salmonsimplesum <- tximport(files = bf, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)

save(txi_salmonsimplesum, file="txi_ncrnaFIX.RData")
#load("txi_salmonsimplesum.RData")

#Calculate offset
cts <- txi_salmonsimplesum$counts
normMat <- txi_salmonsimplesum$length 
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))
files2 = as.data.frame(files2)
y <- DGEList(cts, samples = files2[,2], group = files2[,3])
y$offset <- t(t(log(normMat)) + o) 

y <- calcNormFactors(y)

#Filter lowly expressed genes
y <- y[which(rowSums(cpm(y) > 0.5) > 2), , keep.lib.sizes = FALSE] 
y <- calcNormFactors(y)

##Pre-Plots ####---------------------------------------------------------------------------------------------------------

#Set up factors for analysis and plots
cond <- factor(files2[,3])
cond_forplots <- factor(files[,3])
cond_forplots_fix <- gsub("NORMAL_SSA", "NORMAL.SSA", cond_forplots)
cond_forplots_fix <- gsub("SSA$", "SSA/P", cond_forplots_fix)
cond_forplots_fix <- factor(gsub("NORMAL$", "NORMAL.ADENOMA", cond_forplots_fix), 
                           levels = c("NORMAL.ADENOMA", "ADENOMA", "NORMAL.SSA/P", "SSA/P"))
gender <- factor(files2[,4])
age <- as.character(files2[,5])
agegroup <- ifelse(age <= 50, "<50", "50-70")
agegroup[which(age >= 70)] <- ">70"

# Set up condition colours
tissue.colors = c("#AADEDB","#38678f","#EAABC8","#CD3278")
#Light blue: #AADEDB
#Blue: #38678f
#light pink: #EAABC8
#Violet red: #CD3278

gender.colors <- c("gray56", "gray12")

##MDS plot ####------------------------------------------------------------------------------------------

#2D MDS plot
m1 = plotMDS(y, labels = y$samples$samples)$cmdscale.out
df <- data.frame(dim1 = m1[,1], dim2 = m1[,2])
label = y$samples$samples
ggplot(df)+
  scale_shape_identity() +
  geom_point(data=df, aes(x=dim1, y=dim2, shape = gender, colour = cond_forplots_fix), size=5) +
  geom_text(data=, aes(x=dim1, y=dim2-0.1, label = label), size=4)+
  scale_colour_manual(values=tissue.colors)+
  scale_shape_manual(values=c(19,17))+
  labs(x = "1", y = "2", title = "")+
  labs(shape = "Gender", colour = "Tissue")+
  theme_bw()

#3D MDS plot
m3 <- plotMDS(y, dim.plot=c(1,3), labels = y$samples$samples)
m3cmds = m3$cmdscale.out
par3d(cex=1)
plot3d(x = m3cmds[,1], y = m3cmds[,2], z = m3cmds[,3], xlab = "1", ylab = "2", zlab = "3", 
       col = tissue.colors[cond_forplots_fix], type="p", size=8)
#text3d( x = m3cmds[,1], y = m3cmds[,2], z = m3cmds[,3],cex = 0.8, font = 1, text = y$samples$samples, 
#        adj = 2, color=gender.colors[gender], family="sans") 

subid <- currentSubscene3d()
rglwidget(elementId="simplesum_RNAseq")

##Set up design matrix for GLM ####--------------------------------------------------------------------------------------
#this is an additive linear model with patient as the blocking factor

patient <- factor(y$samples$samples)
cond <- factor(files2[,3])
cond <- relevel(cond, "NORMAL")
design_byP <- model.matrix(~cond+patient) ##block by patient
#design <- model.matrix(~0+cond) ##or not
colnames(design_byP) <- c(levels(cond), levels(y$samples$samples)[-1])
#colnames(design) <- levels(cond)

#GLM estimates of dispersion
y <- estimateGLMCommonDisp(y, design_byP)
y <- estimateGLMTrendedDisp(y, design_byP) 
y <- estimateGLMTagwiseDisp(y, design_byP)

##Grab gene identifiers ####-------------------------------------------------------------------------------------------

#Grab all information from ensembl IDs
ids <- gsub("\\.\\d+", "", as.character(rownames(y)))

edb=EnsDb.Hsapiens.v75
df <- select(edb, keys = ids, keytype = "GENEID", 
             columns = c("GENENAME", "SEQNAME", 
                         "GENESEQSTART","GENESEQEND", "SEQLENGTH", "SEQSTRAND"))

temp = data.frame(round(cpm(y),1))
colnames(temp) = paste0(y$samples$samples,".",y$samples$group)
strand = ifelse(df$SEQSTRAND == -1, "-", "+")

y$genes = data.frame(df[,1:2], paste0("chr", df$SEQNAME), df[,4:6], strand, temp, stringsAsFactors=FALSE)

colnames(y$genes) = c("ENSEMBLID", "GENESYMBOL", 
                      "CHROM", "START", "END", "WIDTH", "STRAND", 
                      paste0(y$samples$samples,".",y$samples$group))

##Statistical testing, fit gene-wise glms ####-------------------------------------------------------------------------
fit <- glmFit(y, design_byP) #blocking 

#save(y,edb,files,files2,cond,design_byP,fit, file="DGE_obj_Block_filtGenes.RData")
#save(y,edb,files,files2,cond,design,f, file="DGE_obj_notBlock.RData")
#save(y,edb,files,files2,cond,design_byP,fit, file="DGE_obj_Block.RData")
#save(y,df,files,files2,cond,design_byP,fit, file="DGE_obj_Block_withVM.RData")
save(y,df,files,files2,cond,design_byP,fit,cond_forplots_fix, 
     age,agegroup,tissue.colors, gender.colors, gender, file="DGE_obj_Block_withVM_moreGenes_novEdit.RData")
#save(y,df,files,files2,cond,design_byP,fit,cond_forplots_fix, 
#     age,agegroup,tissue.colors, gender.colors, gender, file="DGE_obj_Block_withVM_moreGenes_novEdit_ncrnaFIX.RData")

load(file="DGE_obj_Block_withVM_moreGenes_novEdit.RData")

#compare all the pairwise differences between treatments
stub=c("ADENOMA-NORMAL","SSA-NORMAL","SSA-ADENOMA")
mc <- makeContrasts(contrasts=stub, levels=levels(cond))


#Adenoma - normal
lrt1 <- glmLRT(fit, coef=2)
de1 <- decideTestsDGE(lrt1, p.value=0.05, lfc=0.5) 

#plotMD(lrt1)
#abline(h=c(-1, 1), col="blue")

#SSA - normal
lrt2 <- glmLRT(fit, coef=3)
de2 <- decideTestsDGE(lrt2, p.value=0.05, lfc=0.5)

#SSA - Adenoma
lrt3 <- glmLRT(fit, contrast=c(0,-1,1,rep(0,31)))
de3 <- decideTestsDGE(lrt3, p.value=0.05, lfc=0.5)

des_b <- cbind(de1,de2,de3)
des_b <- as.data.frame(des_b)
colnames(des_b) <- c("ADENOMA-NORMAL","SSA-NORMAL","SSA-ADENOMA")
save(lrt1,lrt2,lrt3,de1,de2,de3,des_b, file="DGE_testsBlock_filtGenes_withVM_moreGenes_novEdit.Rdata")
load(file="DGE_testsBlock_filtGenes_withVM_moreGenes_novEdit.Rdata")

##Venn diagram ####---------------------------------------------------------------------------------------

pdf("VennDiagrams.pdf")

vennDiagram(des_b, include=c("up", "down" ), counts.col = c("gray30", "gray50"), 
            cex=1.2, circle.col = c("#38678f","#CD3278","#ffa500"), mar = rep(0,4), lwd = 5)
vennDiagram(des_b, counts.col = "gray30", 
            cex=1.2, circle.col = c("#38678f","#CD3278","#ffa500"), mar = rep(0,4), lwd = 5)

#vennDiagram(des_b[,1:2], include=c("up", "down" ), counts.col = c("gray30", "gray50"), 
#            cex=1.2, circle.col = c("#38678f","#CD3278"), mar = rep(0,4), lwd = 5)
#dev.off()

#Another option
pairVenn = function(n1, n2, counts, name, comp, colors){
  #pdf(paste0("Venn_",name,".pdf"))
  grid.newpage()
  draw.pairwise.venn(area1 = n1, 
                     area2 = n2, 
                     cross.area = counts[4,3], 
                     category = comp, 
                     lty = rep("blank",2), 
                     fill = colors, 
                     alpha = rep(0.5, 2),
                     cex = rep(1.5,3),
                     label.col = rep("gray30", 3),
                     fontfamily = rep("sans", 3),
                     cat.pos = c(-5,5), 
                     cat.dist = rep(0.025, 2),
                     cat.col = rep("gray30", 2),
                     cat.fontfamily = rep("sans", 2))
  #dev.off()
}

par(mfrow=c(2,2))
#Total
countsANSN = vennCounts(des_b[,1:2])
countsANAS = vennCounts(des_b[,c(1,3)])
countsSNAS = vennCounts(des_b[,2:3])
AN = length(which(des_b$`ADENOMA-NORMAL` != 0))
SN = length(which(des_b$`SSA-NORMAL` != 0))
AS = length(which(des_b$`SSA-ADENOMA` != 0))
pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))

#down
countsANSN = vennCounts(des_b[,1:2]) #choose only down genes
AN = length(which(des_b$`ADENOMA-NORMAL` != 0))
SN = length(which(des_b$`SSA-NORMAL` != 0))
pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))

#up
pairVenn(SN, AN, countsANSN, "ANSN", c("SSA/P Vs Normal", "Adenoma Vs Normal"), c("#CD3278", "#38678f"))
#pairVenn(AN, AS, countsANAS, "ANAS", c("Adenoma Vs Normal", "SSA/P Vs Adenoma"), c("#38678f","#ffa500"))
#pairVenn(SN, AS, countsSNAS, "SNAS", c("SSA/P Vs Normal", "SSA/P Vs Adenoma"), c("#CD3278","#ffa500"))

counts = vennCounts(des_b)
grid.newpage()
draw.triple.venn(area1 =AN, 
                 area2 = SN, 
                 area3 = AS, 
                 n12 = countsANSN[4,3], 
                 n23 = countsSNAS[4,3], 
                 n13 = countsANAS[4,3], 
                 n123 = counts[8,4], 
                 category = c("Adenoma Vs Normal", "SSA/P Vs Normal", "SSA/P Vs Adenoma"), 
                 lty = "blank", 
                 fill = c("#38678f","#CD3278","#ffa500"),
                 alpha = 0.5,
                 cex = 1.5,
                 label.col = "gray30",
                 fontfamily = "sans", 
                 cat.dist = 0.05,
                 cat.col = "gray30",
                 cat.fontfamily = "sans",
                 cat.cex = 0.8)
dev.off()

###Smear plot ###----------------------------------------------------------------------------
smear = function(lrt, de, name){
  isde = which(de != 0)
  de.genes <- rownames(lrt$table[isde,])
  pdf(paste0("Smear_plot_",name,".pdf"))
  plotSmear(lrt, smooth.scatter = T)
  plotSmear(lrt, de.tags = de.genes, smooth.scatter = T)
  dev.off()
}

smear(lrt1, de1, "AN")
smear(lrt2, de2, "SN")
smear(lrt3, de3, "SA")

###Volcano plots###-----------------------------------------------------------------------
#AD VS norm
volcanoData <- cbind(lrt1$table$logFC, -log10(lrt1$table$PValue))
colnames(volcanoData) <- c("logFC", "negLogPval")
plot(volcanoData, pch=19)

#ggplot version
volcano = function(lrt, name){
  level = ifelse(lrt$table$logFC >= 1, "Up", "--")
  level[which(lrt$table$logFC <= -1)]  = "Down"
  adjustedp = p.adjust(lrt1$table$PValue,method="BH")
  level[which(adjustedp >= 0.05)] = "--" 
  volcanoData = data.frame(logFC = lrt$table$logFC, negLogPval = -log10(adjustedp), level = level)
  ggplot(volcanoData, aes(x=logFC, y=negLogPval, colour = level)) + geom_point() + 
    theme_bw() +
    geom_vline(xintercept=c(-1, 1), col="black", lty=2, lwd=0.5) +
    geom_hline(yintercept=1.30, col="black", lty=2, lwd=0.5) +
    scale_color_manual(values=c("black", "red", "green")) + 
    coord_cartesian(xlim = c(-10, 10), ylim =c(0,130)) 
  ggsave(paste0("Volcano_plot_",name,".pdf"))
}

volcano(lrt1, "AN")
volcano(lrt2, "SN")
volcano(lrt3, "SA")

#Smear version of volcano
volcano.smooth = function(lrt, name){
  level = ifelse(lrt1$table$logFC >= 1, "Up", "--")
  level[which(lrt1$table$logFC <= -1)]  = "Down"
  adjustedp = p.adjust(lrt1$table$PValue,method="BH")
  level[which(adjustedp >= 0.05)] = "--" 
  volcanoData = data.frame(logFC = lrt1$table$logFC, negLogPval = -log10(adjustedp), level = level)
  
  pdf(paste0("SmoothVolcano_plot_",name,".pdf"))
  Lab.palette <- colorRampPalette(c("white", "orange", "red"), space = "Lab")
  smoothScatter(volcanoData$logFC, volcanoData$negLogPval, nbin = 200,
                colramp = Lab.palette)
  abline(v = c(-1,1), h = 1.30, lty = 2)
  dev.off()
}

volcano.smooth(lrt1, "AN")
volcano.smooth(lrt2, "SN")
volcano.smooth(lrt3, "SA")

##Gene clustering ####-------------------------------------------------------------------------------------------------

#biased - DE genes
isde=des_b != 0
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
ann_colors = list(Tissue = c('NORMAL.ADENOMA'="#AADEDB",'ADENOMA'="#38678f",'NORMAL.SSA/P'="#EAABC8",'SSA/P'="#CD3278"), 
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
zz=as.matrix(log2(1+cps))
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

##Make master table ####-----------------------------------------------------------------------------------------------
res <-data.frame(y$genes)

#with block...too lazy to do a loop
df1=data.frame(lrt1$table, Padj=p.adjust(lrt1$table$PValue,method="BH"), de=de1, stringsAsFactors=FALSE)
colnames(df1) <- paste("Aden-Norm", colnames(df1), sep=".")
colnames(df1)[6] = "Aden-Norm"

df2=data.frame(lrt2$table, Padj=p.adjust(lrt2$table$PValue,method="BH"), de=de2, stringsAsFactors=FALSE)
colnames(df2) <- paste("SSA-Norm", colnames(df2), sep=".")
colnames(df2)[6] = "SSA-Norm"

df3=data.frame(lrt3$table, Padj=p.adjust(lrt3$table$PValue,method="BH"), de=de3, stringsAsFactors=FALSE)
colnames(df3) <- paste("SSA-Aden", colnames(df3), sep=".")
colnames(df3)[6] = "SSA-Aden"
res <- cbind(res, df1,df2,df3)

#Joint Table
pvals <- res[, grep("PValue", colnames(res))]
rm <- rowMins(as.matrix(pvals))
o <- order(rm)
res <- res[o,]
write.table(res, "serrated_table_DGE_moreGenes_lfc0.5_nov1_ncrnaFIX.csv", row.names=F, quote=FALSE, sep="\t")

