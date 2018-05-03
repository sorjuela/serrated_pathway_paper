#########################################################################################
# R script to run a normal DGE analysis workflow starting from Salmon counts import
# ending with export of gene tables
#
# RNA-seq data set with 33 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesion with normal tissue, and 16 adenoma lesion with normal tissue
# 
# part of code taken from Charlotte Soneson
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

files <- read.table("Data/new_exp_setup.csv", stringsAsFactors = F, header = T) 

##Make tx2gene table ####-------------------------------------------------------------------------------------------------

## Read the fasta file in R using a function from the Biostrings package
cdna <- readDNAStringSet("Homo_sapiens.GRCh37.cdna.all.fa")

## Go through the sequence names and extract the required information
tx2gene <- data.frame(t(sapply(names(cdna), function(nm) {
  tmp <- strsplit(nm, " ")[[1]]
  tx <- tmp[1]
  gene <- gsub("gene:", "", tmp[grep("^gene:", tmp)])
  c(tx = tx, gene = gene)
})), stringsAsFactors = FALSE)
rownames(tx2gene) <- NULL

##Process to input edgeR, normalize ####---------------------------------------------------------------------------------

txi_salmonsimplesum <- tximport(files = files$FILES, type = "salmon", tx2gene = tx2gene, dropInfReps = TRUE)
#save(txi_salmonsimplesum, file="txi_cdna.RData")

#Calculate offset
cts <- txi_salmonsimplesum$counts
normMat <- txi_salmonsimplesum$length 
normMat <- normMat/exp(rowMeans(log(normMat)))
o <- log(calcNormFactors(cts/normMat)) + log(colSums(cts/normMat))

cond <- gsub("_SSA", "", files$TISSUE)
y <- DGEList(cts, samples = files$PATIENT, group = cond)
y$offset <- t(t(log(normMat)) + o) 

y <- calcNormFactors(y)

#Filter lowly expressed genes
y <- y[which(rowSums(cpm(y) > 0.5) > 1), , keep.lib.sizes = FALSE] 
y <- calcNormFactors(y)

##Pre-Plots ####---------------------------------------------------------------------------------------------------------

#Set up factors for plots
cond_forplots <- factor(files$TISSUE)
cond_forplots_fix <- gsub("NORMAL_SSA", "NORMAL.SSA", cond_forplots)
cond_forplots_fix <- gsub("SSA$", "SSA/P", cond_forplots_fix)
cond_forplots_fix <- factor(gsub("NORMAL$", "NORMAL.ADENOMA", cond_forplots_fix), 
                           levels = c("NORMAL.ADENOMA", "ADENOMA", "NORMAL.SSA/P", "SSA/P"))
gender <- factor(files$SEX)
age <- as.character(files$AGE)
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
m1 <- plotMDS(y)$cmdscale.out
df <- data.frame(dim1 = m1[,1], dim2 = m1[,2])
label <- y$samples$samples

ggplot(df)+
  scale_shape_identity() +
  geom_point(data=df, aes(x=dim1, y=dim2, shape = gender, colour = cond_forplots_fix), size=5) +
  geom_text(data=, aes(x=dim1, y=dim2-0.1, label = label), size=4)+
  scale_colour_manual(values=tissue.colors)+
  scale_shape_manual(values=c(19,17))+
  labs(x = "1", y = "2", title = "")+
  labs(shape = "Gender", colour = "Tissue")+
  theme_bw()
ggsave("MDS2D_cdna.pdf")

##Set up design matrix for GLM ####--------------------------------------------------------------------------------------
#this is an additive linear model with patient as the blocking factor

cond <- relevel(cond, "NORMAL")

design_byP <- model.matrix(~cond+patient) ##block by patient
colnames(design_byP) <- c(levels(cond), levels(y$samples$samples)[-1])

#GLM estimates of dispersion
yg <- estimateGLMCommonDisp(y, design_byP)
yg <- estimateGLMTrendedDisp(yg, design_byP) 
yg <- estimateGLMTagwiseDisp(yg, design_byP)

##Grab gene identifiers ####-------------------------------------------------------------------------------------------

#Grab all information from ensembl IDs
ids <- gsub("\\.\\d+", "", as.character(rownames(yg)))

edb <- EnsDb.Hsapiens.v75
df <- select(edb, keys = ids, keytype = "GENEID", 
             columns = c("GENENAME", "SEQNAME", 
                         "GENESEQSTART","GENESEQEND", "SEQLENGTH", "SEQSTRAND"))

temp <- data.frame(round(cpm(y),1))
colnames(temp) <- paste0(y$samples$samples,".",y$samples$group)
strand <- ifelse(df$SEQSTRAND == -1, "-", "+")

yg$genes <- data.frame(df[,1:2], paste0("chr", df$SEQNAME), df[,4:6], strand, temp, stringsAsFactors=FALSE)

colnames(yg$genes) <- c("ENSEMBLID", "GENESYMBOL", 
                      "CHROM", "START", "END", "WIDTH", "STRAND", 
                      paste0(y$samples$samples,".",y$samples$group))

##Statistical testing, fit gene-wise glms ####-------------------------------------------------------------------------
fit <- glmFit(yg, design_byP)

save(yg,files,cond,design_byP,fit,cond_forplots_fix, 
     age,agegroup,tissue.colors, gender.colors, gender, file="DGE_obj.RData")


#compare all the pairwise differences between treatments

#A-N and S-N
lrt <- list()
de <- matrix(NA, dim(yg$counts)[1], 3)

for(i in 1:2){
  lrt[[i]] <- glmLRT(fit, coef=i+1)
  de[,i] <- decideTestsDGE(lrt[[i]], p.value=0.05, lfc=1) 
}

#And for S-A
lrt[[3]] <- glmLRT(fit, contrast=c(0,-1,1,rep(0,31)))
de[,3] <- decideTestsDGE(lrt[[3]], p.value=0.05, lfc=1)

des_b <- as.data.frame(de)
colnames(des_b)=c("ADENOMA-NORMAL","SSA-NORMAL","SSA-ADENOMA")
save(lrt, de, des_b, file="DGE_tests.Rdata")

##Make master table ####-----------------------------------------------------------------------------------------------
res <- data.frame(yg$genes)

names <- c("Aden-Norm", "SSA-Norm", "SSA-Aden")
for(i in 1:3){
  df <- data.frame(lrt[[i]]$table, 
                   Padj=p.adjust(lrt[[i]]$table$PValue,method="BH"), 
                   de=de[,i], 
                   stringsAsFactors=FALSE)
  colnames(df) <- paste(names[i], colnames(df), sep=".")
  colnames(df)[6] <- names[i] 
  res <- cbind(res, df)
}

#sort table
library(matrixStats)
pvals <- res[, grep("PValue", colnames(res))]
rm <- rowMins(as.matrix(pvals))
o <- order(rm)
res <- res[o,]
write.table(res, "serrated_table_DGE.csv", row.names=F, quote=FALSE, sep="\t")

