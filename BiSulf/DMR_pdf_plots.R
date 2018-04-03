#########################################################################################
# R script to make plots for DMCs and DMRs
#
# BS-seq data set with 31 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 16 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
# And 6 paired samples with diagnosed colorectal cancer, divided in CIMP (3) and nonCIMP (3) 
# depending on the methylation state of MLH1
# 
# Stephany Orjuela, January 2018
#########################################################################################

setwd("/run/user/1000/gvfs/sftp:host=imlssherborne.uzh.ch/home/sorjuela/serrated_pathway_paper/BiSulf")

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

load("SSA.allclusters.trimmed.RData")
SSAsites <- allClusters.trimmed
load("Adenoma.allclusters.trimmed.RData")
AdenSites <- allClusters.trimmed

###Volcano plot ###---------------------------------------------------------------------------------------
volcano <- function(clusters.trimmed, name){
  DMSites <- data.frame(meth.diff=clusters.trimmed$meth.diff, p.li=clusters.trimmed$p.li)
  level <- ifelse(DMSites$meth.diff >= 0.10, "Hyper", "--")
  level[DMSites$meth.diff <= -0.10]  = "Hypo"
  level[clusters.trimmed$p.li > 0.01] = "--" 
  
  #Add numbers to plot
  numbers <- c()
  numbers[1] <- paste0("n = ",dim(clusters.trimmed[level == "Hypo",])[1])
  numbers[2] <- paste0("n = ",dim(clusters.trimmed[clusters.trimmed$meth.diff > -0.10 & clusters.trimmed$meth.diff < 0.10 & clusters.trimmed$p.li <= 0.01,])[1])
  numbers[3] <- paste0("n = ",dim(clusters.trimmed[level == "Hyper",])[1])
  
  ggplot(DMSites, aes(x=meth.diff, y=-log10(p.li), colour = level)) + 
    geom_point() + 
    theme_bw() +
    labs(x="Methylation difference per CpG site", y = "-log10(P-value)") +
    geom_vline(xintercept=c(-0.1, 0.1), col="black", lty=2, lwd=0.5) +
    geom_hline(yintercept=2, col="black", lty=2, lwd=0.5) +
    scale_color_manual(values=c("black", "yellow", "blue")) +
    coord_cartesian(xlim = c(-0.75, 0.75), ylim =c(0,90)) +
    annotate("text", x = c(-0.6,0,0.6), y = c(75,75,75), label = numbers, size = 3)
  ggsave(paste0("Volcano_plot_",name,"_newrun.jpeg"), device = "jpeg", dpi = 300)
}

volcano(SSAsites, "SN")
volcano(AdenSites, "AN")

### Overlay two histograms ###--------------------------------------------------

r1 <- length(filtSSAsites$chr[filtSSAsites$meth.diff >= 0.10])
r2 <- length(filtAdensites$chr[filtAdensites$meth.diff >= 0.10])
l1 <- paste0("SSA/P = ", length(filtSSAsites$chr[filtSSAsites$meth.diff <= -0.10]))
l2 <- paste0("Adenoma = ", length(filtAdensites$chr[filtAdensites$meth.diff <= -0.10]))
c1 <- length(filtSSAsites$chr[filtSSAsites$meth.diff < 0.10 & filtSSAsites$meth.diff > -0.10])
c2 <- length(filtAdensites$chr[filtAdensites$meth.diff < 0.10 & filtAdensites$meth.diff > -0.10])

d <- data.frame(MethDiff = c(filtSSAsites$meth.diff, filtAdensites$meth.diff), 
               Tissue=rep(c("SSA/P", "Adenoma"), c(length(filtSSAsites$meth.diff), length(filtAdensites$meth.diff))))
p1 <- ggplot(d) + 
  geom_density(aes(x=MethDiff, colour=Tissue, fill=Tissue), alpha=0.5) + 
  labs(x="Methylation difference per CpG site", y="Density") + 
  theme_light(base_size = 15) + 
  scale_color_manual(values=c("#38678f","#CD3278")) + 
  scale_fill_manual(values=c("#38678f","#CD3278")) +
  geom_vline(xintercept=c(-0.1,0.1), col="black", lty=2, lwd=0.5) +
  coord_cartesian(xlim = c(-0.75, 0.75), ylim =c(0,4.5)) +
  #geom_rect(aes(xmin=-0.10, xmax=0.10, ymin=0, ymax=5), fill = "gray", alpha=0.01) +
  annotate("text", x = c(0.5,0.5,-0.5,-0.5,0,0), y = c(4,4.2,4,4.2,4,4.2), label = c(r1,r2,l1,l2,c1,c2), size = 3)
p1  
ggsave("overlayDensity_newrun.jpeg", device = "jpeg")
ggsave("overlayDensity_newrun_filt.pdf")

#### Density plot for CRC samples ####---------------------------------------------------

load("CIMP.allclusters.trimmed.RData")
cimpsites <- allClusters.trimmed[allClusters.trimmed$p.li <= 0.01,]
load("nonCIMP.allclusters.trimmed.RData")
nonsites <- allClusters.trimmed[allClusters.trimmed$p.li <= 0.01,]

r1 <- length(cimpsites$chr[cimpsites$meth.diff >= 0.10])
r2 <- length(nonsites$chr[nonsites$meth.diff >= 0.10])
l1 <- paste0("CIMP = ", length(cimpsites$chr[cimpsites$meth.diff <= -0.10]))
l2 <- paste0("Non-CIMP = ", length(nonsites$chr[nonsites$meth.diff <= -0.10]))
c1 <- length(cimpsites$chr[cimpsites$meth.diff < 0.10 & cimpsites$meth.diff > -0.10])
c2 <- length(nonsites$chr[nonsites$meth.diff < 0.10 & nonsites$meth.diff > -0.10])

d <- data.frame(MethDiff = c(cimpsites$meth.diff, nonsites$meth.diff), 
                Tissue=rep(c("CIMP", "nonCIMP"), c(length(cimpsites$meth.diff), length(nonsites$meth.diff))))
p1 <- ggplot(d) + 
  geom_density(aes(x=MethDiff, colour=Tissue, fill=Tissue), alpha=0.5) + 
  labs(x="Methylation difference per CpG site", y="Density") + 
  theme_light(base_size = 15) + 
  scale_color_manual(values=c("#32cd87","#8f6038")) + 
  scale_fill_manual(values=c("#32cd87","#8f6038")) +
  geom_vline(xintercept=c(-0.1,0.1), col="black", lty=2, lwd=0.5) +
  coord_cartesian(xlim = c(-0.75, 0.75), ylim =c(0,4.5)) +
  annotate("text", x = c(0.5,0.5,-0.5,-0.5,0,0), y = c(4,4.2,4,4.2,4,4.2), label = c(r1,r2,l1,l2,c1,c2), size = 3)
p1  
ggsave("overlayDensity_CANCERS.jpeg", device = "jpeg")
ggsave("overlayDensity_CANCERS.pdf")

#### Upset plot ####-------------------------------------------------------------------

refiltSSAsites <- filtSSAsites[abs(filtSSAsites$meth.diff) >= 0.10,]
refiltAdensites <- filtAdensites[abs(filtAdensites$meth.diff) >= 0.10,]

DMSitesSN <- GRanges(seqnames=refiltSSAsites$chr,
                    ranges=IRanges(start=refiltSSAsites$pos, width=1),
                    hyper=ifelse(refiltSSAsites$meth.diff < 0, 0,1),
                    hypo=ifelse(refiltSSAsites$meth.diff < 0, 1,0))

DMSitesAN <- GRanges(seqnames=refiltAdensites$chr,
                    ranges=IRanges(start=refiltAdensites$pos, width=1),
                    hyper=ifelse(refiltAdensites$meth.diff < 0,0,1),
                    hypo=ifelse(refiltAdensites$meth.diff < 0,1,0))

#Find common sites and make empty matrix
fData <- unique(c(DMSitesSN, DMSitesAN))
hyper.diff <- matrix(integer(length = length(fData) * 2), nrow=length(fData))
hypo.diff = matrix(integer(length = length(fData) * 2), nrow=length(fData))

sites_comp = list(DMSitesAN, DMSitesSN)

#fill in matrix with 1s and 0s, and meth diff values
for(i in seq(along=sites_comp)){
  mtch <- findOverlaps(fData, sites_comp[[i]])
  ind1 <- queryHits(mtch)
  ind2 <- subjectHits(mtch)
  hyper.diff[ind1, i] <- mcols(sites_comp[[i]])$hyper[ind2]
  #hyper.diff[ind1, i + 2] <- mcols(sites_comp[[i]])$hyper.val[ind2]
  hypo.diff[ind1, i] <- mcols(sites_comp[[i]])$hypo[ind2]
  #hypo.diff[ind1, i + 2] <- mcols(sites_comp[[i]])$hypo.val[ind2]
}

colnames(hyper.diff)=c("Adenomas", "SSA/Ps")#, "NA-val", "NS-val")#,"ADENOMA-SSA")
colnames(hypo.diff)=c("Adenomas", "SSA/Ps")#, "NA-val", "NS-val")#,"ADENOMA-SSA")

#Make upset

pdf("upset.plots.newrun_filt.pdf")
upset(as.data.frame(hypo.diff), 
      sets = c("Adenomas","SSA/Ps"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMCs (hypomethylated)",
      main.bar.color = "#3b56d8",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 155000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMCs")


upset(as.data.frame(hyper.diff), 
      sets = c("Adenomas","SSA/Ps"),
      keep.order = T,
      mb.ratio = c(0.8, 0.2), 
      mainbar.y.label = "No. of DMCs (hypermethylated)",
      main.bar.color = "#e1cf22",
      point.size = 3,
      line.size = 1,
      matrix.color = "black",
      mainbar.y.max = 155000, 
      sets.bar.color = "gray23",
      text.scale = c(2,1.5,1.5,1,2,2.5),
      sets.x.label = "Total number of DMCs")
dev.off()



#### DMRs plots ###-----------------------------------------------------------
load(file = "SSAvsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")
SSADMRsTable <- as(SSADMRs, "data.frame")

load(file = "AdenvsNorm.DMRs.RData")
AdenDMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25]
mcols(AdenDMRs)$state <- ifelse(AdenDMRs$median.meth.diff >= 0, "hyper", "hypo")
AdenDMRsTable <- as(AdenDMRs, "data.frame")

#To make another plot with only the overlapping ones
# over <- findOverlaps(SSADMRs, AdenDMRs)
# 
# commonAdenDMRs <- as(AdenDMRs[subjectHits(over)], "data.frame")
# commonSSADMRs <- as(SSADMRs[queryHits(over)],"data.frame")

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
  #geom_boxplot(aes(x = state, y=width) ,fill = "white", colour = "#CD3278")

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

#you need all points even the ones without the filter
load(file = "SSAvsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[biseqDMRs$percentageOverlap >= 25]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")
SSADMRsTable <- as(SSADMRs, "data.frame")

load(file = "AdenvsNorm.DMRs.RData")
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
  ggsave(paste0("Volcano_plot_",name,"DMRs.jpeg"), device = "jpeg", dpi = 300)
}

volcano(SSADMRsTable, "SN")
volcano(AdenDMRsTable, "AN")

#### Upset plots DMRs ####----------------------------------------------------------

#These have the full filters

load(file = "SSAvsNorm.DMRs.RData")
SSADMRs <- biseqDMRs[abs(biseqDMRs$median.meth.diff) >= 0.10 & biseqDMRs$percentageOverlap >= 25 & biseqDMRs$median.p <= 0.01]
mcols(SSADMRs)$state <- ifelse(SSADMRs$median.meth.diff >= 0, "hyper", "hypo")

load(file = "AdenvsNorm.DMRs.RData")
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
  #return(hyperm)
  l <- list(uncomA, uncomS)
  #return(uncomS)
}

hypoM <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
hyperM <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)

#Fix to get these stupid numbers:
uniqHypoAs <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
uniqHyperSs <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)
findOverlaps(uniqHyperSs, uniqHypoAs) #276

### Get number of DMRs per 4kb window and count for each comparison. Exclude ones overlapping
annoData <- toGRanges(EnsDb.Hsapiens.v75)
proms <- promoters(annoData, upstream = 2000, downstream = 2000) #64102

#Get all unique ranges for hypo
hypoM <- makeUpsetTable("hypo","hypo", SSADMRs, AdenDMRs)
hyperM <- makeUpsetTable("hyper","hyper", SSADMRs, AdenDMRs)

allM <- c(hypoM, hyperM)

mat <- matrix(0,length(proms), 4)

for(i in 1:4){
  hits <- findOverlaps(proms, allM[[i]])
  mat[rle(queryHits(hits))$values,i] <- rle(queryHits(hits))$lengths
}

mcols(proms) <- mat
colnames(mcols(proms)) <- c("unHypoA", "unHypoS", "unHyperA", "unHyperS")
proms <- proms[rowSums(mat) > 0]
gene_names <- annoData$gene_name[rowSums(mat) > 0]
proms$gene_name <- gene_names

proms[mcols(proms)[,1] > 0 & mcols(proms)[,2] > 0] # 6 close hypos
proms[mcols(proms)[,3] > 0 & mcols(proms)[,4] > 0] # 29 close hyper
proms[mcols(proms)[,1] > 0 & mcols(proms)[,4] > 0] #122 close hypoA and hyperS

promdf <- as(proms, "data.frame")
write.table(promdf, "closeDMRCount.4kbwindowFromTSS.csv", row.names=F, quote=FALSE, sep="\t")

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

#Annotation with proximity to gene body
#several DMRs per gene
#not several genes per DMR

quickanno <- function(){
annoM <- ChIPpeakAnno::annotatePeakInBatch(SSADMRs,
                                           FeatureLocForDistance="start",
                                           AnnotationData = annoData,
                                           output = "nearestLocation",
                                           PeakLocForDistance= "middle")




#function from other script - DMRannotation
#Fun to choose clostes DMR per gene
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

###Annotation cakes with MethylKit/Genomation packages ####------------------

#TRY TO FIX this track like the refseq track cause It's getting on my nerves
ensgenes <- read.table("genes2Ensembl.bed", stringsAsFactors = FALSE)

V11 <- mapply(function(x,y){
  temp <- as.integer(strsplit(x, ",")[[1]])
  res <- temp - y  
  piece <- paste(res, collapse=",")
  paste0(piece, ",")
  }, x = ensgenes$V11, y = ensgenes$V2, USE.NAMES = F)

V12 <- mapply(function(x,y){
  temp <- as.integer(strsplit(x, ",")[[1]])
  res <- temp - y  
  piece <- paste(res, collapse=",")
  paste0(piece, ",")
}, x = ensgenes$V12, y = ensgenes$V2, USE.NAMES = F)

#WRITE.TABLE

ensgenes2 <- cbind(ensgenes[,1:10],V12,V11)
write.table(ensgenes2, file = "genes2Ensembl.Rver.bed", row.names=F, quote=FALSE, sep="\t", col.names=F)
gene.obj <- readTranscriptFeatures("genes2Ensembl.Rver.bed", up.flank=1000,down.flank=500)
#I give up, dont know what sorcery this function does

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

####Gviz tracks ####------------------------------------------------------------------------------
#HUNK
#chr21:33,245,611-33,248,500 zoom in
afrom <- 33245611
ato <- 33250234
window <- GRanges("chr21", IRanges(afrom, ato))

##meth part start, long intron without anything, expr part

#chr21:33,248,515-33,378,000 zoom out
afrom <- 33340560
ato <- 33378000
window <- GRanges("chr21", IRanges(afrom, ato))

#ZIC2
afrom <- 100632207
ato <- 100644719
window <- GRanges("chr13", IRanges(afrom, ato))

#ideoTrack <- IdeogramTrack(genome = "hg19", chromosome = "chr13")
gtrack <- GenomeAxisTrack()

#GeneTrack
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75
gr <- getGeneRegionTrackForGviz(edb, chromosome = "chr13",
                                start = afrom, end = ato)
gr <- gr[11:15] #for ZIC2
gr <- gr[c(6,8:18)] #for HUNK zoom out

geneTrack <- GeneRegionTrack(gr, collapseTranscripts = F, fill = "#00cd00")

#
#AnnotationTrack -------
#CpG islands
islands <- read.csv("probes/CpGIslandTrack.bed", header =F, sep="\t")
islandR <- GRanges(seqnames = islands$V1, ranges = IRanges(islands$V2, islands$V3)) #28691
hits <- findOverlaps(window, islandR)
islandR <- islandR[subjectHits(hits)]

#probes
probes <- read.table(file="probes/130912_HG19_CpGiant_4M_EPI.bed")
probes <- GRanges(probes$V1, IRanges(start = probes$V2, end = probes$V3))
hits <- findOverlaps(window, probes)
probes <- probes[subjectHits(hits)]


isTrack <- AnnotationTrack(islandR, name = "CpG Islands", fill = "#0000cc")
pTrack <- AnnotationTrack(probes, name = "Probes", fill = "#FFA500")


#DataTrack------

#methylation------
load("chromBSraw/chr13BSraw.RData")
load("chromBSraw/chr21BSraw.RData")
#load("chromPredictedmeth/predictedMethchr13.RData")

#Have to make a GRanges

#Grab only SSAs
s <- which(colData(BSrchr)$group == "SSA")
n <- s-1

#Grab only Adenomas
s <- which(colData(BSrchr)$group == "ADENOMA")
n <- s+1

BSred <- BSrchr[,c(s,n)]
BSred <- BSrchr[,c(n,s)]

#BSred <- predictedMeth[,c(s,n)]

#Get a ranges and overlap with window
gr <- rowRanges(BSred)
hits <- findOverlaps(window, gr)

#grab stuff from whatever hits the window
c <- totalReads(BSred)[subjectHits(hits),]
m <- methReads(BSred)[subjectHits(hits),]
grwindow <- gr[subjectHits(hits)]
p <- m / c
#or
#p <- methLevel(BSred)[subjectHits(hits),]

mcols(grwindow) <- p #785...447 in predicted meth
q <- rowSums(c >= 5) >= 28 #34 
grwindow2 <- grwindow[q] #581
#grwindow2 <- grwindow #447

dTrack <- DataTrack(grwindow2, name = "Methylated reads proportion",
                    #groups = c(rep("SSA/P", 17), rep("Normal", 17)),
                    groups = c(rep("Adenoma", 14), rep("Normal", 14)),
                    #type = c("confint", "smooth"),
                    type = c("p"),
                    aggregateGroups = TRUE, #the mean
                    #cex.legend = 0.5,
                    ylim = c(0,0.8),
                    #col = c("#F5D6E4", "#CD3278"),
                    col = c("#B9C9D7","#38678f"),
                    legend = F)

#expression-----
#options(ucscChromosomeNames=FALSE)
path <- "../RNAseq/HUNK/"
path <- "../RNAseq/ZIC2/"

#run bedtools intersect before coming here
#remove files with nothing and add columns with zeros to the final matrix
#ssa
bgFile <- list.files(path = path, pattern = "SSA.bed")
#names <- gsub(".bed","", bgFile)

#Aden
bgFile <- list.files(path = path, pattern = "ADENOMA.bed")
names <- gsub(".ADENOMA.bed","", bgFile)
bgFile <- unlist(sapply(names, function(i){list.files(path = path, pattern = paste0("^",i))}))[-c(6,7)]

#bgFile <- unlist(sapply(names, function(i){list.files(path = path, pattern = i)}))

bgFile <- paste0(path, bgFile)
info <- file.info(bgFile)
empty <- which(info$size == 0)
bgFileR <- bgFile[-empty]

#namesR <- names[-empty]
#namesE <- names[empty]

#Aden
bgFileE <- bgFile[empty]


beds <- sapply(bgFileR, function(files){ #change bgFile
    bed <- scan(files, skip=0, sep="\t", comment.char = "#",
                what=list(NULL,NULL,NULL, character(), integer(), integer(), numeric()))
    GRanges(paste0("chr",bed[[4]]), IRanges(start = bed[[5]], end = bed[[6]]), reads = bed[[7]])
  }  )

fullwindowGR <- GRanges("chr13", #change chrom
                        IRanges(start = afrom:ato, width  = 1))

matfull <- matrix(0,length(fullwindowGR),length(bgFile)) 

for(i in 1:length(bgFileR)){ #change bgFile
  hits <- findOverlaps(fullwindowGR, beds[[i]])
  matfull[queryHits(hits),i] <- mcols(beds[[i]][subjectHits(hits)])$reads
}

mcols(fullwindowGR) <- matfull
#colnames(mcols(fullwindowGR)) <- c(namesR, namesE)
#colnames(mcols(fullwindowGR)) <- c(bgFileR, bgFileE)
colnames(mcols(fullwindowGR)) <- bgFile

groups <- rep("SSA/P", 34)

groups <- rep("Adenoma", 28)
groups[grep("NORM", c(bgFileR, bgFileE))] <- "Normal"
groups[grep("NORM", bgFile)] <- "Normal"


#finally build track
dTrack2 <- DataTrack(range = fullwindowGR, 
                     type = "a", 
                     #name = "Raw expression read counts",
                     groups = groups,
                     #col = c("#F5D6E4", "#CD3278"),
                     col = c("#38678f", "#B9C9D7"), #for adenoma cols
                     ylim = c(0,80),
                     legend = F)
                     #aggregateGroups = F) #the mean
              
#Plot gviz----
#pdf("gviz.ZIC2.pdf", width=400,height=200)
#options(ucscChromosomeNames=TRUE)
plotTracks(list(gtrack, dTrack,pTrack,isTrack, dTrack2, geneTrack), 
           from = afrom, 
           to = ato,
           background.title = "white",
           col.axis="darkgray"
           #sizes = c(0.5,1,5,2,0.5,0.5),
           )

plotTracks(list(geneTrack),
  from = afrom, 
  to = ato
)
