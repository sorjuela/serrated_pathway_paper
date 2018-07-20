#########################################################################################
# R script to make plots for DMCs
#
# BS-seq data set with 32 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
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
library(EnsDb.Hsapiens.v75)
library(Gviz)
library(cowplot)

load("SSAVsNorm.DMCs.RData")
SSAsites <- allClusters.trimmed
load("AdenVsNorm.DMCs.RData")
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

### Density plot for lesions ###--------------------------------------------------

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
hypo.diff <- matrix(integer(length = length(fData) * 2), nrow=length(fData))

sites_comp = list(DMSitesAN, DMSitesSN)

#fill in matrix with 1s and 0s, and meth diff values
for(i in seq(along=sites_comp)){
  mtch <- findOverlaps(fData, sites_comp[[i]])
  ind1 <- queryHits(mtch)
  ind2 <- subjectHits(mtch)
  hyper.diff[ind1, i] <- mcols(sites_comp[[i]])$hyper[ind2]
  hypo.diff[ind1, i] <- mcols(sites_comp[[i]])$hypo[ind2]
}

colnames(hyper.diff)=c("Adenomas", "SSA/Ps")
colnames(hypo.diff)=c("Adenomas", "SSA/Ps")

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


