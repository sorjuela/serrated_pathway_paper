#########################################################################################
# R script to build Gviz tracks for RNAseq counts and Methylation counts
#
# BS-seq data set with 31 paired samples (collaboration with Hannah Parker and Giancarlo Marra): 
# 17 SSA/Ps lesions with normal tissue, and 15 adenoma lesions with normal tissue
# And 6 paired samples with diagnosed colorectal cancer, divided in CIMP (3) and nonCIMP (3) 
# depending on the methylation state of MLH1
# 
# Stephany Orjuela, January 2018
#########################################################################################


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
