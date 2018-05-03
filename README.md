### The proto CpG island methylator phenotype of sessile serrated adenoma/polyps
## Code and data analysis

This repository contains the code for the following paper:

* H Parker,S Orjuela, A Martinho Oliveira, F Cereatti, M Sauter, G Tanzi, A Weber, P Komminoth, S Vavricka, L Albanese, F Buffoli, MD Robinson and G Marra: The proto CpG island methylator phenotype of sessile serrated adenoma/polyps. In Preparation (2018).

The structure of this repository is as follows:

* `RNAseq/`
Contains all the code in R and shell to obtain differentially expressed genes (DEGs) from bulk RNAseq data.

* `RNAseq/STAR`
Contains the code to align fastq files using __STAR__, followed by the generation of bigwig and bedgraph files for visualization (for example in IGV). This is in a different folder, since these type of files cannot be generated with alignments generated by Salmon.

* `RNAseq/Salmon>`
Contains the main pipeline to identify DEGs, using __Salmon__ to align the fastq files. Pre-processing of data was done in shell, and further analysis in R.

* `RNAseq/Data>`
Contains the DGEList object in an .RData file, as well as the metadata and the final list of differentially expressed genes for each comparison.

* `BiSulf/`
Contains all the code in R and shell to obtain differentially methylated sites (DMCs) and regions (DMRs) from targeted bisulfite sequencing data. 

* `BiSulf/Data`
Contains 6 lists: DMCs and DMRs identified in both Adenomas and SSA/Ps (4), list of annotated DMRs, and list of  genes with their corresponding methylation and matched expression (for figure 4.D).


