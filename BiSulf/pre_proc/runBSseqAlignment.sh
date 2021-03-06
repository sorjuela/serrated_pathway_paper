#!/bin/bash -x

#The script runs the entire methylation extraction pipeline from fastq files (from pre_proc folder)

### Make a genome for alignment: The probes in the SeqCap Epi CpGGiant protocol were designed under the hg19 build

#wget ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.*.fa.gz
#mkdir all_chroms / mv * all_chroms
#cat all_chroms/Homo_sapiens.GRCh37.dna.chromosome.{{1..22},X,Y,MT}.fa.gz > GRCh37.91.fa.gz
#gzip -d GRCh37.91.fa.gz
#/home/sorjuela/bismark_v0.18.1/bismark_genome_preparation --verbose .

### Run alignment and extract methylation
while read line; do
	file2=$(echo "$line" | sed 's/_R1/_R2/' )
	name=$(echo "$line" | egrep -o '[A-Z0-9]{2,4}_[A-Za-z]{3,}')
	header=$(echo "$line" | egrep -o '[0-9]{8}\.[A-Z]{1}\-')  #FASTQs/20160819.A-
	
	#Check quality	
	fastqc -o Trimmomatic/preFASTQC -t 2  "$line" "$file2" 
	
	#Trim

	#TrimGalore: PREFERED
	/home/sorjuela/TrimGalore/TrimGalore-0.4.5/trim_galore --fastqc -o TrimGalore/ --paired "$line" "$file2" --clip_R1 5 --clip_R2 20 --three_prime_clip_R1 5 --three_prime_clip_R2 20 --length 30

	#Trimmomatic
	java -jar /home/sorjuela/trimmomatic-0.35.jar PE -threads 5 -phred33 "$line" "$file2" -baseout Trimmomatic/outfiles/"$name".Tr.fastq.gz LEADING:20 HEADCROP:4 TRAILING:20 SLIDINGWINDOW:4:15 MINLEN:30
	fastqc -o Trimmomatic/FASTQC_outfiles Trimmomatic/outfiles/"$name".Tr_*P.fastq.gz
	
	#run Bismark
	/home/sorjuela/bismark_v0.18.1/bismark -p 10 -o bismark/ -B "$name" --genome /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91 -1 TrimGalore/"$header""$name"_R1_val_1.fq.gz -2 TrimGalore/"$header""$name"_R2_val_2.fq.gz

	#deduplicate
	/home/sorjuela/bismark_v0.18.1/deduplicate_bismark -p --bam bismark/"$name"_pe.bam


	#Extract methylation

	#for some reason path specification does not work at this step, so instead access the folder
	cd bismark

	/home/sorjuela/bismark_v0.18.1/bismark_methylation_extractor -p --bedGraph --comprehensive --genome_folder /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91 --no_overlap --ignore_r2 5 --ignore_3prime_r2 2 --ignore_3prime 2 --multicore 6 -o . "$name"_pe.deduplicated.bam

	#Collapse strand information
	/home/sorjuela/bismark_v0.18.1/coverage2cytosine --dir . -o "$name".c2c --genome_folder /home/Shared_taupo/data/annotation/Human/GRCH37/Bisulfite_Genome.release91 --merge_CpG "$name"_pe.deduplicated.bismark.cov.gz


done < "file_names.txt"
