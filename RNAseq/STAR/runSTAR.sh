#!/bin/bash -x

#The script runs STAR for a list of files specified in another file

#Create reference for seq length of 100
#STAR --runMode genomeGenerate --genomeDir /home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.71/STAR/sjdbOverhang150 --genomeFastaFiles /home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.71/genome/ensembl_Homo_sapiens.GRCh37.71.fa --sjdbGTFfile /home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.71/gtf/Homo_sapiens.GRCh37.71.gtf --runThreadN 10 --sjdbOverhang 150


#done out the script, get forward files and paths from folder, and put them in a list
#~/RNAseq/STAR
#I'll try out first the files in this folder: 
#find ../FASTQs/MAR_2017/*.gz | grep _R1 > File_names
#../FASTQs/All_files/20170622.A-VM_Adenoma2_R1.fastq.gz
while read line; do
	header=$(echo "$line" | cut -d"_" -f1-3)
	#header=$(echo "$line" | cut -d"." -f1-2 | sed s/_R1//)
	#name=$(echo "$line" | cut -d"." -f2 | cut -d"_" -f1,2)
	name=$(echo "$line" | cut -d"." -f4 | cut -d"_" -f1-2)
	#echo "$header"
	#echo "$name"
	#if [ -f out_files/"$name"_STAR.bam ]; then
	#	echo "File exists, this sample was merged"		
#STAR --genomeDir /home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.71/STAR/sjdbOverhang150/ --readFilesIn "$line" "$header"."$name"_R2.fastq.gz --runThreadN 10 --outFilterMultimapNmax 1 --readFilesCommand zcat
	#	samtools view -@ 10 -Sb Aligned.out.sam > out_files/"$name"_2_STAR.bam && rm Aligned.out.sam 
	#	mv Log.final.out out_files/"$name"_2_STAR.bam.log
	#	samtools sort -@ 10 -m 20G -O bam -T _tmp -o out_files/"$name"_2_STAR_s.bam out_files/"$name"_2_STAR.bam
	#	samtools index out_files/"$name"_2_STAR_s.bam
	#else
	STAR --genomeDir /home/Shared_taupo/data/annotation/Human/Ensembl_GRCh37.71/STAR/sjdbOverhang150/ --readFilesIn "$line" "$header"_R2.fastq.gz --runThreadN 6 --outFilterMultimapNmax 1 --readFilesCommand zcat 
	samtools view -@ 6 -Sb Aligned.out.sam > out_files/"$name"_STAR.bam && rm Aligned.out.sam 
	mv Log.final.out out_files/"$name"_STAR.bam.log
	samtools sort -@ 6 -m 20G -O bam -T _tmp -o out_files/"$name"_STAR_s.bam out_files/"$name"_STAR.bam
	samtools index out_files/"$name"_STAR_s.bam
	#fi
done < "File_names2"


#extra stuff
#scp runSTAR.bash sorjuela@imlssherborne:RNAseq/STAR
#chmod +x runSTAR.bash
