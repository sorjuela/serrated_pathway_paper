#!/bin/bash -x

#The script runs STAR for a list of files specified in another file

#Create reference for seq length of 100
STAR --runMode genomeGenerate --genomeDir sjdbOverhang150/ --genomeFastaFiles Homo_sapiens.GRCh37.cdna.fa --sjdbGTFfile Homo_sapiens.GRCh37.gtf --runThreadN 10 --sjdbOverhang 150

while read line; do
	header=$(echo "$line" | cut -d"_" -f1-3)
	name=$(echo "$line" | cut -d"." -f4 | cut -d"_" -f1-2)

	STAR --genomeDir sjdbOverhang150/ --readFilesIn "$line" "$header"_R2.fastq.gz --runThreadN 6 --outFilterMultimapNmax 1 --readFilesCommand zcat 
	samtools view -@ 6 -Sb Aligned.out.sam > out_files/"$name"_STAR.bam && rm Aligned.out.sam 
	mv Log.final.out out_files/"$name"_STAR.bam.log
	samtools sort -@ 6 -m 20G -O bam -T _tmp -o out_files/"$name"_STAR_s.bam out_files/"$name"_STAR.bam
	samtools index out_files/"$name"_STAR_s.bam

done < "file_names.txt"

