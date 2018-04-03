#!/bin/bash -x

#The script runs Salmon for a list of files specified in another file


#Download reference and create index
wget ftp://ftp.ensembl.org/pub/grch37/release-83/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.fa.gz
salmon index -t Homo_sapiens.GRCh37.cdna.fa -i Homo_sapiens.GRCh37.83.cdna.salmon0.8.2.sidx -p 6

#Loop through fasta files
while read line; do
	header=$(echo "$line" | cut -d"." -f1-3)
	name=$(echo "$line" | cut -d"." -f4 | cut -d"_" -f1-4 | sed s/_R1//)
	
	salmon quant -i Homo_sapiens.GRCh37.83.cdna.salmon0.8.2.sidx -l A -1 "$line" -2 "$header"."$name"_R2.fastq.gz -o out_files_cdna83/"$name"_quant -p 6

done < "file_names.txt"




