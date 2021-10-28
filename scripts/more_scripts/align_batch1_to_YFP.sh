#!/bin/bash

for i in *R1*fastq.gz;
do
	sample=$(echo "$i" | sed "s/_R1\_001\.fastq\.gz//g")
	
	echo $sample;

	#read mapping using bwa
	bwa mem -M -t 14 /project/NAS-data/scRNA-seq/references/YFP.fa /project/NAS-data/scRNA-seq/raw_data/FASTQ_batch1/$sample* > /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_aligned.sam

	#convert sam to bam
	samtools view -Sb /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_aligned.sam > /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_aligned.bam
	
	#statistics
	samtools flagstat /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_aligned.bam > /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_stats.txt
	bam_stat.py -i /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_aligned.bam 2>&1 | tee /project/NAS-data/scRNA-seq/mapping/YFP/batch1/$sample\_stats2.txt

done


