#!/bin/bash

for i in *_Aligned.out.bam;
do
	#sample=$(echo "$i" | sed "s/_R1\_001\.fastq\.gz//g")
	#echo $sample;
	
	echo $i
	
	bedtools coverage -abam $i -b ../Isl1/Isl1_target_file.bed -d -counts > ../Isl1/$i\_Isl1_coverage.txt	
	
done
