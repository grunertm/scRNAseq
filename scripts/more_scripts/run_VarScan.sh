#!/bin/bash

for i in *.rg_addd.bam;
do
	sample=$(echo "$i" | sed "s/\.Isl1\_3UTR\.bam//g")
	
	echo $sample;
	
	samtools mpileup -f ../../../references/ucsc/mm10.fa $sample > $sample\.mpileup
	java -jar ../../../tools/VarScan.v2.3.9.jar pileup2snp $sample\.mpileup --min-coverage 2 --min-reads2 1 > $sample\.SNVs
	java -jar ../../../tools/VarScan.v2.3.9.jar pileup2indel $sample\.mpileup --min-coverage 2 --min-reads2 1 > $sample\.INDELs
	
done

