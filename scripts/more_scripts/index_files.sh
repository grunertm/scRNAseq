#!/bin/bash

for i in *_Aligned.Isl1_3UTR.bam;
do
	sample=$(echo "$i" | sed "s/\.Isl1\_3UTR\.bam//g")
	
	echo $sample;
	#echo $i

	java -jar ../../tools/picard.jar AddOrReplaceReadGroups I=$i O=$sample\.Isl1_3UTR.rg_addd.bam SO=coordinate RGID=id RGLB=library RGPL=platform RGPU=machine RGSM=sample
	../../tools/samtools-1.9/samtools index $sample\.Isl1_3UTR.rg_addd.bam
	
done

