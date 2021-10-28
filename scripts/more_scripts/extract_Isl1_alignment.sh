#!/bin/bash

for i in *_Aligned.out.bam;
do
	sample=$(echo "$i" | sed "s/\.out\.bam//g")
	
	echo $sample;
	#echo $i
	
	../tools/samtools-1.9/samtools sort $i -o $sample\.sorted.bam
	../tools/samtools-1.9/samtools index $sample\.sorted.bam
	../tools/samtools-1.9/samtools view -h $sample\.sorted.bam "chr13:116298281-116298572" > $sample\.Isl1_3UTR.bam
	
done

