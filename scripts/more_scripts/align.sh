#!/bin/bash

for i in *R1*fastq.gz;
do
	sample=$(echo "$i" | sed "s/_R1\_001\.fastq\.gz//g")
	
	echo $sample;

	#read mapping usig STAR
	# !!! CHANGE BATCH number & overhang !!!
	/project/NAS-data/scRNA-seq/tools/STAR-2.6.0a/bin/Linux_x86_64/STAR --genomeDir /project/NAS-data/scRNA-seq/references/ucsc/star_indices_overhang150 --sjdbGTFfile /project/NAS-data/scRNA-seq/annotations/gencode.vM20.annotation.gtf --readFilesIn /project/NAS-data/scRNA-seq/raw_data/FASTQ_batch2/$sample* --readFilesCommand zcat --outSAMtype BAM Unsorted --outSAMunmapped Within --quantMode TranscriptomeSAM --twopassMode Basic --runThreadN 24 --outFileNamePrefix /project/NAS-data/scRNA-seq/mapping/batch2/$sample\_

	#calculate gene and transcipt expression levels using RSEM, which uses as input genomic alignments in transcriptomic coordiantes (i.e. Aligned.toTranscriptome.out.bam)
	#forward-prob: probability of generating a read from the forward strand of a transcript. 1: strand-specific protocol where all (upstream) reads are derived from the forward strand; 0: strand-specific protocol where all (upstream) read are derived from the reverse strand; 0.5: non-strand-specific protocol.
	# !!! CHANGE BATCH number !!!
	/project/NAS-data/scRNA-seq/tools/RSEM-1.3.1/rsem-calculate-expression --bam --no-bam-output -p 24 --paired-end --forward-prob 1 /project/NAS-data/scRNA-seq/mapping/batch2/$sample\_Aligned.toTranscriptome.out.bam /project/NAS-data/scRNA-seq/GENOME_data/rsem/rsem_mm10 /project/NAS-data/scRNA-seq/gene_counts/batch2/$sample
	
done


