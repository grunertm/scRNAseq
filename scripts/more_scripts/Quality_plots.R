batch1 <- read.delim(file = "FASTQ_batch1/OUT/multiqc_data/multiqc_fastqc.txt", as.is = T, sep="\t", header=T)
batch2 <- read.delim(file = "FASTQ_batch2/OUT/multiqc_data/multiqc_fastqc.txt", as.is = T, sep="\t", header=T)

batch1_R1 <- batch1[grep("_R1_", batch1$Sample),]
batch2_R1 <- batch2[grep("_R1_", batch2$Sample),]

#paired-end reads
E8_p1 <- batch1_R1[grep("E8_1",batch1_R1$Sample),]
E8_p2 <- batch1_R1[grep("E8_2",batch1_R1$Sample),]
E9_ant_p1 <- batch1_R1[grep("E9_anterior_1",batch1_R1$Sample),]
E9_ant_p2 <- batch2_R1[grep("E9_anterior_2",batch2_R1$Sample),]
E9_post_p3 <-batch2_R1[grep("E9_posterior_3",batch2_R1$Sample),]
E9_post_p4 <- batch1_R1[grep("E9_posterior_4",batch1_R1$Sample),]

jpeg("read.counts.jpg", width = 680, height = 480)
boxplot(E8_p1$Total.Sequences,
        E8_p2$Total.Sequences, 
        E9_ant_p1$Total.Sequences,
        E9_ant_p2$Total.Sequences,
        E9_post_p3$Total.Sequences,
        E9_post_p4$Total.Sequences,
        names = c("E8.5 (#1)", "E8.5 (#2)", "E9.5 ant (#1)", "E9.5 ant (#2)", "E9.5 post (#3)", "E9.5 post (#4)"),
        main="Reads per cell/well",
        ylab="Number of read pairs")

points(1,mean(E8_p1$Total.Sequences), pch=18, col="red")
points(2,mean(E8_p2$Total.Sequences), pch=18, col="red")
points(3,mean(E9_ant_p1$Total.Sequences), pch=18, col="red")
points(4,mean(E9_ant_p2$Total.Sequences), pch=18, col="red")
points(5,mean(E9_post_p3$Total.Sequences), pch=18, col="red")
points(6,mean(E9_post_p4$Total.Sequences), pch=18, col="red")
legend("topleft", "mean", pch=18, col="red")
dev.off()

