data <- read.delim(file="multiqc_fastqc.txt", as.is = T, sep="\t",header = T)

data_R1 <- data[grep("_R1_", data$Sample),]

data_R1_E8_1 <- data[grep("E8_1", data$Sample),]
data_R1_E8_2 <- data[grep("E8_2", data$Sample),]
data_R1_E9_a1 <- data[grep("E9_anterior_1", data$Sample),]
data_R1_E9_p4 <- data[grep("E9_posterior_4", data$Sample),]


boxplot(data_R1$Total.Sequences, ylab="PE reads", main="Total number of reads per cell (n=380)")

jpeg("read_couts.jpg")
boxplot(data_R1_E8_1$Total.Sequences, data_R1_E8_2$Total.Sequences, data_R1_E9_a1$Total.Sequences, data_R1_E9_p4$Total.Sequences, 
        names=c("E8.5_#1","E8.5_#2","E9.5_a_#1","E9.5_p_#4"), main="Reads per cell/well", ylab="Number of read pairs")

points(1,mean(data_R1_E8_1$Total.Sequences), pch=18, col="red")
points(2,mean(data_R1_E8_2$Total.Sequences), pch=18, col="red")
points(3,mean(data_R1_E9_a1$Total.Sequences), pch=18, col="red")
points(4,mean(data_R1_E9_p4$Total.Sequences), pch=18, col="red")

dev.off()
