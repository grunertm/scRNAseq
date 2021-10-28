coverage <- read.delim("Isl1_coverage.txt", as.is = T, header = T, sep="\t")


pdf("Isl1_scRNAseq_coverage.pdf")
max_y <- max(coverage$E8.1,coverage$E8.2, coverage$E0_ant, coverage$E9_pos)

plot(coverage$X.Pos, coverage$E8.1, pch=".", main="Isl1 scRNAseq Coverage", xlab="Position (mm10)", ylab="Mean read coverage (n=95 for each set)", ylim=c(0,max_y))
lines(coverage$X.Pos,coverage$E8.1, col=c("black"))
lines(coverage$X.Pos,coverage$E8.2, col=c("darkolivegreen"))
lines(coverage$X.Pos,coverage$E9_pos, col=c("cornflowerblue"))
lines(coverage$X.Pos,coverage$E0_ant, col=c("coral3"))

legend("topright", c("E8_1 (Set A)", "E8_2 (Set C)", "E9 anterior (Set D)", "E9 posterior (Set B)"), lty=1, col = c("black", "darkolivegreen", "cornflowerblue", "coral3"), cex = 0.8)


#Isl1
lines(c(116309689,116298281),c(-1.5,-1.5),col=c("black"))

#rect(116309689,-2,116309395,-1, col=c("gray"))
rect(116309689,-1.75,116309395,-1.25, col=c("gray"))
rect(116309395,-2,116309395-28,-1, col=c("gray"))
text(116309900,-1.5, "5'", cex=0.7)

rect(116308463,-2,116308274,-1, col=c("gray"))
rect(116305477,-2,116305218,-1, col=c("gray"))
rect(116303332,-2,116303046,-1, col=c("gray"))
rect(116301778,-2,116301611,-1, col=c("gray"))

#rect(116299596,-2,116298281,-1, col=c("gray"))
rect(116299596,-1.75,116298281,-1.25, col=c("gray")) 
rect(116299596-117,-2,116299596,-1, col=c("gray"))
text(116298100,-1.5, "3'", cex=0.7)
dev.off()



pdf("Isl1_scRNAseq_coverage_3UTR.pdf")
max_y <- max(coverage$E8.1,coverage$E8.2, coverage$E0_ant, coverage$E9_pos)

plot(coverage$X.Pos, coverage$E8.1, pch=".", main="Isl1 scRNAseq Coverage", xlab="Position (mm10)", ylab="Mean read coverage (n=95 for each set)", xlim=c(116298100,116300000), ylim=c(0,max_y))
lines(coverage$X.Pos,coverage$E8.1, col=c("black"))
lines(coverage$X.Pos,coverage$E8.2, col=c("darkolivegreen"))
lines(coverage$X.Pos,coverage$E9_pos, col=c("cornflowerblue"))
lines(coverage$X.Pos,coverage$E0_ant, col=c("coral3"))

legend("topright", c("E8_1 (Set A)", "E8_2 (Set C)", "E9 anterior (Set D)", "E9 posterior (Set B)"), lty=1, col = c("black", "darkolivegreen", "cornflowerblue", "coral3"), cex = 0.8)


#Isl1
lines(c(116309689,116298281),c(-1.5,-1.5),col=c("black"))

#rect(116309689,-2,116309395,-1, col=c("gray"))
rect(116309689,-1.75,116309395,-1.25, col=c("gray"))
rect(116309395,-2,116309395-28,-1, col=c("gray"))
text(116309900,-1.5, "5'", cex=0.7)

rect(116308463,-2,116308274,-1, col=c("gray"))
rect(116305477,-2,116305218,-1, col=c("gray"))
rect(116303332,-2,116303046,-1, col=c("gray"))
rect(116301778,-2,116301611,-1, col=c("gray"))

#rect(116299596,-2,116298281,-1, col=c("gray"))
rect(116299596,-1.75,116298281,-1.25, col=c("gray")) 
rect(116299596-117,-2,116299596,-1, col=c("gray"))
text(116298100,-1.5, "3'", cex=0.7)

abline(v=116298281, col=c("red"))
abline(v=116298572, col=c("red"))

dev.off()
