genes <- read.delim("genes.TPMs.noZeros.txt", as.is = T, header=T, row.names = 1, sep="\t")
expressed_genes <- genes[rowSums(genes>1)>1,]

#expressed genes (TPM>1)
num_expr_genes <- sapply(expressed_genes, function(x) length(x[x>1]))
mean(num_expr_genes)
median(num_expr_genes)

num_expr_genes_E8_1 <- sapply(expressed_genes[,grep("E8_1", colnames(expressed_genes))], function(x) length(x[x>1]))
num_expr_genes_E8_2 <- sapply(expressed_genes[,grep("E8_2", colnames(expressed_genes))], function(x) length(x[x>1]))
num_expr_genes_E9_anterior_1 <- sapply(expressed_genes[,grep("E9_anterior_1", colnames(expressed_genes))], function(x) length(x[x>1]))
num_expr_genes_E9_anterior_2 <- sapply(expressed_genes[,grep("E9_anterior_2", colnames(expressed_genes))], function(x) length(x[x>1]))
num_expr_genes_E9_posterior_3 <- sapply(expressed_genes[,grep("E9_posterior_3", colnames(expressed_genes))], function(x) length(x[x>1]))
num_expr_genes_E9_posterior_4 <- sapply(expressed_genes[,grep("E9_posterior_4", colnames(expressed_genes))], function(x) length(x[x>1]))

jpeg("gene_couts.jpg", height = 600, width = 800)
boxplot(num_expr_genes_E8_1, num_expr_genes_E8_2, 
        num_expr_genes_E9_anterior_1,
        num_expr_genes_E9_anterior_2,
        num_expr_genes_E9_posterior_3,
        num_expr_genes_E9_posterior_4, ylab="Number of expressed genes (TPM>1)", names=c("E8.5 (#1)","E8.5 (#2)","E9.5 ant (#1)","E9.5 ant (#2)","E9.5 post (#3)","E9.5 post (#4)"))
dev.off()


isoforms <- read.delim("isoforms.TPMs.noZeros.txt", as.is = T, header=T, row.names = 1, sep="\t")
expressed_isoforms <- isoforms[rowSums(isoforms>1)>1,]

#expressed isoforms (TPM>1)
num_expr_isoforms <- sapply(expressed_isoforms, function(x) length(x[x>1]))
mean(num_expr_isoforms)
median(num_expr_isoforms)

num_expr_isoforms_E8_1 <- sapply(expressed_isoforms[,grep("E8_1", colnames(expressed_isoforms))], function(x) length(x[x>1]))
num_expr_isoforms_E8_2 <- sapply(expressed_isoforms[,grep("E8_2", colnames(expressed_isoforms))], function(x) length(x[x>1]))
num_expr_isoforms_E9_anterior_1 <- sapply(expressed_isoforms[,grep("E9_anterior_1", colnames(expressed_isoforms))], function(x) length(x[x>1]))
num_expr_isoforms_E9_anterior_2 <- sapply(expressed_isoforms[,grep("E9_anterior_2", colnames(expressed_isoforms))], function(x) length(x[x>1]))
num_expr_isoforms_E9_posterior_3 <- sapply(expressed_isoforms[,grep("E9_posterior_3", colnames(expressed_isoforms))], function(x) length(x[x>1]))
num_expr_isoforms_E9_posterior_4 <- sapply(expressed_isoforms[,grep("E9_posterior_4", colnames(expressed_isoforms))], function(x) length(x[x>1]))

jpeg("isoform_couts.jpg", height = 600, width = 800)
boxplot(num_expr_isoforms_E8_1, num_expr_isoforms_E8_2, 
        num_expr_isoforms_E9_anterior_1,
        num_expr_isoforms_E9_anterior_2,
        num_expr_isoforms_E9_posterior_3,
        num_expr_isoforms_E9_posterior_4, ylab="Number of expressed isoforms (TPM>1)", names=c("E8.5 (#1)","E8.5 (#2)","E9.5 ant (#1)","E9.5 ant (#2)","E9.5 post (#3)","E9.5 post (#4)"))
dev.off()

