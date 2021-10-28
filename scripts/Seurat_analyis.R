### Analysis using Seurat ###

library(Seurat)
library(cowplot)
library(dplyr)

#load data based on gene counts
data <- read.delim("genes.expected_read_counts.noZeros.gene_names.txt", row.names = 1, header=T, as.is=T, sep = "\t")

#analyisis based on: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html

#Initialize the Seurat object with the raw (non-normalized data)
heart <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
#An object of class Seurat 
#21769 features across 572 samples within 1 assay 
#Active assay: RNA (21769 features)

#QC and selectigng cells for further analysis
heart[["percent.mt"]] <- PercentageFeatureSet(object = heart, pattern = "^mt-")

#samples with multiple cells
samples_multiple_cells <- c(
  "E8_1_Set_A_A1_S1",
  "E8_1_Set_A_A2_S9",
  "E8_1_Set_A_A3_S17",
  "E8_1_Set_A_B1_S2",
  "E8_1_Set_A_B2_S10",
  "E8_1_Set_A_B3_S18",
  "E8_1_Set_A_C1_S3",
  "E8_1_Set_A_C2_S11",
  "E8_1_Set_A_C3_S19",
  "E9.anterior_1_Set_D_A1_S286",
  "E9.anterior_1_Set_D_A2_S294",
  "E9.anterior_1_Set_D_A3_S302",
  "E9.anterior_1_Set_D_B1_S287",
  "E9.anterior_1_Set_D_B2_S295",
  "E9.anterior_1_Set_D_B3_S303",
  "E9.anterior_1_Set_D_C1_S288",
  "E9.anterior_1_Set_D_C2_S296",
  "E9.anterior_1_Set_D_C3_S304",
  "E9.posterior_3_Set_A_A1_S1",
  "E9.posterior_3_Set_A_A2_S9",
  "E9.posterior_3_Set_A_A3_S17",
  "E9.posterior_3_Set_A_B1_S2",
  "E9.posterior_3_Set_A_B2_S10",
  "E9.posterior_3_Set_A_B3_S18",
  "E9.posterior_3_Set_A_C1_S3",
  "E9.posterior_3_Set_A_C2_S11",
  "E9.posterior_3_Set_A_C3_S19")

#Visualize QC metrics as a violin plot
#png("QC_metrics.png", height = 600, width = 800)
VlnPlot(heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#dev.off()

#png("FeatureScatter.png", height = 600, width = 800)
#plot1 <- FeatureScatter(object = heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object = heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1,plot2))
#dev.off()

#sample with multiple cells
#heart2 <- SubsetData(heart, cells = rownames(heart@meta.data) %in% samples_multiple_cells)
#png("QC_metrics.samples_with_multiple_cells.png", height = 600, width = 800)
#VlnPlot(heart2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
#dev.off()

#png("FeatureScatter.samples_with_multiple_cells.png", height = 600, width = 800)
#plot1 <- FeatureScatter(object = heart2, feature1 = "nCount_RNA", feature2 = "percent.mt")
#plot2 <- FeatureScatter(object = heart2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
#CombinePlots(plots = list(plot1,plot2))
#dev.off()

#filter cells based on QC metrics
#21769 features across 465 samples within 1 assay 
heart <- subset(x=heart, subset = nFeature_RNA >3000 & nFeature_RNA <10000 & percent.mt<15)
#21769 features across 462 samples within 1 assay 
heart <- SubsetData(heart, cells = !rownames(heart@meta.data) %in% samples_multiple_cells)

#normalizing the data
heart <- NormalizeData(object = heart)

#identify higly variable features
heart <- FindVariableFeatures(object = heart, selection.method = "vst", nfeatures = 2000)
top20 <- head(x=VariableFeatures(object = heart),20)
#png("Variable_features.png", height = 600, width = 800)
#plot1 <- VariableFeaturePlot(object = heart)
#plot2 <- LabelPoints(plot = plot1, points = top20, repel = TRUE)
#CombinePlots(plots = list(plot1,plot2))
#dev.off()

#scaling the data
all.genes <- rownames(x=heart)
heart <- ScaleData(object = heart, features = all.genes)

#Perform linear dimensional reduction
heart <- RunPCA(object = heart, features = VariableFeatures(object = heart))

#png("PCA.png")
DimPlot(object = heart, reduction = "pca")
#dev.off()

#Heatmap based on 1.PC
#png("PCA_heatmap_dim1.png")
DimHeatmap(object = heart, dims = 1, cells = 500, balanced = TRUE)
#dev.off()

#Heatmap based on 1-15.PCs
png("PCA_heatmap_dim1_15.png")
DimHeatmap(object = heart, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

#determine the dimensionaly of the dataset
heart <- JackStraw(object = heart, num.replicate = 100)
heart <- ScoreJackStraw(object = heart, dims = 1:20)

#png("JackStraPlot.png")
JackStrawPlot(object = heart, dims = 1:20)
#dev.off()

#png("ElbowPlot.png")
ElbowPlot(object = heart)
#dev.off()

#Cluster the cells
heart <- RunUMAP(object = heart, reduction="pca", dims = 1:10)
heart <- FindNeighbors(object = heart, reduction = "pca", dims=1:10)
heart <- FindClusters(heart, resolution = 0.5)

#Non-linear dimensional reduction (UMAP/tSNE)
#UMAP
heart <- RunUMAP(object = heart, dims = 1:10)
#png("UMAP_groups.png",height = 600, width = 800)
DimPlot(object = heart, reduction="umap", group.by = "orig.ident", label = TRUE, pt.size=3)
#dev.off()
#png("UMAP_cluster.png",height = 600, width = 800)
DimPlot(object = heart, reduction="umap", label = TRUE, pt.size=3)
#dev.off()

#tSNE
heart <- RunTSNE(object = heart, dims = 1:10)
#png("plots/tSNE_clusters.png",height = 600, width = 800)
DimPlot(object = heart, reduction="tsne", label = TRUE, pt.size=3)
#dev.off()
#png("plots/tSNE_groups.png",height = 600, width = 800)
DimPlot(object = heart, reduction="tsne",group.by = "orig.ident", pt.size=3, label=TRUE)
#dev.off()

pdf("plots/tSNE_clusters.pdf")
DimPlot(object = heart, reduction="tsne", label = TRUE, pt.size=3)
dev.off()
pdf("plots/tSNE_groups.pdf")
DimPlot(object = heart, reduction="tsne",group.by = "orig.ident", pt.size=3, label=TRUE)
dev.off()

# pdf("PCA_heatmap_dim1.pdf")
# DimHeatmap(object = heart, dims = 1, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim2.pdf")
# DimHeatmap(object = heart, dims = 2, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim3.pdf")
# DimHeatmap(object = heart, dims = 3, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim4.pdf")
# DimHeatmap(object = heart, dims = 4, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim5.pdf")
# DimHeatmap(object = heart, dims = 5, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim6.pdf")
# DimHeatmap(object = heart, dims = 6, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim7.pdf")
# DimHeatmap(object = heart, dims = 7, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim8.pdf")
# DimHeatmap(object = heart, dims = 8, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim9.pdf")
# DimHeatmap(object = heart, dims = 9, cells = 500, balanced = TRUE)
# dev.off()
# pdf("PCA_heatmap_dim10.pdf")
# DimHeatmap(object = heart, dims = 10, cells = 500, balanced = TRUE)
# dev.off()
  
#some test plots for marking specific genes
#DimPlot(object = heart, reduction="tsne", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
#FeatureScatter(heart, feature1 = "Tbx1", feature2 = "Tbx5")
#DotPlot(heart, features = c("Tbx1", "Tbx5"))
#LabelClusters(plot_data, id="ident", color="red") + points(plot_data$data$tSNE_1,plot_data$data$tSNE_2)

#78 cells with Tbx1/Tbx5 coexpression
cells_with_Tbx1_Tbx5 <- WhichCells(object = heart, expression = Tbx1 >0 & Tbx5 >0)

#tSNE plot based on cluster with marked cells (Tbx1>0 and Tbx5>0)
plot_data = DimPlot(object = heart, reduction="tsne")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red", pdata2$colour),]
require(scales)
identities <- levels(heart$seurat_clusters)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_cluster_marked_cells_Tbx1_Tbx5.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Marked cell with Tbx1>0 & Tbx5>0")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$seurat_clusters), pch=19, col=my_colors)
dev.off()

#tSNE plot based on cluster with marked cells (Tbx1>0 and Tbx5>0) - clustered based on origins
plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_origin_marked_cells_Tbx1_Tbx5.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Marked cell with Tbx1>0 & Tbx5>0")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()

#plots to check for batch effects
plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = grep("E9.anterior_1_", rownames(heart@meta.data)))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_origins_marked_cells_batch1_E9_ant_plate1.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Batch 1 E9.5 anterior plate 1")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()

plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = grep("E9.anterior_2_", rownames(heart@meta.data)))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_origins_marked_cells_batch2_E9_ant_plate2.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Batch 2 E9.5 anterior plate 2")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()

plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = grep("E9.posterior_3_", rownames(heart@meta.data)))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_origins_marked_cells_batch2_E9_post_plate3.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Batch 2 E9.5 posterior plate 3")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()

plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = grep("E9.posterior_4_", rownames(heart@meta.data)))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))
pdf("tSNE_origins_marked_cells_batch1_E9_post_plate4.pdf")
plot(pdata$x,pdata$y, col=pdata$colour, pch=19, cex=0.6, xlab="tSNE_1", ylab="tSNE_2", main="Batch 1 E9.5 posterior plate 4")
points(GoI$x,GoI$y)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()


# Cluster biomarkers
cluster0.markers <- FindMarkers(object = heart, ident.1 = 0, min.pct = 0.25)
x<-head(x=cluster0.markers, n=6)
VlnPlot(object = heart, features = rownames(x), slot = "counts", log=TRUE)

cluster1.markers <- FindMarkers(object = heart, ident.1 = 1, min.pct = 0.25)
head(x=cluster1.markers, n=10)
cluster2.markers <- FindMarkers(object = heart, ident.1 = 2, min.pct = 0.25)
head(x=cluster2.markers, n=10)
cluster3.markers <- FindMarkers(object = heart, ident.1 = 3, min.pct = 0.25)
head(x=cluster3.markers, n=10)
cluster4.markers <- FindMarkers(object = heart, ident.1 = 4, min.pct = 0.25)
head(x=cluster4.markers, n=10)
cluster5.markers <- FindMarkers(object = heart, ident.1 = 5, min.pct = 0.25)
x<-head(x=cluster5.markers, n=6)
VlnPlot(object = heart, features = rownames(x), slot = "counts", log=TRUE)


#find markers for every cluster compared to all remaining cells, report only the postive ones
heart.markers <- FindAllMarkers(object = heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)
write.table(heart.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC), "top20_genes_per_cluster.txt", quote=F, sep="\t")

#plotting the top 10 markers for each cluster
top1 <-  heart.markers %>% group_by(cluster) %>% top_n(n=1, wt=avg_logFC)
png("Top1_marker_cluster.png")
FeaturePlot(object = heart, features = c("Hist1h2ap", "Rps2-ps13", "Itm2a", "Nnat", "Actc1", "Sat1"))
dev.off()

png("Top1_marker_cluster.violine.png")
VlnPlot(object = heart, features = c("Hist1h2ap", "Rps2-ps13", "Itm2a", "Nnat", "Actc1", "Sat1"), slot = "counts", log=TRUE)
dev.off()

png("Top5_marker_cluster.heatmap.png")
top5 <-  heart.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
DoHeatmap(object = heart, features = top5$gene) + NoLegend()
dev.off()

DoHeatmap(object = heart, features = top5$gene) + NoLegend()

heart.markers <- FindAllMarkers(object = heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)



### UCSC Cell Browser output ###

heart@misc$markers <- FindAllMarkers(heart)
saveRDS(heart, "heart_small.rds")




########################################################################
################### Genes of interest: Cardiohead ######################
########################################################################


pdf("plots/Cardiohead_genes.UMAP.pdf", width = 12, height = 12)
FeaturePlot(object = heart, features = c("Mef2c","Tbx1", "Fgf10", "Lhx2", "Tbx5", "Hoxb1", "Trp53", "Nkx2-5", "Pitx2"))
dev.off()
pdf("plots/Cardiohead_genes.tSNE.pdf", width = 12, height = 12)
FeaturePlot(object = heart, reduction="tsne", features = c("Mef2c","Tbx1", "Fgf10", "Lhx2", "Tbx5", "Hoxb1", "Trp53", "Nkx2-5", "Pitx2"))
dev.off()


########################################################################
#################### Genes of interest: R. Kelly #######################
########################################################################

pdf("plots/Gene_list_Kelly_I.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Acta1", "Ebf1", "Edn1", "Epha2", "Epha7", "Fgf8", "Meox1", "Mmp9", "Myf5"), reduction="tsne")
FeaturePlot(object = heart, features = c("Myod1", "Pax3", "Prdm1", "Ret", "Sema3c", "Sfrp5", "Six2", "Smoc2", "Tcf21"), reduction="tsne")
FeaturePlot(object = heart, features = c("Tlx1", "Tnc", "Prdm1"), reduction="tsne")
dev.off()

pdf("plots/Gene_list_Kelly_I.UMAP.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Acta1", "Ebf1", "Edn1", "Epha2", "Epha7", "Fgf8", "Meox1", "Mmp9", "Myf5"), reduction="umap")
FeaturePlot(object = heart, features = c("Myod1", "Pax3", "Prdm1", "Ret", "Sema3c", "Sfrp5", "Six2", "Smoc2", "Tcf21"), reduction="umap")
FeaturePlot(object = heart, features = c("Tlx1", "Tnc", "Prdm1"), reduction="umap")
dev.off()

#Lefty2, Pax8, and Wnt5 are not expressed

pdf("plots/Gene_list_Kelly_II.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Aldh1a2","Aldh1a3","Arg1","Cdh3","Cdh4","Clu","Crabp1","Crabp2","Cyp26a1"), reduction="tsne")
FeaturePlot(object = heart, features = c("Cyp26b1","Cyp26b1","Cxcl12","Cxcr4","Dach1","Dlk1","Dsc2","Epha3","Foxc1"), reduction="tsne")
FeaturePlot(object = heart, features = c("Hand2","Hoxa1","Hoxa3","Hoxb4","Gm13715","Irx2","Irx3","Isl1","Kazald1"), reduction="tsne")
FeaturePlot(object = heart, features = c("Lefty1","Lhx6","Mmp2","Moxd1","Msx1","Nkx2-6","Nodal","Nrp2","Osr1"), reduction="tsne")
FeaturePlot(object = heart, features = c("Osr2","Pcdh7","Prdm1","Prrx1","Prrx2","Prtg","Ptpre","Rdh10","Shox2"), reduction="tsne")
FeaturePlot(object = heart, features = c("Snai1","Snai2","Sox9","Sox10","Tbx3","Tbx15","Tcf4","Tdrd12","Tgfb2"), reduction="tsne")
FeaturePlot(object = heart, features = c("Twist1"), reduction="tsne")
dev.off()

pdf("plots/Gene_list_Kelly_II.UMAP.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Aldh1a2","Aldh1a3","Arg1","Cdh3","Cdh4","Clu","Crabp1","Crabp2","Cyp26a1"), reduction="umap")
FeaturePlot(object = heart, features = c("Cyp26b1","Cyp26b1","Cxcl12","Cxcr4","Dach1","Dlk1","Dsc2","Epha3","Foxc1"), reduction="umap")
FeaturePlot(object = heart, features = c("Hand2","Hoxa1","Hoxa3","Hoxb4","Gm13715","Irx2","Irx3","Isl1","Kazald1"), reduction="umap")
FeaturePlot(object = heart, features = c("Lefty1","Lhx6","Mmp2","Moxd1","Msx1","Nkx2-6","Nodal","Nrp2","Osr1"), reduction="umap")
FeaturePlot(object = heart, features = c("Osr2","Pcdh7","Prdm1","Prrx1","Prrx2","Prtg","Ptpre","Rdh10","Shox2"), reduction="umap")
FeaturePlot(object = heart, features = c("Snai1","Snai2","Sox9","Sox10","Tbx3","Tbx15","Tcf4","Tdrd12","Tgfb2"), reduction="umap")
FeaturePlot(object = heart, features = c("Twist1"), reduction="umap")
dev.off()


########################################################################
#################### Genes of interest: TOF genes ######################
########################################################################

#Hcne2, Tceb2 and Wbscr16 are not expressed

pdf("plots/TOF_and_affected_genes.UMAP.pdf", width = 12, height = 12)
FeaturePlot(object = heart, features = c("Barx1", "Bccip", "Dag1", "Edn1", "Fancl", "Fancm", "Fmr1", "Foxk1"), reduction="umap")
FeaturePlot(object = heart, features = c("Myom2","Pex6","Rock1", "Trp53bp2", "Ttn", "Acads", "Mybpc3", "Arvcf", "Sco2"), reduction="umap")
dev.off()

pdf("plots/TOF_and_affected_genes.tSNE.pdf", width = 12, height = 12)
FeaturePlot(object = heart, features = c("Barx1", "Bccip", "Dag1", "Edn1", "Fancl", "Fancm", "Fmr1", "Foxk1"), reduction="tsne")
FeaturePlot(object = heart, features = c("Myom2","Pex6","Rock1", "Trp53bp2", "Ttn", "Acads", "Mybpc3", "Arvcf", "Sco2"), reduction="tsne")
dev.off()


########################################################################
#################### Analysis based on Transcripts #####################
########################################################################

transcripts <- read.delim("isoforms.expected_read_counts.noZeros.txt", row.names = 1, header=T, as.is=T, sep = "\t")
heart2 <- CreateSeuratObject(transcripts, min.cells = 3, min.features = 200)
#76367 features across 572 samples within 1 assay 

heart2[["percent.mt"]] <- PercentageFeatureSet(object = heart2, pattern = "^mt-")

png("QC_metrics.transcripts.png", height = 600, width = 800)
VlnPlot(heart2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()

png("FeatureScatter.transcripts.png", height = 600, width = 800)
plot1 <- FeatureScatter(object = heart2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = heart2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
dev.off()

#samples with multiple cells
heart3 <- SubsetData(heart2, cells = rownames(heart2@meta.data) %in% samples_multiple_cells)

png("QC_metrics.samples_with_multiple_cells.transcripts.png", height = 600, width = 800)
VlnPlot(heart3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol=3)
dev.off()

png("FeatureScatter.samples_with_multiple_cells.transcript.png", height = 600, width = 800)
plot1 <- FeatureScatter(object = heart3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(object = heart3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1,plot2))
dev.off()

#filter cells
#76367 features across 469 samples within 1 assay 
heart2 <- subset(x=heart2, subset = nFeature_RNA >5000 & nFeature_RNA <20000 & percent.mt<15)
#76367 features across 460 samples within 1 assay
heart2 <- SubsetData(heart2, cells = !rownames(heart2@meta.data) %in% samples_multiple_cells)

#normalize data
heart2 <- NormalizeData(object = heart2)

#identify higly variable features
heart2 <- FindVariableFeatures(object = heart2, selection.method = "vst", nfeatures = 2000)
top20 <- head(x=VariableFeatures(object = heart2),20)
png("Variable_features.transcripts.png", height = 600, width = 800)
plot1 <- VariableFeaturePlot(object = heart2)
LabelPoints(plot = plot1, points = top20, repel = TRUE)
dev.off()

#find markers for every cluster compared to all remaining cells, report only the postive ones
heart.markers2 <- FindAllMarkers(object = heart2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart.markers2 %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)
write.table(heart.markers2 %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC), "top20_transcripts_per_cluster.txt", quote=F, sep="\t")

#scaling the data
all.transcripts <- rownames(x=heart2)
heart2 <- ScaleData(object = heart2, features = all.transcripts)

#Perform linear dimensional reduction
heart2 <- RunPCA(object = heart2, features = VariableFeatures(object = heart2))

png("PCA.transcripts.png")
DimPlot(object = heart2, reduction = "pca")
dev.off()


#determine the dimensionaly of the dataset
heart2 <- JackStraw(object = heart2, num.replicate = 100)
heart2 <- ScoreJackStraw(object = heart2, dims = 1:20)

png("JackStraPlot.transcripts.png")
JackStrawPlot(object = heart2, dims = 1:20)
dev.off()

ElbowPlot(object = heart2)

#Cluster the cells
heart2 <- RunUMAP(object = heart2, reduction="pca", dims = 1:10)
heart2 <- FindNeighbors(object = heart2, reduction = "pca", dims=1:10)
heart2 <- FindClusters(heart2, resolution = 0.5)

#UMAP
heart2 <- RunUMAP(object = heart2, dims = 1:10)
png("UMAP_groups.transcripts.png",height = 600, width = 800)
DimPlot(object = heart2, reduction="umap", group.by = "orig.ident")
dev.off()
png("UMAP_cluster.transcripts.png",height = 600, width = 800)
DimPlot(object = heart2, reduction="umap")
dev.off()

#tSNE
heart2 <- RunTSNE(object = heart2, dims = 1:10)
png("tSNE_clusters.transcripts.png",height = 600, width = 800)
DimPlot(object = heart2, reduction="tsne")
dev.off()
png("tSNE_groups.transcripts.png",height = 600, width = 800)
DimPlot(object = heart2, reduction="tsne", group.by = "orig.ident")
dev.off()


################################################################################################################################################
########################################## Pseudotime/trajectory inference using Monocle #######################################################
################################################################################################################################################

library(monocle)

  #function to import Seurat object 
  newimport <- function(otherCDS, import_all = TRUE) {
    if(class(otherCDS)[1] == 'Seurat') {
      requireNamespace("Seurat")
      data <- otherCDS@assays$RNA@counts
      
      if(class(data) == "data.frame") {
        data <- as(as.matrix(data), "sparseMatrix")
      }
      
      pd <- tryCatch( {
        pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
        pd
      }, 
      #warning = function(w) { },
      error = function(e) { 
        pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
        pd <- new("AnnotatedDataFrame", data = pData)
        
        message("This Seurat object doesn't provide any meta data");
        pd
      })
      
      # remove filtered cells from Seurat
      if(length(setdiff(colnames(data), rownames(pd))) > 0) {
        data <- data[, rownames(pd)]  
      }
      
      fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
      fd <- new("AnnotatedDataFrame", data = fData)
      lowerDetectionLimit <- 0
      
      if(all(data == floor(data))) {
        expressionFamily <- negbinomial.size()
      } else if(any(data < 0)){
        expressionFamily <- uninormal()
      } else {
        expressionFamily <- tobit()
      }
      
      valid_data <- data[, row.names(pd)]
      
      monocle_cds <- newCellDataSet(data,
                                    phenoData = pd, 
                                    featureData = fd,
                                    lowerDetectionLimit=lowerDetectionLimit,
                                    expressionFamily=expressionFamily)
      
      if(import_all) {
        if("Monocle" %in% names(otherCDS@misc)) {
          otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
          otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
          
          monocle_cds <- otherCDS@misc$Monocle
          mist_list <- otherCDS
          
        } else {
          # mist_list <- list(ident = ident) 
          mist_list <- otherCDS
        }
      } else {
        mist_list <- list()
      }
      
      if(1==1) {
        var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
        
      }
      monocle_cds@auxClusteringData$seurat <- mist_list
      
    } else if (class(otherCDS)[1] == 'SCESet') {
      requireNamespace("scater")
      
      message('Converting the exprs data in log scale back to original scale ...')    
      data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
      
      fd <- otherCDS@featureData
      pd <- otherCDS@phenoData
      experimentData = otherCDS@experimentData
      if("is.expr" %in% slotNames(otherCDS))
        lowerDetectionLimit <- otherCDS@is.expr
      else 
        lowerDetectionLimit <- 1
      
      if(all(data == floor(data))) {
        expressionFamily <- negbinomial.size()
      } else if(any(data < 0)){
        expressionFamily <- uninormal()
      } else {
        expressionFamily <- tobit()
      }
      
      if(import_all) {
        # mist_list <- list(iotherCDS@sc3,
        #                   otherCDS@reducedDimension)
        mist_list <- otherCDS 
        
      } else {
        mist_list <- list()
      }
      
      monocle_cds <- newCellDataSet(data,
                                    phenoData = pd, 
                                    featureData = fd,
                                    lowerDetectionLimit=lowerDetectionLimit,
                                    expressionFamily=expressionFamily)
      # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
      # monocle_cds@auxOrderingData$scran <- mist_list
      
      monocle_cds@auxOrderingData$scran <- mist_list
      
    } else {
      stop('the object type you want to export to is not supported yet')
    }
    
    return(monocle_cds)
  }

#create monocle object from Seurat object
monocle <- newimport(heart)
  
View(pData(monocle))
View(fData(monocle))

#trajectory step1: choose genes that defince a cell's process
var_genes <- heart[["RNA"]]@var.features
ordering_genes <- var_genes  

monocle <- setOrderingFilter(monocle, ordering_genes)  
print(dim(exprs(monocle)))  


#trajectory step2: reduce data dimensionality
#reduce the space down to one with two dimensions
monocle <- reduceDimension(monocle, reduction_method = "DDRTree", max_components = 2)

#trajectory step3: order cells along the trajectory
monocle <- orderCells(monocle)

pdf("monocle_clusters.Tbx1.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers="Tbx1")
dev.off()
pdf("monocle_clusters.Tbx5.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers="Tbx5")
dev.off()

pdf("monocle_origins.Tbx1.pdf")
plot_cell_trajectory(monocle, color_by = 'orig.ident', markers="Tbx1")
dev.off()
pdf("monocle_origins.Tbx5.pdf")
plot_cell_trajectory(monocle, color_by = 'orig.ident', markers="Tbx5")
dev.off()

#visualize the trajectory in the reduced dimensional space
#the tracjectory has a tree-like structure
#pdf("monocle_clusters.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters')
#dev.off()
#pdf("monocle_origins.pdf")
plot_cell_trajectory(monocle, color_by = 'orig.ident')
#dev.off()

#Facet the trajetories
#pdf("monocle_clusters_facet.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters') +
  facet_wrap(~seurat_clusters, nrow = 1)
#dev.off()

#pdf("monocle_origins_facet.pdf")
plot_cell_trajectory(monocle, color_by = 'orig.ident') +
  facet_wrap(~orig.ident, nrow = 1)
#dev.off()

#plot kinetic trends for a couple of marker genes such as "Mef2c","Tbx1", "Fgf10", "Lhx2", "Tbx5", "Hoxb1", "Trp53", "Nkx2-5", "Pitx2"
my_genes <- row.names(subset(fData(monocle),gene_short_name %in% c("Tbx1","Tbx5", "Nkx2-5")))
monocle_subset <- monocle[my_genes,]
plot_genes_in_pseudotime(monocle_subset, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% c("Mef2c","Fgf10", "Lhx2")))
monocle_subset <- monocle[my_genes,]
plot_genes_in_pseudotime(monocle_subset, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% c("Hoxb1","Trp53", "Pitx2")))
monocle_subset <- monocle[my_genes,]
plot_genes_in_pseudotime(monocle_subset, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% c("Tbx1","Tbx5")))
#pdf("monocle_Tbx1_Tbx5.pdf")
plot_genes_branched_pseudotime(monocle[my_genes,], branch_point = 1, color_by = "orig.ident", ncol = 1)
#dev.off()
#pdf("monocle_Tbx1_Tbx5.clusters.pdf")
plot_genes_branched_pseudotime(monocle[my_genes,], branch_point = 1, color_by = "seurat_clusters", ncol = 1)
#dev.off()


#Alternative choices for ordering genes
#http://cole-trapnell-lab.github.io/monocle-release/docs/#alternative-choices-for-ordering-genes

#select superset of feature genes as genes expressed in at least 5% of all the cells
monocle <- detectGenes(monocle, min_expr = 0.1)
fData(monocle)$use_for_ordering <- fData(monocle)$num_cells_expressed > 0.05 * ncol(monocle)

#expressed genes
monocle_expressed_genes <- row.names(subset(fData(monocle), use_for_ordering==TRUE))
#monocle_expressed_genes <- row.names(subset(fData(monocle), num_cells_expressed >= 10))

#perform differential gene expression test as a way to extract the genes that distinguish the clusters
clustering_DEG_genes <- differentialGeneTest(monocle[monocle_expressed_genes,], fullModelFormulaStr = '~seurat_clusters', cores = 8)

# We will then select the top 1000 significant genes as the ordering genes. 
monocle_ordering_genes <- row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:200]

monocle <- setOrderingFilter(monocle, ordering_genes = monocle_ordering_genes)
monocle <- reduceDimension(monocle, method = 'DDRTree', max_components = 2)
monocle <- orderCells(monocle)

# plot_cell_trajectory(monocle, color_by = "seurat_clusters") + theme(legend.position = "right")
# plot_cell_trajectory(monocle, color_by = 'seurat_clusters') + facet_wrap(~seurat_clusters, nrow = 1)
# plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers=c("Tbx1","Tbx5"), markers_linear = TRUE)
# plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers="Tbx5", markers_linear = TRUE, show_state_number = TRUE)

#constructing single cell trajectories based on the top 200 DEGs and mark intermediate states
#During development, in response to stimuli, and througout life, cells transition from one functional "state" to another. Cells in different states express different sets of genes.
pdf("plots/trajectories_DEG200.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = TRUE, cell_size = 3)
plot_cell_trajectory(monocle, color_by = 'orig.ident', show_state_number = TRUE, cell_size = 3)
dev.off()

pdf("plots/trajectories_DEG200_noStates.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters')
plot_cell_trajectory(monocle, color_by = 'orig.ident')
dev.off()

state_cells <- monocle@phenoData[['State']] 
sample_cells <- monocle@phenoData[['seurat_clusters']]
df <- data.frame(state_cells, sample_cells)
table(df)
#               sample_cells
# state_cells   0   1   2   3   4   5
#           1   0  61   0   9   0   4
#           2   0  17   1   2   5   0
#           3   1  16   2   4   0   2
#           4 112  26  38  46   0   7
#           5   3   0   2   0   0   0
#           6  17   3  45  14  11   3
#           7   0   9   0   0   2   0
sample_cells <- monocle@phenoData[['orig.ident']]
df <- data.frame(state_cells, sample_cells)
table(df)
#              sample_cells
# state_cells  E8 E9.anterior E9.posterior
#           1  74           0            0
#           2  21           4            0
#           3  21           2            2
#           4  24          70          135
#           5   0           1            4
#           6   5          56           32
#           7  10           0            1



#Finding Genes that Change as a Function of Pseudotime
#Once we have a trajectory, we can use differentialGeneTest() to find genes that 
#have an expression pattern that varies according to pseudotime.
#https://davetang.org/muse/2017/10/01/getting-started-monocle/

#In Monocle, a single cell trajectory is the inferred developmental timeline of single cells. 
#The “unit” used for the trajectory is pseudotime; a cell at the beginning of the trajectory, i.e. starting state, 
#will have a lower pseudotime than cells at the end of the trajectory, i.e. end state. Specifically, it is:
#"the distance between a cell and the start of the trajectory, measured along the shortest path. The trajectory’s total length 
#is defined in terms of the total amount of transcriptional change that a cell undergoes as it moves from the starting state to the end state." (Monocle tutorial)

#The trajectory analysis consists of three stages:
# 1. Choose genes that define progress
# 2. Reduce the dimensionality of the data
# 3. Order cells in pseudotime
# For the trajectory analysis one can use only a subset of all cells (e.g., cluster 1, 4, and 6) and genes that are expressed in at least 10 cells.

#Finding Genes that Change as a Function of Pseudotime
#Once we have a trajectory, we can use differentialGeneTest() to find genes that have an expression pattern that varies according to pseudotime.
head(pData(monocle))
monocle_de <- differentialGeneTest(monocle, fullModelFormulaStr = "~Pseudotime", cores = 3)
monocle_de %>% arrange(qval) %>% head(n=10)
# status family pval qval gene_short_name use_for_ordering num_cells_expressed
# 1      OK  tobit    0    0            Agrp            FALSE                   4
# 2      OK  tobit    0    0        B4galnt2            FALSE                   4
# 3      OK  tobit    0    0           Hnf4a            FALSE                   3
# 4      OK  tobit    0    0          Slc6a4            FALSE                   4
# 5      OK  tobit    0    0          Slc5a7            FALSE                   6
# 6      OK  tobit    0    0         Ifi202b            FALSE                   4
# 7      OK  tobit    0    0          Hsd3b6            FALSE                   6
# 8      OK  tobit    0    0           Plbd1            FALSE                   5
# 9      OK  tobit    0    0         Creb3l3            FALSE                   4
# 10     OK  tobit    0    0          Prss53            FALSE                   4

# save the top 10 genes
monocle_de %>% arrange(qval) %>% head(n=10) %>% select(gene_short_name) -> monocle_de_genes
monocle_de_genes <- monocle_de_genes$gene_short_name
plot_genes_in_pseudotime(monocle[monocle_de_genes,])

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% c("Agrp","B4galnt2","Hnf4a","Slc6a4","Slc5a7","Ifi202b","Hsd3b6","Plbd1","Creb3l3","Prss53")))
monocle_subset <- monocle[my_genes,]
pdf("top10_DEGs_in_pseudotime.pdf")
plot_genes_in_pseudotime(monocle_subset, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')
dev.off()

#Analyzing Branches in Single-Cell Trajectories
#In our trajectory, we have two branches, which represents cells that have alternative gene expression patterns. 
#These represent cells that have supposedly gone through different developmental paths. Monocle provides functions that 
#allows you to identify the genes that differ at a particular branch point. 

#The BEAM() function takes a CellDataSet that has been ordered with orderCells() and a branch point in the trajectory. 
#A table of genes is returned with significance values that indicate whether genes have expression patterns that are branch dependent.

BEAM_res <- BEAM(monocle, branch_point = 1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]

head(BEAM_res)

table(BEAM_res$qval < 1e-4)
# FALSE  TRUE 
# 15196   999 
table(BEAM_res$qval ==0)
# FALSE  TRUE 
# 15900   295 

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% row.names(subset(BEAM_res, qval ==0 ))))
monocle_subset <- monocle[my_genes,]
my_branched_heatmap <- plot_genes_branched_heatmap(monocle_subset, branch_point = 1, num_clusters = 6, cores = 3, use_gene_short_name = TRUE, show_rownames = TRUE, return_heatmap = TRUE)

#We can return genes that belong to specific clusters that were identified by BEAM().
head(my_branched_heatmap$annotation_row)
# 
# Cluster
# Pih1d2         1
# Serpinf1       2
# Txnrd3         2
# Hpn            2
# Nom1           3
# Fbxo7          2

dim(my_branched_heatmap$annotation_row)
# [1] 253   1

table(my_branched_heatmap$annotation_row$Cluster)
# 1  2  3  4  5  6 
# 25 46 81 32 44 25 

my_row <- my_branched_heatmap$annotation_row
my_row <- data.frame(cluster = my_row$Cluster,
                     gene = row.names(my_row),
                     stringsAsFactors = FALSE)

head(my_row[my_row$cluster == 4,'gene'])
# [1] "Klf5"   "Hspb7"  "Cd34"   "Map4k2" "Tbx20"  "Fam81a"

my_gene <- row.names(subset(fData(monocle_subset), gene_short_name %in% head(my_row[my_row$cluster == 4,'gene'])))

# plot genes that are expressed in a branch dependent manner
#The trend lines show the expression pattern of genes along the paths formed by branch point 1.
pdf("branch_point_1_cluster4.pdf")
plot_genes_branched_pseudotime(monocle_subset[my_gene,],
                               branch_point = 1,
                               ncol = 1)
dev.off()




############################################# Tbx1 & Tbx5 ############################################# 

pdf("plots/trajectories_Tbx1_Tbx5.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers = c("Tbx1","Tbx5"), use_color_gradient = FALSE, markers_linear = TRUE, show_state_number = FALSE, cell_size = 3)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', markers = c("Tbx1","Tbx5"), use_color_gradient = FALSE, markers_linear = TRUE, show_state_number = TRUE, cell_size = 3)
plot_cell_trajectory(monocle, color_by = 'orig.ident', markers = c("Tbx1","Tbx5"), use_color_gradient = FALSE, markers_linear = TRUE, show_state_number = FALSE, cell_size = 3)
plot_cell_trajectory(monocle, color_by = 'orig.ident', markers = c("Tbx1","Tbx5"), use_color_gradient = FALSE, markers_linear = TRUE, show_state_number = TRUE, cell_size = 3)
dev.off()

#462 cells: 78 with Tbx1>0 and Tbx5>0 (22xE8; 11x E9 ant; 46x E9 post)
#22 cells with Tbx1>0.2 and Tbx5>0.2 (4xE8;18x E9 post)
#cells_with_Tbx1_Tbx5[grep("E8",cells_with_Tbx1_Tbx5)]
cells_with_Tbx1_Tbx5 <- WhichCells(object = heart, expression = Tbx1>0 & Tbx5>0)
cells_no_Tbx1_Tbx5 <- WhichCells(object = heart, expression = !(Tbx1>0 & Tbx5>0))

pdf("plots/Tbx1_Tbx5_coexpression.pdf")
DimPlot(object = heart, reduction="tsne", cells = c(cells_with_Tbx1_Tbx5))
DimPlot(object = heart, reduction="tsne", cells = c(cells_with_Tbx1_Tbx5), group.by = "orig.ident")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5)
dev.off()

#n=462 cells
#48 Tbx1==0 & Tbx5==0
#307 Tbx1>0 & Tbx5==0: 236 Tbx1<=1; 71 Tbx1>1
#29 Tbx1==0 & Tbx5>0: 16 Tbx5<0.5; 13 Tbx5>0 
#78 Tbx1>0 & Tbx5>0
cells_with_Tbx5_only <- WhichCells(object = heart, expression = Tbx1==0 & Tbx5>0) #29 cells
cells_with_strong_Tbx1_only <- WhichCells(object = heart, expression = Tbx1>1 & Tbx5==0) #71 cells
cells_with_strong_Tbx5_only <- WhichCells(object = heart, expression = Tbx1==0 & Tbx5>0.5) #13 cells
cells_with_moderate_Tbx1_only <- WhichCells(object = heart, expression = Tbx1>0 & Tbx1<=1 & Tbx5==0) #236 cells
cells_with_moderate_Tbx5_only <- WhichCells(object = heart, expression = Tbx1==0 & Tbx5>0 & Tbx5<=0.5 ) #16 cells
cells_without_Tbx1_Tbx5 <- WhichCells(object = heart, expression = Tbx1==0 & Tbx5==0) #48 cells
cells_with_Tbx1_Tbx5 <- WhichCells(object = heart, expression = Tbx1>0 & Tbx5>0) #78 cells

category <- colnames(heart)
category[category %in% cells_without_Tbx1_Tbx5] <- "no_Tbx1_Tbx5"
category[category %in% cells_with_moderate_Tbx1_only] <- "moderate_Tbx1_only"
category[category %in% cells_with_strong_Tbx1_only] <- "strong_Tbx1_only"
category[category %in% cells_with_moderate_Tbx5_only] <- "moderate_Tbx5_only"
category[category %in% cells_with_strong_Tbx5_only] <- "strong_Tbx5_only"
category[category %in% cells_with_Tbx1_Tbx5] <- "Tbx1_Tbx5"

heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "Tbx1_Tbx5_expression"
)

#Finds markers (differentially expressed genes) for identity classes
cat1.markers <- FindMarkers(object = heart, ident.1 = "no_Tbx1_Tbx5", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat1.genes <- row.names(head(x=cat1.markers, n=5))
cat2.markers <- FindMarkers(object = heart, ident.1 = "moderate_Tbx1_only", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat2.genes <- row.names(head(x=cat2.markers, n=1))
cat3.markers <- FindMarkers(object = heart, ident.1 = "strong_Tbx1_only", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat3.genes <- row.names(head(x=cat3.markers, n=5))
cat4.markers <- FindMarkers(object = heart, ident.1 = "moderate_Tbx5_only", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat4.genes <- row.names(head(x=cat4.markers, n=5))
cat5.markers <- FindMarkers(object = heart, ident.1 = "strong_Tbx5_only", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat5.genes <- row.names(head(x=cat5.markers, n=5))
cat6.markers <- FindMarkers(object = heart, ident.1 = "Tbx1_Tbx5", min.pct = 0.25, group.by = "Tbx1_Tbx5_expression")
cat6.genes <- row.names(head(x=cat6.markers, n=5))

cat.genes <- unique(c(cat1.genes, cat2.genes, cat3.genes, cat4.genes, cat5.genes, cat6.genes))
pdf("plots/heatmap_Tbx1_Tbx5_expression.pdf")
DoHeatmap(object = heart, features = cat.genes, group.by = "Tbx1_Tbx5_expression", size = 3)
dev.off()

#cat.heart <- SetIdent(object = heart, value = "Tbx1_Tbx5_expression")
#cat.markers <- FindAllMarkers(object = cat.heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#x <- cat.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
#write.table(x, "cat.markers.txt", sep="\t", quote = FALSE)
#DoHeatmap(object = cat.heart, features = x$gene, size=3)
#View(x)




pdf("plots/Tbx1_Tbx5_expression.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_only, cols.highlight = "blue") + ggtitle(label = "Expression of Tbx1>0 & Tbx5==0") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5, cols.highlight = "red") + ggtitle(label = "Expression of Tbx1>0 & Tbx5>0") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx5_only, cols.highlight = "darkgreen") + ggtitle(label = "Expression of Tbx1==0 & Tbx5>0") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_without_Tbx1_Tbx5, cols.highlight = "black") + ggtitle(label = "Expression of Tbx1==0 & Tbx5==0") + NoLegend()
dev.off()

png("plots/Tbx1_only.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_only) + ggtitle(label = "Expression of Tbx1>0 & Tbx5==0") + NoLegend() + scale_color_manual(values = c("transparent", "blue") )
dev.off()
png("plots/Tbx5_only.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx5_only) + ggtitle(label = "Expression of Tbx1==0 & Tbx5>0") + NoLegend() + scale_color_manual(values = c("transparent", "aquamarine3") )
dev.off()
png("plots/Tbx1_Tbx5.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5) + ggtitle(label = "Expression of Tbx1>0 & Tbx5>0") + NoLegend() + scale_color_manual(values = c("transparent", "red") )
dev.off()
png("plots/no_Tbx1_Tbx5.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_without_Tbx1_Tbx5) + ggtitle(label = "Expression of Tbx1==0 & Tbx5==0") + NoLegend() + scale_color_manual(values = c("transparent", "black") )
dev.off()

png("plots/strong_Tbx1_only.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_strong_Tbx1_only) + ggtitle(label = "Expression of Tbx1>1 & Tbx5==0") + NoLegend() + scale_color_manual(values = c("transparent", "cyan") )
dev.off()
png("plots/strong_Tbx5_only.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_strong_Tbx5_only) + ggtitle(label = "Expression of Tbx1==0 & Tbx5>0.5") + NoLegend() + scale_color_manual(values = c("transparent", "chartreuse") )
dev.off()

#Tbx1, Tbx5 and Mef2c-cre co-expressing cells
Cre <- read.delim("../Cre_YFP_analysis/Cre_alignment_results.txt", sep="\t", as.is = T, header=T, row.names = 1)
Cell_Cre_10reads <- WhichCells(heart,cells =rownames(Cre[Cre$reads>10,])) #n=461/462
Cell_Cre_100reads <- WhichCells(heart,cells =rownames(Cre[Cre$reads>100,])) #n=387/462
Cell_Cre_1000reads <- WhichCells(heart,cells =rownames(Cre[Cre$reads>1000,])) #n=114/462

cells_with_Tbx1_Tbx5_Cre1000 <- WhichCells(object = heart, expression = Tbx1>0 & Tbx5>0, cells=Cell_Cre_1000reads) #8/78 cells
cells_with_Tbx1_Tbx5_Cre100 <- WhichCells(object = heart, expression = Tbx1>0 & Tbx5>0, cells=Cell_Cre_100reads) #59/78

pdf("plots/Tbx1_Tbx5_Cre_expression.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5, cols.highlight = "black", pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Tbx1>0 & Tbx5>0") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5_Cre100, cols.highlight = "black", pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Tbx1>0 & Tbx5>0 & Cre>100reads") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5_Cre1000, cols.highlight = "black", pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Tbx1>0 & Tbx5>0 & Cre>1000reads") + NoLegend()
dev.off()


########################################## Hoxb1 & Tbx5 ########################################## 

#n=173 cells with Hoxb1 expression (out of 462 cells)
cells_with_Hoxb1 <- WhichCells(object = heart, expression = Hoxb1>0)
Hoxb1_subset <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Hoxb1)
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Hoxb1, cols.highlight = "blue") + ggtitle(label = "Expression of Hoxb1>0") + NoLegend()
pdf("plots/Hoxb1.pdf")
FeaturePlot(object = heart, features = c("Hoxb1"), reduction="tsne")
dev.off()

x = FetchData(object = Hoxb1_subset, vars = c('seurat_clusters', 'ident', "Tbx5", "Hoxb1") )
median(x$Hoxb1) #0.349527
mean(x$Hoxb1) #0.4154911
median(x$Tbx5) #0
mean(x$Tbx5) #0.0674181

#n=42
cells_with_Hoxb1_Tbx5 <- WhichCells(object = heart, expression = Hoxb1>0 & Tbx5>0)
Hoxb1_subset1 <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Hoxb1_noTbx5)
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Hoxb1_Tbx5, cols.highlight = "red") + ggtitle(label = "Expression of Hoxb1>0") + NoLegend()
table(Hoxb1_subset1@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 71            6           54 
table(Hoxb1_subset1@meta.data[['seurat_clusters']])
#0  1  2  3  4  5 
#45 66  9  6  4  1

#n=131
cells_with_Hoxb1_noTbx5 <- WhichCells(object = heart, expression = Hoxb1>0 & Tbx5==0)
Hoxb1_subset2 <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Hoxb1_Tbx5)
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Hoxb1_noTbx5, cols.highlight = "blue") + ggtitle(label = "Expression of Hoxb1>0") + NoLegend()
table(Hoxb1_subset2@meta.data[['orig.ident']])
#E8  E9.anterior E9.posterior 
#20            1           21
table(Hoxb1_subset2@meta.data[['seurat_clusters']])
#0  1  2  3  4  5 
#13 19  9  1  0  0

png("plots/Hoxb1_Tbx5.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Hoxb1_Tbx5) + ggtitle(label = "Expression of Hoxb1>0 & Tbx5>0") + NoLegend() + scale_color_manual(values = c("transparent", "blue") )
dev.off()
png("plots/Hoxb1_noTbx5.png", bg = "transparent")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Hoxb1_noTbx5) + ggtitle(label = "Expression of Hoxb1>0 & Tbx5==0") + NoLegend() + scale_color_manual(values = c("transparent", "aquamarine3") )
dev.off()
png("plots/noHoxb1_noTbx5.png", bg = "transparent")
cells_with_noHoxb1_noTbx5 <- WhichCells(object = heart, expression = Hoxb1==0 & Tbx5==0)
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_noHoxb1_noTbx5) + ggtitle(label = "Expression of Hoxb1>0 & Tbx5==0") + NoLegend() + scale_color_manual(values = c("transparent", "grey") )
dev.off()

#DEGs in Hoxb1>0 & Tbx5>0 versus Hoxb1>0 & Tbx5==0
#Identifies differentially expressed genes between two groups of cells using a Wilcoxon Rank Sum test (default)

#Stores unnormalized data such as raw counts or TPMs
Hoxb1.raw.data <- as.matrix(GetAssayData(heart, slot = "counts")[,  WhichCells(object = heart, expression = Hoxb1>0)])
DEGs_Hoxb1 = FindMarkers(Hoxb1.raw.data, cells.1 = cells_with_Hoxb1_Tbx5, cells.2 = cells_with_Hoxb1_noTbx5, verbose = T)
# p_val    avg_logFC pct.1 pct.2    p_val_adj
# Tbx5          2.360752e-38  112.2623304 1.000 0.000 5.139121e-34
# Wnt2          5.244056e-14  341.1375277 0.810 0.221 1.141579e-09
# Gm43050       1.008187e-08    7.2889061 0.238 0.000 2.194722e-04
# Tnfrsf19      1.017585e-08   51.1433316 0.786 0.328 2.215180e-04
# Gm266         5.786829e-08 -118.8624723 0.714 0.939 1.259735e-03
# Cxcl13        5.999780e-08   53.2623304 0.214 0.000 1.306092e-03
# Rspo2         7.230156e-07          Inf 0.381 0.076 1.573933e-02
# Tbx1          1.019376e-06 -202.5556195 0.762 0.962 2.219080e-02

pdf("plots/DEGS_of_Hoxb1_Tbx5_pos_cells.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Tbx5","Wnt2","Gm43050","Tnfrsf19","Gm266","Cxcl13","Rspo2","Tbx1"), reduction="tsne")
dev.off()

#Normalized data matrix
Hoxb1.data <- as.matrix(GetAssayData(heart, slot = "data")[,  WhichCells(object = heart, expression = Hoxb1>0)])
DEGs_Hoxb1_data = FindMarkers(Hoxb1.data, cells.1 = cells_with_Hoxb1_Tbx5, cells.2 = cells_with_Hoxb1_noTbx5, verbose = T)
# p_val    avg_logFC pct.1 pct.2    p_val_adj
# Tbx5                       2.484167e-38  0.3158806 1.000 0.000 5.407783e-34
# Wnt2                       6.469280e-14  0.6236596 0.810 0.221 1.408298e-09
# Gm266                      2.727806e-09 -0.4287265 0.714 0.939 5.938162e-05
# Tbx1                       2.131141e-08 -0.4210730 0.762 0.962 4.639280e-04
# Rspo2                      6.545764e-07  0.4964418 0.381 0.076 1.424947e-02

#Scaled data matrix
Hoxb1.scaled_data <- as.matrix(GetAssayData(heart, slot = "scale.data")[,  WhichCells(object = heart, expression = Hoxb1>0)]) 
DEGs_Hoxb1_scaled_data = FindMarkers(Hoxb1.scaled_data, cells.1 = cells_with_Hoxb1_Tbx5, cells.2 = cells_with_Hoxb1_noTbx5, verbose = T)
# p_val    avg_logFC pct.1 pct.2    p_val_adj
# Tbx5          2.484167e-38  2.2996010 0.667 0.000 5.407783e-34
# Wnt2          6.469280e-14  1.7508273 0.619 0.099 1.408298e-09
# Gm266         2.727806e-09 -1.2207284 0.190 0.611 5.938162e-05
# Tnfrsf19      6.686295e-09  1.1828436 0.571 0.137 1.455540e-04
# Gm43050       1.008911e-08  3.3096320 0.238 0.000 2.196298e-04
# Tbx1          2.131141e-08 -0.9973145 0.167 0.626 4.639280e-04
# Cxcl13        6.000491e-08  4.5379832 0.190 0.000 1.306247e-03
# Rspo2         6.545764e-07  4.0885682 0.286 0.038 1.424947e-02

DEGs_Hoxb1_negbinom = FindMarkers(Hoxb1.raw.data, cells.1 = cells_with_Hoxb1_Tbx5, cells.2 = cells_with_Hoxb1_noTbx5, test.use = "negbinom", verbose = T, pseudocount.use = 0.1)
write.table(DEGs_Hoxb1_negbinom,"DEGs_Hoxb1_negbinom.txt", sep = "\t", quote = FALSE)

#Heatmap 
cells_with_Hoxb1_Tbx5 <- WhichCells(object = heart, expression = Hoxb1>0 & Tbx5>0) #42
cells_with_Hoxb1_noTbx5 <- WhichCells(object = heart, expression = Hoxb1>0 & Tbx5==0) #131
cells_with_noHoxb1_noTbx5 <- WhichCells(object = heart, expression = Hoxb1==0 & Tbx5==0) #224
cells_with_noHoxb1_Tbx5 <- WhichCells(object = heart, expression = Hoxb1==0 & Tbx5>0) #65

category <- colnames(heart)
category[category %in% cells_with_Hoxb1_Tbx5] <- "Hoxb1_Tbx5"
category[category %in% cells_with_Hoxb1_noTbx5] <- "Hoxb1_but_no_Tbx5"
category[category %in% cells_with_noHoxb1_noTbx5] <- "no_Hoxb1_no_Tbx5"
category[category %in% cells_with_noHoxb1_Tbx5] <- "no_Hoxb1_but_Tbx5"

heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "Hoxb1_Tbx5_expression"
)

#Finds markers (differentially expressed genes) for identity classes
cat1.markers <- FindMarkers(object = heart, ident.1 = "Hoxb1_Tbx5", min.pct = 0.25, group.by = "Hoxb1_Tbx5_expression")
cat1.genes <- row.names(head(x=cat1.markers, n=4))
cat2.markers <- FindMarkers(object = heart, ident.1 = "Hoxb1_but_no_Tbx5", min.pct = 0.25, group.by = "Hoxb1_Tbx5_expression")
cat2.genes <- row.names(head(x=cat2.markers, n=5))
cat3.markers <- FindMarkers(object = heart, ident.1 = "no_Hoxb1_no_Tbx5", min.pct = 0.25, group.by = "Hoxb1_Tbx5_expression")
cat3.genes <- row.names(head(x=cat3.markers, n=5))
cat4.markers <- FindMarkers(object = heart, ident.1 = "no_Hoxb1_but_Tbx5", min.pct = 0.25, group.by = "Hoxb1_Tbx5_expression")
cat4.genes <- row.names(head(x=cat4.markers, n=5))

cat.genes <- unique(c(cat1.genes, cat2.genes, cat3.genes, cat4.genes))
pdf("plots/heatmap_Hoxb1_Tbx5_expression.pdf")
DoHeatmap(object = heart, features = cat.genes, group.by = "Hoxb1_Tbx5_expression", size = 3)
dev.off()



#################################### RA target and pathay genes ###########################################

#RA signaling 
#Identify RA target and pathway genes (GO term Retinoic acid signaling pathway)
#Mine the transcriptome of Raldh2 (Aldh1a2) expressing cells, which incude the RA catabolic enzyme Cyp26 that antagnonizes RA induction in the anterior SHF.
#Cyp26 (Cyp26a1,Cyp26c1,Cyp26b1) and Tbx5 exclusive expression


#There 30 genes associated with the Go term 'retinoic acid signaling pathway'
#http://www.informatics.jax.org/vocab/gene_ontology/GO:0048384

#Akr1c18 and Ptf1a are not expressed and thus not part of the following list
GO_0048384_genes <- c("Actn4","Aldh1a2","Aldh1a3","Asxl1","Calr","Cnot1","Crabp2","Crkl","Ctbp2","Cyp26a1","Cyp26b1","Dhrs3","Esrrg","Ezh2","Nr1h2","Nr2c1","Pml",
                      "Rara","Rarb","Rarg","Rxra","Rxrb","Rxrg","Snw1","Tbx1","Tgif1","Trim16","Zfp536")

#Akr1c18 and Ptf1a are not expressed
pdf("plots/Gene_list_RAgenes.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Actn4","Aldh1a2","Aldh1a3","Asxl1","Calr","Cnot1","Crabp2","Crkl","Ctbp2"), reduction="tsne")
FeaturePlot(object = heart, features = c("Cyp26a1","Cyp26b1","Dhrs3","Esrrg","Ezh2","Nr1h2","Nr2c1","Pml","Rara"), reduction="tsne")
FeaturePlot(object = heart, features = c("Rarb","Rarg","Rxra","Rxrb","Rxrg","Snw1","Tbx1", "Tgif1", "Trim16"), reduction="tsne")
FeaturePlot(object = heart, features = c("Zfp536"), reduction="tsne")
dev.off()

pdf("plots/heatmap_RA_genes.pdf")
DoHeatmap(object = heart, features = GO_0048384_genes, group.by = "seurat_clusters", size = 3)
DoHeatmap(object = heart, features = GO_0048384_genes, group.by = "orig.ident", size = 3)
dev.off()

#Raldh2 (Aldh1a2)
pdf("plots/Raldh2_featureplot.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Aldh1a2"), reduction="tsne", pt.size = 3)
dev.off()
cells_with_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2>0) #n=186
Aldh1a2_subset <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Aldh1a2)
pdf("plots/Raldh2_dimplot.tSNE.pdf", height = 12, width = 12)
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_with_Aldh1a2, cols.highlight = "blue", pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Aldh1a2(Raldh2)>0") + NoLegend()
dev.off()
table(Aldh1a2_subset@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 77           12           97
table(Aldh1a2_subset@meta.data[['seurat_clusters']])
# 0  1  2  3  4  5 
# 78 70 19  6  8  5

cells_with_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2>1) #n=70
Aldh1a2_subset2 <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Aldh1a2)
pdf("plots/Raldh2_highExp_dimplot.tSNE.pdf", height = 12, width = 12)
DimPlot(object = heart, reduction="tsne", cells.highlight =  WhichCells(object = heart, expression = Aldh1a2>1), cols.highlight = "blue", pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Aldh1a2(Raldh2)>0") + NoLegend()
dev.off()
table(Aldh1a2_subset2@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 24           0           46
table(Aldh1a2_subset2@meta.data[['seurat_clusters']])
# 0  1  2  3  4  5 
# 39 25 6  0  0  0


#Stores unnormalized data such as raw counts or TPMs
cells_with_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2>0) #n=186
cells_with_strong_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2>1) #n=70
cells_with_moderate_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2>0 & Aldh1a2<=1) #n=116
cells_with_no_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2==0) #n=276
cells_with_less_Aldh1a2 <- WhichCells(object = heart, expression = Aldh1a2<=1) #n=392

#DEGs 70 vs. 392 cells
Aldh1a2.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))
DEGs_Hoxb1 = FindMarkers(Aldh1a2.raw.data, cells.1 = cells_with_strong_Aldh1a2, cells.2 = cells_with_less_Aldh1a2, verbose = T)

head(DEGs_Hoxb1, 20)
pdf("plots/heatmap_top20_DEGs_of_strong_Aldh1a2_expressing_cells.pdf")
DoHeatmap(object = heart, features = row.names(head(DEGs_Hoxb1, 20)), group.by = "seurat_clusters", size = 3)
DoHeatmap(object = heart, features = row.names(head(DEGs_Hoxb1, 20)), group.by = "orig.ident", size = 3)
dev.off()

pdf("plots/Cyp26_featureplot.tSNE.pdf", height = 12, width = 12)
FeaturePlot(object = heart, features = c("Cyp26a1","Cyp26b1","Cyp26c1"), reduction="tsne")
dev.off()

#138 cells that express Cyp26
cells_with_Cyp26 <- WhichCells(object = heart, expression = Cyp26a1>0 | Cyp26b1 >0 | Cyp26c1>0)
#26/138 cells also express Tbx5 (19/91 with Cyp26a1>0; 10/41 with Cyp26b1>0; 2/41 with Cyp26c1>0)
cells_with_Cyp26_Tbx5 <- WhichCells(object = heart, expression = (Cyp26a1>0 | Cyp26b1>0 | Cyp26c1>0) & Tbx5>0)
pdf("plots/Cyp26_Tbx5_dimplot.tSNE.pdf", height = 12, width = 12)
DimPlot(object = heart, reduction="tsne", cells.highlight =  cells_with_Cyp26_Tbx5, pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Expression of Cyp26>0 & Tbx5>0") + NoLegend() + scale_color_manual(values = c("lightgrey", "black") ) 
dev.off()
Cyp26_Tbx5_subset <- SubsetData(heart, cells = rownames(heart@meta.data) %in% cells_with_Cyp26_Tbx5)
table(Cyp26_Tbx5_subset@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 11            2           13 
table(Cyp26_Tbx5_subset@meta.data[['seurat_clusters']])
# 0  1  2  3  4  5 
# 4 11  8  1  2  0

DoHeatmap(object = Cyp26_Tbx5_subset, group.by = "orig.ident", size = 3)
DoHeatmap(object = Cyp26_Tbx5_subset, group.by = "seurat_clusters", size = 3)

cells_with_Cyp26_Tbx5 <- WhichCells(object = heart, expression = (Cyp26a1>0 | Cyp26b1>0 | Cyp26c1>0) & Tbx5>0)
cells_with_no_Cyp26_Tbx5 <- WhichCells(object = heart, expression = !((Cyp26a1>0 | Cyp26b1>0 | Cyp26c1>0) & Tbx5>0))

#15 DEGs
Cyp26.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))
DEGs_Cyp26 = FindMarkers(Cyp26.raw.data, cells.1 = cells_with_Cyp26_Tbx5, cells.2 = cells_with_no_Cyp26_Tbx5, verbose = T)
View(DEGs_Cyp26)

pdf("plots/heatmap_DEGs_Cyp26.pdf")
DoHeatmap(object = Cyp26_Tbx5_subset, features = row.names(head(DEGs_Cyp26,15)), group.by = "seurat_clusters", size = 3)
dev.off()

#Cyp26.markers <- FindAllMarkers(object = Cyp26_Tbx5_subset, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#x <- Cyp26.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_logFC)
#DoHeatmap(object = Cyp26_Tbx5_subset, features = x$gene, size=3)
#View(x)


my_genes <- row.names(subset(fData(monocle),gene_short_name %in% c("Cyp26a1","Cyp26b1","Cyp26c1","Tbx5")))
monocle_subset <- monocle[my_genes,]
plot_genes_in_pseudotime(monocle, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')

my_genes <- row.names(subset(fData(monocle), gene_short_name %in% c("Mef2c","Fgf10", "Lhx2")))
monocle_subset <- monocle[my_genes,]
plot_genes_in_pseudotime(monocle_subset, color_by = 'seurat_clusters')
plot_genes_in_pseudotime(monocle_subset, color_by = 'orig.ident')

#Raldh2 -> Aldh1a2
#(Aldehyde dehydrogenase 1 family, member A2, also known as ALDH1A2 or retinaldehyde dehydrogenase 2 (RALDH2))

#https://www.sciencedirect.com/science/article/pii/S0012160614004709
#known RA targets
# Cyp26a1 (Loudig et al., 2000)
# Hoxa1 (Marshall et al., 1996)
# Cdx1 (Houle et al., 2000)
# Arg1 (Chang et al., 2013)
# Dhrs3 (Feng et al., 2010)
# Lhx1 (Hunter and Rhodes, 2005, Cartry et al., 2006)
# Ptprz (Paschaki et al., 2013)
# Gcnf1 (Heinzer et al., 1998, Barreto et al., 2003)
# Mrg1 (Oulad-Abdelghani et al., 1997)
# Nrip1 (Kerley et al., 2001). 
# Wnt5a (Kumar and Duester, 2010)
# T (Iulianella et al., 1999


#################################### New cell populations in the development of trapezius muscle ###########################################

#derives from the most posterior SHF domain 
#Cell must be posterior expressing Tbx1, Isl1, Ebf but not Tbx5 and not Mef2c-Cre (probably)
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4321263/
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5468682/
#Ebf: Ebf1, Ebf2, Ebf3 --> Collier, Olf And EBF Transcription Factors

#Tbx1, Tbx5 and Mef2c-cre co-expressing cells
Cre <- read.delim("../Cre_YFP_analysis/Cre_alignment_results.txt", sep="\t", as.is = T, header=T, row.names = 1)
Cell_Cre_less_100reads <- WhichCells(heart,cells =rownames(Cre[Cre$reads<100,])) #n=75/462
DimPlot(object = heart, reduction="tsne", cells.highlight =  Cell_Cre_less_100reads) + ggtitle(label = "Cre neg") + NoLegend() + scale_color_manual(values = c("lightgrey", "black") ) 

pdf("plots/Trapezius_muccle_cells_dimplot.pdf")

#n=193
trapezius_muscle_cells2 <- WhichCells(object = heart, expression = Tbx1>0 & Isl1 >0 & (Ebf1>0 | Ebf2>0 | Ebf3>0) & Tbx5==0)
DimPlot(object = heart, reduction="tsne", cells.highlight =  trapezius_muscle_cells2, pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Trapezius muscle cells") + NoLegend() + scale_color_manual(values = c("lightgrey", "black") ) 
TM_subset <- SubsetData(heart, cells = rownames(heart@meta.data) %in% trapezius_muscle_cells2)
table(TM_subset@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 68           64           61 
table(TM_subset@meta.data[['seurat_clusters']])
# 0  1  2  3  4  5 
# 50 56 28 47  5  7

#n=14
trapezius_muscle_cells2 <- WhichCells(object = heart, expression = (Tbx1>0 & Isl1 >0 & (Ebf1>0 | Ebf2>0 | Ebf3>0) & Tbx5==0), cells = Cell_Cre_less_100reads)
DimPlot(object = heart, reduction="tsne", cells.highlight =  trapezius_muscle_cells2, pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Trapezius muscle cells (Cre neg)") + NoLegend() + scale_color_manual(values = c("lightgrey", "black") ) 
TM_subset <- SubsetData(heart, cells = rownames(heart@meta.data) %in% trapezius_muscle_cells2)
table(TM_subset@meta.data[['orig.ident']])
# E8  E9.anterior E9.posterior 
# 11            2            1
table(TM_subset@meta.data[['seurat_clusters']])
# 0  1  2  3  4  5 
# 2 11  0  0  0  1 

dev.off()

#n=29
pdf("plots/Trapezius_muscle_cells_dimplot_cre_pos.pdf")
trapezius_muscle_cells2 <- WhichCells(object = heart, expression = Tbx1>0 & Isl1 >0 & Ebf1>0 & Ebf2>0 & Ebf3>0 & Tbx5==0)
DimPlot(object = heart, reduction="tsne", cells.highlight =  trapezius_muscle_cells2, pt.size = 3, sizes.highlight = 3) + ggtitle(label = "Trapezius muscle cells (Cre pos)") + NoLegend() + scale_color_manual(values = c("lightgrey", "black") ) 
dev.off()

#Heatmap of DEGs for 14 cells
trapezius_muscle_cells2 <- WhichCells(object = heart, expression = (Tbx1>0 & Isl1 >0 & Ebf1>0 & Ebf2>0 & Ebf3>0 & Tbx5==0))
no_trapezius_muscle_cells2 <- WhichCells(object = heart, expression = !(Tbx1>0 & Isl1 >0 & Ebf1>0 & Ebf2>0 & Ebf3>0 & Tbx5==0))
trapezius.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))
DEGs_trapezius = FindMarkers(trapezius.raw.data, cells.1 = trapezius_muscle_cells2, cells.2 = no_trapezius_muscle_cells2, verbose = T)
View(DEGs_trapezius)

pdf("plots/heatmap_DEGs_trapezius.pdf")
DoHeatmap(object = heart, features = row.names(head(DEGs_trapezius,47)), group.by = "seurat_clusters", size = 3)
DoHeatmap(object = heart, features = row.names(head(DEGs_trapezius,47)), group.by = "orig.ident", size = 3)
dev.off()



################################################################################################################################################
####################################################### Trajectory analysis using Slingshot ####################################################
################################################################################################################################################

library("slingshot")

## get tSNE embedings
cca.data <- FetchData(heart, vars = c("tSNE_1","tSNE_2"))

## get cluster ids
cnames <- colnames(heart@meta.data)
clus <- heart@meta.data[,grep("res",cnames)]

## run slingshot
sce <- slingshot(cca.data,clusterLabels=clus)

## line
#lin1 <- getLineages(cca.data, clus, start.clus= '0', end.clus = '5')
lin1 <- getLineages(cca.data, clus, start.clus= '0')
lin1 <- getCurves(lin1)


#plot 1
plot_data = DimPlot(object = heart, reduction="tsne")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$seurat_clusters)
my_colors <- hue_pal()(length(identities))

pdf("slingshot.clusters.pdf")
plot(cca.data, col = pdata$colour, pch=19, cex=0.6, asp = 1)
lines(lin1, lwd = 3, show.constraints = TRUE)
legend("bottomleft", legend=levels(heart$seurat_clusters), pch=19, col=my_colors)
dev.off()

pdf("slingshot.clusters.curved4.pdf")
plot(cca.data, col = pdata$colour, pch=19, cex=0.6)
lines(lin1, lwd = 3, type = 'c')
legend("bottomleft", legend=levels(heart$seurat_clusters), pch=19, col=my_colors)
dev.off()

#plot 2
plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))

pdf("slingshot.origins.curved.pdf")
plot(cca.data, col = pdata$colour, pch=19, cex=0.6, asp = 1)
lines(lin1, lwd = 3, show.constraints = TRUE)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()

lin0 <- getLineages(cca.data, clus, start.clus= '2')
lin0 <- getCurves(lin0)

plot_data = DimPlot(object = heart, reduction="tsne", group.by = "orig.ident")
pbuild <- ggplot2::ggplot_build(plot_data)
pdata <- pbuild$data[[1]]
plot_data2 <- DimPlot(object = heart, reduction="tsne", group.by = "orig.ident", label = TRUE, cells.highlight = WhichCells(heart, expression = Tbx1 >0 & Tbx5 >0))
pbuild2 <- ggplot2::ggplot_build(plot_data2)
pdata2 <- pbuild2$data[[1]]
GoI <- pdata2[grep("red",pdata2$colour),]
require(scales)
identities <- levels(heart$orig.ident)
my_colors <- hue_pal()(length(identities))

pdf("slingshot.origins.started_at_cl2.curved.pdf")
plot(cca.data, col = pdata$colour, pch=19, cex=0.6, asp = 1)
lines(lin0, lwd = 3, show.constraints = TRUE)
legend("bottomleft", legend=levels(heart$orig.ident), pch=19, col=my_colors)
dev.off()


#################################################################################################################################
####################################################### Analysis @ Marseille ####################################################
#################################################################################################################################

##############################################
### sub-clustering for cells in cluster 0 ####
##############################################
heart_small.c0 <- subset(heart, idents = 0)
rownames(heart_small.c0@meta.data)
cluster0 <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
cluster0 <- SubsetData(cluster0, cells = rownames(cluster0@meta.data) %in% rownames(heart_small.c0@meta.data))
cluster0 <- NormalizeData(object = cluster0)
cluster0 <- FindVariableFeatures(object = cluster0, selection.method = "vst", nfeatures = 500)
all.genes0 <- rownames(x=cluster0)
cluster0 <- ScaleData(object = cluster0, features = all.genes0)
cluster0 <- RunPCA(object = cluster0, features = VariableFeatures(object = cluster0))
cluster0 <- JackStraw(object = cluster0, num.replicate = 100)
cluster0 <- ScoreJackStraw(object = cluster0, dims = 1:20)
JackStrawPlot(object = cluster0, dims = 1:20)
cluster0 <- RunUMAP(object = cluster0, reduction="pca", dims = 1:4)
cluster0 <- FindNeighbors(object = cluster0, reduction = "pca", dims=1:4)
cluster0 <- FindClusters(cluster0, resolution = 0.5)
cluster0 <- RunUMAP(object = cluster0, dims = 1:4)
cluster0 <- RunTSNE(object = cluster0, dims = 1:4)

pdf("Subclusering_of_cluster0.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight =  rownames(heart_small.c0@meta.data), pt.size = 2, sizes.highlight = 2) + 
  ggtitle(label = "Cluster 0") + NoLegend() + scale_color_manual(values = c("lightgrey", "red") ) 
DimPlot(object = cluster0, reduction="tsne", label = TRUE, pt.size=2)
FeaturePlot(object = cluster0, features = c("Tbx5", "Tbx1", "Scx", "Gdnf"), reduction = "tsne", pt.size = 2)
FeaturePlot(object = cluster0, features = c("Aldh1a2", "Wnt4", "Shox2", "Klf16"), reduction = "tsne", pt.size = 2)
DimHeatmap(object = cluster0, dims = 1, cells = 133, balanced = TRUE)
dev.off()

#cell browser output
cluster0@misc$markers <- FindAllMarkers(cluster0)
saveRDS(cluster0, "cluster0.rds")




pdf("top9_genes_subclustering.pdf")
# Cluster biomarkers
cluster0.markers <- FindMarkers(object = cluster0, ident.1 = 0, min.pct = 0.25)
x<-head(x=cluster0.markers, n=9)
VlnPlot(object = cluster0, features = rownames(x), slot = "counts", log=TRUE)

cluster1.markers <- FindMarkers(object = cluster0, ident.1 = 1, min.pct = 0.25)
x<-head(x=cluster1.markers, n=9)
VlnPlot(object = cluster0, features = rownames(x), slot = "counts", log=TRUE)

cluster2.markers <- FindMarkers(object = cluster0, ident.1 = 2, min.pct = 0.25)
x<-head(x=cluster2.markers, n=9)
VlnPlot(object = cluster0, features = rownames(x), slot = "counts", log=TRUE)

cluster3.markers <- FindMarkers(object = cluster0, ident.1 = 3, min.pct = 0.25)
x<-head(x=cluster3.markers, n=9)
VlnPlot(object = cluster0, features = rownames(x), slot = "counts", log=TRUE)
dev.off()

#find markers for every cluster compared to all remaining cells, report only the postive ones
heart.markers <- FindAllMarkers(object = cluster0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)
write.table(heart.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC), "top20_genes_per_cluster.txt", quote=F, sep="\t")

png("Top10_marker_cluster.heatmap.png")
top5 <-  heart.markers %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)
DoHeatmap(object = cluster0, features = top5$gene) + NoLegend()
dev.off()













##############################################
### sub-clustering for cells in cluster 3 ####
##############################################
heart_small.c3 <- subset(heart, idents = 3)
rownames(heart_small.c3@meta.data)
cluster3 <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
cluster3 <- SubsetData(cluster3, cells = rownames(cluster3@meta.data) %in% rownames(heart_small.c3@meta.data))
cluster3 <- NormalizeData(object = cluster3)
cluster3 <- FindVariableFeatures(object = cluster3, selection.method = "vst", nfeatures = 500)
all.genes3 <- rownames(x=cluster3)
cluster3 <- ScaleData(object = cluster3, features = all.genes3)
cluster3 <- RunPCA(object = cluster3, features = VariableFeatures(object = cluster3))
cluster3 <- JackStraw(object = cluster3, num.replicate = 100)
cluster3 <- ScoreJackStraw(object = cluster3, dims = 1:20)
JackStrawPlot(object = cluster3, dims = 1:20)
cluster3 <- RunUMAP(object = cluster3, reduction="pca", dims = 1:2)
cluster3 <- FindNeighbors(object = cluster3, reduction = "pca", dims=1:2)
cluster3 <- FindClusters(cluster3, resolution = 0.5)
cluster3 <- RunUMAP(object = cluster3, dims = 1:2)
#cluster3 <- RunTSNE(object = cluster3, dims = 1:2)

pdf("Subclusering_of_cluster3.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight =  rownames(heart_small.c3@meta.data), pt.size = 3, sizes.highlight = 3) + 
  ggtitle(label = "Cluster 3") + NoLegend() + scale_color_manual(values = c("lightgrey", "red") ) 
DimPlot(object = cluster3, reduction="pca", label = TRUE, pt.size=2)
DimHeatmap(object = cluster3, dims = 1, cells = 133, balanced = TRUE)
dev.off()

#cell browser output
cluster3@misc$markers <- FindAllMarkers(cluster3)
saveRDS(cluster3, "cluster3.rds")


####################################################################
### Scleraxis analysis: Scx+ cells in cluster 0 versus cluster 2 ###
####################################################################
Scx_pos_cells_cluster0_handselected <- c("E9.posterior_3_Set_A_A8_S57", "E9.posterior_3_Set_A_A11_S81", "E9.posterior_3_Set_A_H8_S64","E9.posterior_4_Set_B_B1_S97",
                                         "E9.posterior_3_Set_A_H10_S80", "E9.posterior_4_Set_B_C9_S161", "E9.posterior_3_Set_A_G4_S31", "E9.posterior_4_Set_B_B9_S160",  "E9.posterior_4_Set_B_H9_S166", "E9.anterior_2_Set_B_E1_S101")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster0_handselected, pt.size = 1, sizes.highlight = 1)

Scx_neg_cells_cluster0_handselected <- c("E9.posterior_4_Set_B_F1_S101","E9.posterior_4_Set_B_H4_S127","E9.posterior_4_Set_B_A2_S104","E9.posterior_3_Set_A_A6_S41",
                                         "E9.posterior_3_Set_A_H11_S88","E9.posterior_3_Set_A_A7_S49","E9.posterior_4_Set_B_A12_S183","E9.posterior_4_Set_B_G1_S102",
                                         "E9.posterior_3_Set_A_H9_S72","E9.posterior_4_Set_B_B8_S153","E9.posterior_4_Set_B_H3_S119", "E9.posterior_3_Set_A_C7_S51","E9.anterior_2_Set_B_H2_S112")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster0_handselected, pt.size = 1, sizes.highlight = 1)

Scx_pos_cells_cluster2_handselected <- c("E9.anterior_2_Set_B_A12_S185","E9.posterior_3_Set_A_H5_S40","E9.anterior_2_Set_B_H6_S144","E9.posterior_3_Set_A_F2_S14",
                                         "E9.posterior_3_Set_A_F11_S86", "E9.posterior_4_Set_B_C5_S130", "E9.anterior_2_Set_B_B7_S146","E9.anterior_1_Set_D_D2_S297","E9.anterior_2_Set_B_C4_S123","E9.anterior_2_Set_B_F9_S166")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster2_handselected, pt.size = 1, sizes.highlight = 1)

Scx_neg_cells_cluster2_handselected <- c("E9.anterior_1_Set_D_D11_S368","E9.anterior_1_Set_D_B11_S366","E9.posterior_3_Set_A_E9_S69","E9.posterior_3_Set_A_E3_S21",
                                         "E9.posterior_3_Set_A_G6_S47", "E9.posterior_3_Set_A_B9_S66","E9.posterior_3_Set_A_D10_S76","E9.posterior_4_Set_B_C7_S146",
                                         "E9.posterior_3_Set_A_D11_S84","E9.anterior_1_Set_D_H9_S356","E9.anterior_2_Set_B_H9_S168","E9.anterior_2_Set_B_B2_S106",
                                         "E9.anterior_2_Set_B_F5_S134","E9.anterior_1_Set_D_D12_S376","E9.anterior_1_Set_D_E12_S377")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster2_handselected, pt.size = 1, sizes.highlight = 1)

pdf("Scleraxis_selected_cells.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster0_handselected, pt.size = 1, sizes.highlight = 1) + ggtitle(label = "Hand-selected cells+ (cluster0)") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster0_handselected, pt.size = 1, sizes.highlight = 1) + ggtitle(label = "Hand-selected cells- (cluster0)") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster2_handselected, pt.size = 1, sizes.highlight = 1) + ggtitle(label = "Hand-selected cells+ (cluster2)") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster2_handselected, pt.size = 1, sizes.highlight = 1) + ggtitle(label = "Hand-selected cells- (cluster2)") + NoLegend()
dev.off()
  
heart.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))
DEGs_Scx_pos_neg_cluster0 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster0_handselected, cells.2 = Scx_neg_cells_cluster0_handselected, verbose = T)
View(DEGs_Scx_pos_neg_cluster0)
DEGs_Scx_pos_neg_cluster2 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster2_handselected, cells.2 = Scx_neg_cells_cluster2_handselected, verbose = T)
View(DEGs_Scx_pos_neg_cluster2)  
DEGs_Scx_pos_cluster0_vs_2 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster0_handselected, cells.2 = Scx_pos_cells_cluster2_handselected, verbose = T)
View(DEGs_Scx_pos_cluster0_vs_2)

Scx_pos_cells_cluster0 <- WhichCells(object = heart, expression = (Scx>0), idents = 0)
Scx_pos_cells_cluster2 <- WhichCells(object = heart, expression = (Scx>0), idents = 2)
DEGs_Scx_cluster0_vs_2 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster0, cells.2 = Scx_pos_cells_cluster2, verbose = T)
View(DEGs_Scx_cluster0_vs_2)  

Scx_pos_cells_cluster0 <- WhichCells(object = heart, expression = (Scx>0), idents = 0)
Scx_neg_cells_cluster0 <- WhichCells(object = heart, expression = (Scx==0), idents = 0)
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster0, pt.size = 1, sizes.highlight = 1)
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster0, pt.size = 1, sizes.highlight = 1)
DEGs_Scx_cluster0 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster0, cells.2 = Scx_neg_cells_cluster0, verbose = T)
View(DEGs_Scx_cluster0)  

Scx_pos_cells_cluster2 <- WhichCells(object = heart, expression = (Scx>0), idents = 2)
Scx_neg_cells_cluster2 <- WhichCells(object = heart, expression = (Scx==0), idents = 2)
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_pos_cells_cluster2, pt.size = 1, sizes.highlight = 1)
DimPlot(object = heart, reduction="tsne", cells.highlight =  Scx_neg_cells_cluster2, pt.size = 1, sizes.highlight = 1)
DEGs_Scx_cluster2 = FindMarkers(heart.raw.data, cells.1 = Scx_pos_cells_cluster2, cells.2 = Scx_neg_cells_cluster2, verbose = T)
View(DEGs_Scx_cluster2)  


#########################################
### Tbx1 to Tbx5 boundary/transititon ###
#########################################
# Selction based on Gdnf, Klf16, Cldn11, Foxf2 and Aldh1a2 (n=10 cells with certain overlaps)

Genes_pos <- c("E9.posterior_4_Set_B_C6_S138", "E9.posterior_4_Set_B_H1_S103", "E9.posterior_4_Set_B_A11_S175", 
         "E9.posterior_3_Set_A_A9_S65", "E9.posterior_4_Set_B_A4_S120", "E9.posterior_4_Set_B_H2_S111", 
         "E9.posterior_3_Set_A_H4_S32","E9.posterior_4_Set_B_C1_S98", 
         "E9.posterior_3_Set_A_H3_S24", "E9.posterior_3_Set_A_C11_S83", 
         "E9.posterior_3_Set_A_C6_S43", "E9.posterior_3_Set_A_F10_S78", "E9.posterior_3_Set_A_D4_S28")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Genes_pos, pt.size = 1, sizes.highlight = 1, cols.highlight = c("red"))

Genes_neg <- c("E9.anterior_2_Set_B_F10_S174", "E9.posterior_3_Set_A_H2_S16", "E9.posterior_3_Set_A_A10_S73", "E9.posterior_4_Set_B_A1_S96",
            "E9.posterior_4_Set_B_D1_S99", "E9.posterior_4_Set_B_B6_S137", "E9.posterior_3_Set_A_G12_S95",
             "E9.posterior_3_Set_A_D2_S12", "E9.posterior_3_Set_A_E11_S85", "E9.posterior_3_Set_A_A4_S25", "E9.posterior_3_Set_A_E2_S13",
             "E9.posterior_3_Set_A_G5_S39", "E9.posterior_3_Set_A_B12_S90", "E9.posterior_3_Set_A_C4_S27",
             "E9.posterior_4_Set_B_A7_S144", "E9.posterior_3_Set_A_E1_S5",
             "E9.posterior_4_Set_B_C12_S185","E9.posterior_4_Set_B_A10_S167","E9.posterior_3_Set_A_H1_S8","E9.posterior_4_Set_B_H12_S190",
             "E9.posterior_4_Set_B_H5_S135", "E9.posterior_3_Set_A_B11_S82", "E9.posterior_4_Set_B_C11_S177", "E9.posterior_3_Set_A_F9_S70",
             "E9.posterior_3_Set_A_F3_S22", "E9.posterior_3_Set_A_G11_S87", "E9.posterior_3_Set_A_F7_S54", "E9.posterior_3_Set_A_E6_S45",
             "E9.posterior_4_Set_B_D11_S178", "E9.posterior_4_Set_B_H10_S174",
             "E9.posterior_3_Set_A_C12_S91", "E9.posterior_3_Set_A_G3_S23", "E9.posterior_4_Set_B_E12_S187")
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg, pt.size = 1, sizes.highlight = 1, cols.highlight = "black")

Genes_neg_top <- c("E9.posterior_3_Set_A_A4_S25", "E9.posterior_3_Set_A_E2_S13","E9.posterior_3_Set_A_G5_S39", "E9.posterior_3_Set_A_B12_S90", "E9.posterior_3_Set_A_C4_S27",
                   "E9.posterior_4_Set_B_A7_S144", "E9.posterior_3_Set_A_E1_S5","E9.posterior_4_Set_B_C12_S185","E9.posterior_4_Set_B_A10_S167","E9.posterior_3_Set_A_H1_S8",
                   "E9.posterior_4_Set_B_D11_S178", "E9.posterior_4_Set_B_H10_S174","E9.posterior_4_Set_B_E12_S187","E9.posterior_4_Set_B_H12_S190","E9.posterior_3_Set_A_E6_S45",
                   "E9.posterior_3_Set_A_F7_S54")
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg_top, pt.size = 1, sizes.highlight = 1, cols.highlight = "black")

Genes_neg_bottom <- c("E9.anterior_2_Set_B_F10_S174", "E9.posterior_3_Set_A_H2_S16", "E9.posterior_3_Set_A_A10_S73", "E9.posterior_4_Set_B_A1_S96",
                      "E9.posterior_4_Set_B_D1_S99", "E9.posterior_4_Set_B_B6_S137", "E9.posterior_3_Set_A_G12_S95","E9.posterior_3_Set_A_D2_S12",
                      "E9.posterior_4_Set_B_H5_S135", "E9.posterior_3_Set_A_B11_S82", "E9.posterior_4_Set_B_C11_S177", "E9.posterior_3_Set_A_F9_S70",
                      "E9.posterior_3_Set_A_C12_S91", "E9.posterior_3_Set_A_G3_S23", "E9.posterior_3_Set_A_E11_S85", "E9.posterior_3_Set_A_F3_S22",
                      "E9.posterior_3_Set_A_G11_S87")
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg_bottom, pt.size = 1, sizes.highlight = 1, cols.highlight = "black")

pdf("cells_selected_for_Tbx1_Tbx5_transition.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight =  Genes_pos, pt.size = 1, sizes.highlight = 1, cols.highlight = c("red")) + ggtitle(label = "Cells in boundary") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg, pt.size = 1, sizes.highlight = 1, cols.highlight = "black") + ggtitle(label = "Cells in transition (top/bottom)") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg_top, pt.size = 1, sizes.highlight = 1, cols.highlight = "black") + ggtitle(label = "Cells in transition (top)") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = Genes_neg_bottom, pt.size = 1, sizes.highlight = 1, cols.highlight = "black") + ggtitle(label = "Cells in transition (bottom)") + NoLegend()
dev.off()

heart.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))

#DEGs in cell+ vs. cell- (top & bottom)
##DEGs: Cldn11, Gm5127
DEGs_GoI_pos_neg_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg, verbose = T)
View(DEGs_GoI_pos_neg_wilcox)  
write.table(DEGs_GoI_pos_neg_wilcox, "genes_of_interest.wilcox_analysis.txt", sep="\t", quote = F)
DEGs_GoI_pos_neg_ROC = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg, verbose = T, test.use = "roc")
View(DEGs_GoI_pos_neg_ROC)  
write.table(DEGs_GoI_pos_neg, "genes_of_interest.ROC_analysis.txt", sep="\t", quote = F)

#cell+ vs. cell- (top)
#DEGs: -
DEGs_GoI_pos_neg_top_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg_top, verbose = T)
View(DEGs_GoI_pos_neg_top_wilcox)  

#cell+ vs. cell- (bottom)
#DEGs: Foxf1
DEGs_GoI_pos_neg_bottom_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg_bottom, verbose = T)
View(DEGs_GoI_pos_neg_bottom_wilcox)  

#assign cells to categories (i.e., cells in the boundary and cells in transition (top/bottom))
category <- colnames(heart)
category[category %in% Genes_pos] <- "cells+"
category[category %in% Genes_neg_top] <- "cell- top"
category[category %in% Genes_neg_bottom] <- "cell- bottom"

heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "Tbx1_Tbx5_selection"
)

pdf("Heatmap_DEGs_Tbx15_transition_top20_top_middle_bottom.pdf")
DoHeatmap(object = heart, features = c("Tbx1", 
                                       rownames(DEGs_GoI_pos_neg_wilcox)[1:20],
                                       rownames(DEGs_GoI_pos_neg_top_wilcox)[1:20],
                                       rownames(DEGs_GoI_pos_neg_bottom_wilcox)[1:20],
                                       "Tbx5"), 
                                        cells = c(Genes_pos,Genes_neg), 
                                        group.by = "Tbx1_Tbx5_selection", 
                                        size = 3) + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

#Heatmap of differentially expression genes (top 20 of pos/neg top/bottom selected cells) as well as Tbx1 and Tbx5 (in total 48 unique genes according overlaps). 
#The selected cells are ordered according their Tbx1 expression (from high to low) and the genes are hierarchically clustered.  
library(gplots)
library(RColorBrewer)

genes_of_interest <- unique(c("Tbx1", 
                              rownames(DEGs_GoI_pos_neg_wilcox)[1:20],
                              rownames(DEGs_GoI_pos_neg_top_wilcox)[1:20],
                              rownames(DEGs_GoI_pos_neg_bottom_wilcox)[1:20],
                              "Tbx5") )
heart.subset <- as.matrix(GetAssayData(heart, slot = "counts"))


#Cells ordered according Tbx1 expression
mat <- heart.subset[genes_of_interest, c(Genes_pos, Genes_neg_top,Genes_neg_bottom)]
col.table <- cbind(colnames(mat),c(rep("red", 13),rep("blue", 16),rep("green", 17)))
mat_ordered <- order(mat["Tbx1",], decreasing = TRUE)
mat_cell_names <- colnames(mat[,mat_ordered])

idx <- sapply(mat_cell_names, function(x) {
  which(col.table[,1] == x)
})

pdf("Heatmap_Tbx1Tbx5_transition_cells_ordered_according_Tbx1.pdf")
genes_of_interest <- unique(c("Tbx1", 
                              rownames(DEGs_GoI_pos_neg_wilcox)[1:20],
                              rownames(DEGs_GoI_pos_neg_top_wilcox)[1:20],
                              rownames(DEGs_GoI_pos_neg_bottom_wilcox)[1:20],
                              "Tbx5") )
mat2 <- heart.subset[genes_of_interest, mat_cell_names]
heatmap.2(mat2, Colv=NULL, scale="row", col=colorRampPalette(c("darkblue", "white","red")),
          sepcolor="gray", sepwidth=c(0.05,0.05), key=TRUE, trace="none", colsep=1:ncol(mat),
          rowsep=1:nrow(mat), cexRow=1, cexCol=1, colCol=c(col.table[idx,2]))
legend(y=1.1, x=.25, xpd=TRUE, legend = c("positive selected cells (n=13)","negative selected cells above (n=13)","negative selected cells below (n=17)"), col = c("red", "blue", "green"), lty= 1, lwd = 5, cex=.7)
dev.off()


#Cells clustered according boundary and top/bottom transition
pdf("Heatmap_Tbx1Tbx5_transition_cells_ordered_according_transition.pdf")
mat2 <- heart.subset[genes_of_interest, c(Genes_neg_bottom, Genes_pos, Genes_neg_top)]
heatmap.2(mat2, Colv=NULL, scale="row", col=colorRampPalette(c("darkblue", "white","red")),
          sepcolor="gray", sepwidth=c(0.05,0.05), key=TRUE, trace="none", colsep=1:ncol(mat),
          rowsep=1:nrow(mat), cexRow=1, cexCol=0.8, colCol = c(rep("green", 17), rep("red", 13),rep("blue", 16)))
legend(y=1.1, x=.25, xpd=TRUE, legend = c("Cells in boundary (n=13)","Cells in top transition (n=13)","Cells in bottom transition (n=17)"), col = c("red", "blue", "green"), lty= 1, lwd = 5, cex=.7)
dev.off()


#######################################################
############### RA analysis (Caudio Cortes) ###########
#######################################################
#Rara, Rarb & Rarg cells+ at E8
heart_E8 <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E8")
RAR_cells <- WhichCells(object = heart_E8, expression = Rara>0 & Rarb>0 & Rarg >0)  
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()

### DEGs analysis: 12 Rar+ cells at E8 vs. 114 E9.posterior Tbx5- cells (no Tbx5-)
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_E9post_noTbx5 <- WhichCells(object = heart_E9.posterior, expression = (Tbx5 == 0))  
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_noTbx5, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx- cells at E9.5 posterior") + NoLegend()
DEGs_Rar_vs_noTbx5 = FindMarkers(heart.raw.data, cells.1 = RAR_cells, cells.2 = cells_E9post_noTbx5, verbose = T)
View(DEGs_Rar_vs_noTbx5)  
write.table(DEGs_Rar_vs_noTbx5, "DEGs_Rar_vs_noTbx5.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells] <- "Rar+"
category[category %in% cells_E9post_noTbx5] <- "Tbx5-"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

#122 DEGs with agj p<0.05
DEGs <- DEGs_Rar_vs_noTbx5[DEGs_Rar_vs_noTbx5$p_val_adj<0.05,]
pdf("Rar_cells_E8_vs_noTbx5E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_noTbx5, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- cells at E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells,cells_E9post_noTbx5), size = 3)
dev.off()

### DEGs analysis: 12 Rar+ cells at E8 vs. 26 E9.posterior Tbx1+ cells (strong expression)
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_E9post_Tbx1 <- WhichCells(object = heart_E9.posterior, expression = (Tbx1 > 1))  
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Tbx1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ (Tbx1>1) cells at E9.5 posterior") + NoLegend()
DEGs_Rar_vs_Tbx1 = FindMarkers(heart.raw.data, cells.1 = RAR_cells, cells.2 = cells_E9post_Tbx1, verbose = T)
View(DEGs_Rar_vs_Tbx1)  
write.table(DEGs_Rar_vs_Tbx1, "DEGs_Rar_vs_strong_Tbx1.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells] <- "Rar+"
category[category %in% cells_E9post_Tbx1] <- "Tbx1+"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

#28 DEGs with agj p<0.05
DEGs <- DEGs_Rar_vs_Tbx1[DEGs_Rar_vs_Tbx1$p_val_adj<0.05,]
pdf("Rar_cells_E8_vs_Tbx1E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Tbx1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ (Tbx1>1) cells at E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(DEGs), group.by = "RA_analysis", cells = c(RAR_cells,cells_E9post_Tbx1), size = 3)
dev.off()

### DEGs analysis: 12 Rar+ cells at E8 vs. 156 E9.posterior Tcf21+ cells
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_E9post_Tcf21 <- WhichCells(object = heart_E9.posterior, expression = (Tcf21 > 0))  
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Tcf21, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tcf21+ cells at E9.5 posterior") + NoLegend()
DEGs_Rar_vs_Tcf21 = FindMarkers(heart.raw.data, cells.1 = RAR_cells, cells.2 = cells_E9post_Tcf21, verbose = T)
View(DEGs_Rar_vs_Tcf21)  
write.table(DEGs_Rar_vs_Tcf21, "DEGs_Rar_vs_Tcf21.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells] <- "Rar+"
category[category %in% cells_E9post_Tcf21] <- "Tcf21+"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

#154 DEGs with agj p<0.05
DEGs <- DEGs_Rar_vs_Tcf21[DEGs_Rar_vs_Tcf21$p_val_adj<0.05,]
pdf("Rar_cells_E8_vs_Tcf21E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Tbx1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tcf21+ cells at E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells,cells_E9post_Tcf21), size = 3)
dev.off()

### DEGs analysis: 12 Rar+ cells at E8 vs. E9.posterior (42)/anterior(26) Dhrs3+/Stra6+ cells
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_E9post_Dhrs3_Stra6 <- WhichCells(object = heart_E9.posterior, expression = (Dhrs3 > 0 & Stra6 >0))  
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Dhrs3_Stra6, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Dhrs3+/Stra6+ cells at E9.5 posterior") + NoLegend()
heart_E9.anterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.anterior")
cells_E9ant_Dhrs3_Stra6 <- WhichCells(object = heart_E9.anterior, expression = (Dhrs3 > 0 & Stra6 >0))  
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9ant_Dhrs3_Stra6, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Dhrs3+/Stra6+ cells at E9.5 anterior") + NoLegend()

#42 DEGs at E9.posterior
DEGs_Rar_vs_Dhrs3_Stra6 = FindMarkers(heart.raw.data, cells.1 = RAR_cells, cells.2 = cells_E9post_Dhrs3_Stra6, verbose = T)
DEGs <- DEGs_Rar_vs_Dhrs3_Stra6[DEGs_Rar_vs_Dhrs3_Stra6$p_val_adj<0.05,]
write.table(DEGs_Rar_vs_Dhrs3_Stra6, "DEGs_Rar_vs_Dhrs3_Stra6.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells] <- "Rar+"
category[category %in% cells_E9post_Dhrs3_Stra6] <- "Dhrs3+/Stra6+"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Rar_cells_E8_vs_Dhrs3_Stra6_E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9post_Dhrs3_Stra6, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Dhrs3+/Stra6+ cells at E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(DEGs), group.by = "RA_analysis", cells = c(RAR_cells,cells_E9post_Dhrs3_Stra6), size = 3)
dev.off()

#10 DEGs at E9.anterior
DEGs_Rar_vs_Dhrs3_Stra6 = FindMarkers(heart.raw.data, cells.1 = RAR_cells, cells.2 = cells_E9ant_Dhrs3_Stra6, verbose = T)
DEGs <- DEGs_Rar_vs_Dhrs3_Stra6[DEGs_Rar_vs_Dhrs3_Stra6$p_val_adj<0.05,]
write.table(DEGs_Rar_vs_Dhrs3_Stra6, "DEGs_Rar_vs_Dhrs3_Stra6.anterior.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells] <- "Rar+"
category[category %in% cells_E9ant_Dhrs3_Stra6] <- "Dhrs3+/Stra6+"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Rar_cells_E8_vs_Dhrs3_Stra6_E9ant.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara, Rarb & Rarg cells+ at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_E9ant_Dhrs3_Stra6, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Dhrs3+/Stra6+ cells at E9.5 anterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(DEGs), group.by = "RA_analysis", cells = c(RAR_cells,cells_E9ant_Dhrs3_Stra6), size = 3)
dev.off()


### 130 RAR_cells_less_string cells at E8 vs. 174 E9 posterior cells
heart_E8 <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E8")
RAR_cells_less_string <- WhichCells(object = heart_E8, expression = Rara>0 | Rarb>0 | Rarg>0)
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")

#906 DEGs
DEGs_lessRarE8_vs_E9post = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = colnames(heart_E9.posterior), verbose = T)
DEGs <- DEGs_lessRarE8_vs_E9post[DEGs_lessRarE8_vs_E9post$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_E9post, "DEGs_lessRarE8_vs_E9post.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+"
category[category %in% colnames(heart_E9.posterior)] <- "E9.posterior"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = colnames(heart_E9.posterior), pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,colnames(heart_E9.posterior)), size = 3)
dev.off()


### 130 RAR_cells_less_string cells at E8 vs. 114 Tbx5- E9 posterior cells
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_noTbx5_E9.posterior <- WhichCells(object = heart_E9.posterior, expression = Tbx5==0)

#815 DEGs
DEGs_lessRarE8_vs_noTbx5E9post = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_noTbx5_E9.posterior, verbose = T)
DEGs <- DEGs_lessRarE8_vs_noTbx5E9post[DEGs_lessRarE8_vs_noTbx5E9post$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_noTbx5E9post, "DEGs_lessRarE8_vs_noTbx5E9post.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_noTbx5_E9.posterior] <- "Tbx5- E9.posterior"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_noTbx5_E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_noTbx5_E9.posterior, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_noTbx5_E9.posterior), size = 3)
dev.off()

### RAR_cells_less_string cells at E8 vs. 103 E9 posterior Tbx5- Tcf21+ cells
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_noTbx5_Tcf21_E9.posterior <- WhichCells(object = heart_E9.posterior, expression = (Tbx5==0 & Tcf21>0))

#742 DEGs
DEGs_lessRarE8_vs_noTbx5_Tcf21E9post = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_noTbx5_Tcf21_E9.posterior, verbose = T)
DEGs <- DEGs_lessRarE8_vs_noTbx5_Tcf21E9post[DEGs_lessRarE8_vs_noTbx5_Tcf21E9post$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_noTbx5_Tcf21E9post, "DEGs_lessRarE8_vs_noTbx5_Tcf21E9post.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_noTbx5_Tcf21_E9.posterior] <- "Tbx5- Tcf21+ E9.posterior"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_noTbx5_Tcf21_E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_noTbx5_Tcf21_E9.posterior, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- Tcf21+ E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_noTbx5_Tcf21_E9.posterior), size = 3)
dev.off()


### RAR_cells_less_string cells at E8 vs. 156 E9 posterior Tcf21+ cells
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_Tcf21_E9.posterior <- WhichCells(object = heart_E9.posterior, expression = (Tcf21>0))

#888 DEGs
DEGs_lessRarE8_vs_Tcf21E9post = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_Tcf21_E9.posterior, verbose = T)
DEGs <- DEGs_lessRarE8_vs_Tcf21E9post[DEGs_lessRarE8_vs_Tcf21E9post$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_Tcf21E9post, "DEGs_lessRarE8_vs_Tcf21E9post.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_Tcf21_E9.posterior] <- "Tcf21+ E9.posterior"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_Tcf21_E9post.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_Tcf21_E9.posterior, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tcf21+ E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_Tcf21_E9.posterior), size = 3)
dev.off()

### RAR_cells_less_string cells at E8 vs. 133 cluster 0 cells
cells_cluster0 <- WhichCells(object = heart, idents = 0)

#1096 DEGs
DEGs_lessRarE8_vs_cluster0 = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_cluster0, verbose = T)
DEGs <- DEGs_lessRarE8_vs_cluster0[DEGs_lessRarE8_vs_cluster0$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_cluster0, "DEGs_lessRarE8_vs_cluster0.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_cluster0] <- "Cluster 0"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_cluster0.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_cluster0, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Cluster 0") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_cluster0), size = 3)
dev.off()

### RAR_cells_less_string cells at E8 vs. 93 cluster 0 Tbx5- cells
cells_cluster0_noTbx5 <- WhichCells(object = heart, expression = (Tbx5==0), idents = 0)

#1007 DEGs
DEGs_lessRarE8_vs_cluster0_noTbx5 = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_cluster0_noTbx5, verbose = T)
DEGs <- DEGs_lessRarE8_vs_cluster0_noTbx5[DEGs_lessRarE8_vs_cluster0_noTbx5$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_cluster0_noTbx5, "DEGs_lessRarE8_vs_cluster0_noTbx5.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_cluster0_noTbx5] <- "Tbx5- cluster 0"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_cluster0_noTbx5.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_cluster0_noTbx5, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- cluster 0") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_cluster0_noTbx5), size = 3)
dev.off()

### RAR_cells_less_string cells at E8 vs. 122 cluster 0 Tcf21+ cells
cells_cluster0_Tcf21 <- WhichCells(object = heart, expression = (Tcf21>0), idents = 0)

#1013 DEGs
DEGs_lessRarE8_vs_cluster0_Tcf21 = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_cluster0_Tcf21, verbose = T)
DEGs <- DEGs_lessRarE8_vs_cluster0_Tcf21[DEGs_lessRarE8_vs_cluster0_Tcf21$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_cluster0_Tcf21, "DEGs_lessRarE8_vs_cluster0_Tcf21.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_cluster0_Tcf21] <- "Tcf21+ cluster 0"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_cluster0_Tcf21.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_cluster0_Tcf21, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tcf21+ cluster 0") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_cluster0_Tcf21), size = 3)
dev.off()

# RAR_cells_less_string cells at E8 vs. 83 cluster 0 Tbx5- Tcf21+ cells
cells_cluster0_noTbx5_Tcf21 <- WhichCells(object = heart, expression = (Tbx5==0 & Tcf21>0), idents = 0)

#906 DEGs
DEGs_lessRarE8_vs_cluster0_noTbx5_Tcf21 = FindMarkers(heart.raw.data, cells.1 = RAR_cells_less_string, cells.2 = cells_cluster0_noTbx5_Tcf21, verbose = T)
DEGs <- DEGs_lessRarE8_vs_cluster0_noTbx5_Tcf21[DEGs_lessRarE8_vs_cluster0_noTbx5_Tcf21$p_val_adj<0.05,]
write.table(DEGs_lessRarE8_vs_cluster0_noTbx5_Tcf21, "DEGs_lessRarE8_vs_cluster0_noTbx5_Tcf21.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% RAR_cells_less_string] <- "Rar+ E8.5"
category[category %in% cells_cluster0_noTbx5_Tcf21] <- "Tbx5- Tcf21+ cluster 0"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "RA_analysis"
)

pdf("Less_Rar_cells_E8_vs_cluster0_noTbx5_Tcf21.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = RAR_cells_less_string, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Rara>0 | Rarb>0 | Rarg>0 cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_cluster0_noTbx5_Tcf21, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- Tcf21+ cluster 0") + NoLegend()
DoHeatmap(object = heart, features = rownames(head(DEGs,50)), group.by = "RA_analysis", cells = c(RAR_cells_less_string,cells_cluster0_noTbx5_Tcf21), size = 3)
dev.off()


#RAR_cells_less_string2 <- WhichCells(object = heart, expression = ( (Rara>0 | Rarb>0 | Rarg>0) & (Rxra>0 | Rxrb>0 | Rxrg>0) ) )
#x <- SubsetData(heart, cells = rownames(heart@meta.data) %in% RAR_cells_less_string2)
#table(x@meta.data[['orig.ident']])


###############################################################################
############# Tbx1 to Tbx 5 trajectory anaylsis (Charlotte & Claudio) #########
###############################################################################

# Tbx1/Tbx5 expression profiles along cluster 0,1,2; combinations of ...
# Set A: E8 Tbx1+ cells (127)
# Set B: E8/9 Tbx1+ cells
# Set C: E9 Tbx5+ cells
# Set D: E9 Tbx1+ and Tbx5+ cells
# Set E: E9 ant Tbx1+ vs. E9 post Tbx5- cells

heart_E8 <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E8")
cells_setA <- WhichCells(object = heart_E8, expression = (Tbx1>0), idents = c(0,1,2)) #104
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setA, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E8.5") + NoLegend()
cells_setB <- WhichCells(object = heart, expression = (Tbx1>0), idents = c(0,1,2)) #291
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setB, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E8.5/E9.5") + NoLegend()

heart_E9 <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.anterior" | heart@meta.data$orig.ident %in% "E9.posterior")
cells_setC <- WhichCells(object = heart_E9, expression = (Tbx5>0), idents = c(0,1,2)) #72
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setC, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5+ cells at E9.5") + NoLegend()
cells_setD <- WhichCells(object = heart_E9, expression = (Tbx5>0 & Tbx1>0), idents = c(0,1,2)) #52
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setD, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ & Tbx5+ cells at E9.5") + NoLegend()

heart_E9.anterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.anterior")
heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior")
cells_setE_1 <- WhichCells(object = heart_E9.anterior, expression = (Tbx1>0), idents = c(0,1,2)) #33
cells_setE_2 <- WhichCells(object = heart_E9.posterior, expression = (Tbx5==0), idents = c(0,1,2)) #112
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E9.5 anterior") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_2, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- cells at E9.5 posterior") + NoLegend()

#44 DEGs
DEGs_setE <- FindMarkers(heart.raw.data, cells.1 = cells_setE_1, cells.2 = cells_setE_2, verbose = T)
DEGs <- DEGs_setE[DEGs_setE$p_val_adj<0.05,]
write.table(DEGs_setE, "DEGs_setE.txt", sep="\t", quote = F)

category <- colnames(heart)
category[category %in% cells_setE_1] <- "Tbx1+ E9.5 anterior"
category[category %in% cells_setE_2] <- "Tbx5- E9.5 posterior"
heart <- AddMetaData(
  object = heart,
  metadata = category,
  col.name = "Tbx_analysis"
)

pdf("Tbx1+E9anteriorvs_Tbx5-E9posterior_cluster012.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E9.5 anterior") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_2, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- cells at E9.5 posterior") + NoLegend()
DoHeatmap(object = heart, features = rownames(DEGs), group.by = "Tbx_analysis", cells = c(cells_setE_1,cells_setE_2), size = 3)
dev.off()

pdf("Tbx1_Tbx5_expression_cluster012.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setA, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E8.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setB, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E8.5/E9.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setC, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5+ cells at E9.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setD, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ & Tbx5+ cells at E9.5") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_1, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx1+ cells at E9.5 anterior") + NoLegend()
DimPlot(object = heart, reduction="tsne", cells.highlight = cells_setE_2, pt.size = 2, sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "Tbx5- cells at E9.5 posterior") + NoLegend()
dev.off()

Tbx1_Tbx5_cells_cluster012 <- unique(c(cells_setA,cells_setB,cells_setC,cells_setD,cells_setE_1, cells_setE_2)) #315 cells
Tbx1_Tbx5_cells_cluster012_object <- SubsetData(heart, cells = rownames(heart@meta.data) %in% Tbx1_Tbx5_cells_cluster012)

heart_cluster012.markers <- FindAllMarkers(object = Tbx1_Tbx5_cells_cluster012_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
heart_cluster012.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC)
write.table(heart_cluster012.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_logFC), "top20_genes_per_cluster.txt", quote=F, sep="\t")


######################################################################
############# Tbx1 to Tbx 5 trajectory anaylsis (Charlotte ) #########
######################################################################
### Candidate genes with similar Tbx1 and Tbx5 expression profile

# 21,769 features
# 462 cells
heart.subset <- as.matrix(GetAssayData(heart, slot = "counts"))
M <- t(heart.subset)
which( colnames(M)=="Tbx1" )
out <- t(cor(M[, 1047], M[, -1047]))

#Pearson correlation
Tbx1_correlations <- as.data.frame(`colnames<-`(cbind(out, out^2), c("cor", "r2")))
head(Tbx1_correlations[order(Tbx1_correlations$cor, decreasing = T),],10)

# cor        r2
# Nnat    0.6447345 0.4156826
# Gm266   0.6263209 0.3922778
# Col13a1 0.6150451 0.3782805
# Ebf1    0.5984943 0.3581955
# Syt17   0.5693655 0.3241771
# Eif2b2  0.5680130 0.3226388
# Pgf     0.5513607 0.3039986
# Pknox2  0.5206642 0.2710912
# Serinc2 0.5131601 0.2633333
# Slit1   0.5112133 0.2613391

pdf("Top3_candidates_similar_to_Tbx1.pdf")
plot(M[,"Tbx1"],M[,"Nnat"], xlab = "Tbx1 expression", ylab = "Nnat expression")
lines(lowess(M[,"Tbx1"],M[,"Nnat"], f=.2), col="red")
legend("topleft", "R=0.645")
FeaturePlot(object = heart, features = c("Tbx1","Nnat"), reduction="tsne",pt.size = 1.5)
plot(M[,"Tbx1"],M[,"Gm266"], xlab = "Tbx1 expression", ylab = "Gm266 expression")
lines(lowess(M[,"Tbx1"],M[,"Gm266"], f=.2), col="red")
legend("topleft", "R=0.626")
FeaturePlot(object = heart, features = c("Tbx1","Gm266"), reduction="tsne",pt.size = 1.5)
plot(M[,"Tbx1"],M[,"Col13a1"], xlab = "Tbx1 expression", ylab = "Col13a1 expression")
lines(lowess(M[,"Tbx1"],M[,"Col13a1"], f=.2), col="red")
legend("topleft", "R=0.615")
FeaturePlot(object = heart, features = c("Tbx1","Col13a1"), reduction="tsne",pt.size = 1.5)
dev.off()


which( colnames(M)=="Tbx5" )
out <- t(cor(M[, 1613], M[, -1613]))

#Pearson correlation
Tbx5_correlations <- as.data.frame(`colnames<-`(cbind(out, out^2), c("cor", "r2")))
head(Tbx5_correlations[order(Tbx5_correlations$cor, decreasing = T),],10)

# cor        r2
# Wnt2    0.6691759 0.4477964
# Sfrp5   0.4891877 0.2393046
# Osr1    0.4844789 0.2347198
# Gata4   0.4788422 0.2292899
# Fitm1   0.4761432 0.2267123
# Gm5463  0.4686606 0.2196427
# Fbxo32  0.4575572 0.2093586
# Bhlhe41 0.4325392 0.1870902
# Gm38158 0.4252043 0.1807987
# Gm16046 0.4238539 0.179652

pdf("Top3_candidates_similar_to_Tbx5.pdf")
plot(M[,"Tbx5"],M[,"Wnt2"], xlab = "Tbx5 expression", ylab = "Nnat expression")
lines(lowess(M[,"Tbx5"],M[,"Wnt2"], f=.2), col="red")
legend("top", "R=0.669")
FeaturePlot(object = heart, features = c("Tbx5","Wnt2"), reduction="tsne",pt.size = 1.5)
plot(M[,"Tbx5"],M[,"Sfrp5"], xlab = "Tbx5 expression", ylab = "Sfrp5 expression")
lines(lowess(M[,"Tbx5"],M[,"Sfrp5"], f=.2), col="red")
legend("topleft", "R=0.489")
FeaturePlot(object = heart, features = c("Tbx1","Sfrp5"), reduction="tsne",pt.size = 1.5)
plot(M[,"Tbx5"],M[,"Osr1"], xlab = "Tbx5 expression", ylab = "Osr1 expression")
lines(lowess(M[,"Tbx5"],M[,"Osr1"], f=.2), col="red")
legend("topleft", "R=0.484")
FeaturePlot(object = heart, features = c("Tbx5","Osr1"), reduction="tsne",pt.size = 1.5)
dev.off()

### Trajectory analysis based on cluster 0,1,2 (n=353 cells)

library(monocle)

#function to import Seurat object 
newimport <- function(otherCDS, import_all = TRUE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- otherCDS@assays$RNA@counts
    
    if(class(data) == "data.frame") {
      data <- as(as.matrix(data), "sparseMatrix")
    }
    
    pd <- tryCatch( {
      pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
      pd
    }, 
    #warning = function(w) { },
    error = function(e) { 
      pData <- data.frame(cell_id = colnames(data), row.names = colnames(data))
      pd <- new("AnnotatedDataFrame", data = pData)
      
      message("This Seurat object doesn't provide any meta data");
      pd
    })
    
    # remove filtered cells from Seurat
    if(length(setdiff(colnames(data), rownames(pd))) > 0) {
      data <- data[, rownames(pd)]  
    }
    
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    lowerDetectionLimit <- 0
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    valid_data <- data[, row.names(pd)]
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
        
      } else {
        # mist_list <- list(ident = ident) 
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    }
    
    if(1==1) {
      var.genes <- setOrderingFilter(monocle_cds, otherCDS@assays$RNA@var.features)
      
    }
    monocle_cds@auxClusteringData$seurat <- mist_list
    
  } else if (class(otherCDS)[1] == 'SCESet') {
    requireNamespace("scater")
    
    message('Converting the exprs data in log scale back to original scale ...')    
    data <- 2^otherCDS@assayData$exprs - otherCDS@logExprsOffset
    
    fd <- otherCDS@featureData
    pd <- otherCDS@phenoData
    experimentData = otherCDS@experimentData
    if("is.expr" %in% slotNames(otherCDS))
      lowerDetectionLimit <- otherCDS@is.expr
    else 
      lowerDetectionLimit <- 1
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    if(import_all) {
      # mist_list <- list(iotherCDS@sc3,
      #                   otherCDS@reducedDimension)
      mist_list <- otherCDS 
      
    } else {
      mist_list <- list()
    }
    
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    # monocle_cds@auxClusteringData$sc3 <- otherCDS@sc3
    # monocle_cds@auxOrderingData$scran <- mist_list
    
    monocle_cds@auxOrderingData$scran <- mist_list
    
  } else {
    stop('the object type you want to export to is not supported yet')
  }
  
  return(monocle_cds)
}


heart_cluster012 <- SubsetData(heart, ident.use = c(0,1,2))
DimPlot(object = heart_cluster012, reduction="tsne") + ggtitle(label = "Cluster 0, 1, 2")

#317/353 cells with Tbx1 or Tbx5 expression in cluster 0,1 and 2
cells_with_Tbx1_Tbx5 <- WhichCells(object = heart_cluster012, expression = Tbx1 >0 | Tbx5 >0)
DimPlot(object = heart_cluster012, reduction="tsne", cells.highlight = cells_with_Tbx1_Tbx5, pt.size = 1, sizes.highlight = 1, cols.highlight = "red") + ggtitle(label = "Tbx1/Tbx5 boundary & transition")

heart_Tbx1_Tbx5_cluster012 <- SubsetData(heart, cells = cells_with_Tbx1_Tbx5)
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne")

#create monocle object from Seurat object
monocle <- newimport(heart_Tbx1_Tbx5_cluster012)

var_genes <- heart_Tbx1_Tbx5_cluster012[["RNA"]]@var.features
ordering_genes <- var_genes  
monocle <- setOrderingFilter(monocle, ordering_genes)  
monocle <- reduceDimension(monocle, reduction_method = "DDRTree", max_components = 2)
monocle <- orderCells(monocle)

monocle <- detectGenes(monocle, min_expr = 0.1)
fData(monocle)$use_for_ordering <- fData(monocle)$num_cells_expressed > 0.05 * ncol(monocle)
monocle_expressed_genes <- row.names(subset(fData(monocle), use_for_ordering==TRUE))


Genes_pos <- c("E9.posterior_4_Set_B_C6_S138", "E9.posterior_4_Set_B_H1_S103", "E9.posterior_4_Set_B_A11_S175", 
               "E9.posterior_3_Set_A_A9_S65", "E9.posterior_4_Set_B_A4_S120", "E9.posterior_4_Set_B_H2_S111", 
               "E9.posterior_3_Set_A_H4_S32","E9.posterior_4_Set_B_C1_S98", 
               "E9.posterior_3_Set_A_H3_S24", "E9.posterior_3_Set_A_C11_S83", 
               "E9.posterior_3_Set_A_C6_S43", "E9.posterior_3_Set_A_F10_S78", "E9.posterior_3_Set_A_D4_S28")

Genes_neg <- c("E9.anterior_2_Set_B_F10_S174", "E9.posterior_3_Set_A_H2_S16", "E9.posterior_3_Set_A_A10_S73", "E9.posterior_4_Set_B_A1_S96",
               "E9.posterior_4_Set_B_D1_S99", "E9.posterior_4_Set_B_B6_S137", "E9.posterior_3_Set_A_G12_S95",
               "E9.posterior_3_Set_A_D2_S12", "E9.posterior_3_Set_A_E11_S85", "E9.posterior_3_Set_A_A4_S25", "E9.posterior_3_Set_A_E2_S13",
               "E9.posterior_3_Set_A_G5_S39", "E9.posterior_3_Set_A_B12_S90", "E9.posterior_3_Set_A_C4_S27",
               "E9.posterior_4_Set_B_A7_S144", "E9.posterior_3_Set_A_E1_S5",
               "E9.posterior_4_Set_B_C12_S185","E9.posterior_4_Set_B_A10_S167","E9.posterior_3_Set_A_H1_S8","E9.posterior_4_Set_B_H12_S190",
               "E9.posterior_4_Set_B_H5_S135", "E9.posterior_3_Set_A_B11_S82", "E9.posterior_4_Set_B_C11_S177", "E9.posterior_3_Set_A_F9_S70",
               "E9.posterior_3_Set_A_F3_S22", "E9.posterior_3_Set_A_G11_S87", "E9.posterior_3_Set_A_F7_S54", "E9.posterior_3_Set_A_E6_S45",
               "E9.posterior_4_Set_B_D11_S178", "E9.posterior_4_Set_B_H10_S174",
               "E9.posterior_3_Set_A_C12_S91", "E9.posterior_3_Set_A_G3_S23", "E9.posterior_4_Set_B_E12_S187")

Genes_neg_top <- c("E9.posterior_3_Set_A_A4_S25", "E9.posterior_3_Set_A_E2_S13","E9.posterior_3_Set_A_G5_S39", "E9.posterior_3_Set_A_B12_S90", "E9.posterior_3_Set_A_C4_S27",
                   "E9.posterior_4_Set_B_A7_S144", "E9.posterior_3_Set_A_E1_S5","E9.posterior_4_Set_B_C12_S185","E9.posterior_4_Set_B_A10_S167","E9.posterior_3_Set_A_H1_S8",
                   "E9.posterior_4_Set_B_D11_S178", "E9.posterior_4_Set_B_H10_S174","E9.posterior_4_Set_B_E12_S187","E9.posterior_4_Set_B_H12_S190","E9.posterior_3_Set_A_E6_S45",
                   "E9.posterior_3_Set_A_F7_S54")

Genes_neg_bottom <- c("E9.anterior_2_Set_B_F10_S174", "E9.posterior_3_Set_A_H2_S16", "E9.posterior_3_Set_A_A10_S73", "E9.posterior_4_Set_B_A1_S96",
                      "E9.posterior_4_Set_B_D1_S99", "E9.posterior_4_Set_B_B6_S137", "E9.posterior_3_Set_A_G12_S95","E9.posterior_3_Set_A_D2_S12",
                      "E9.posterior_4_Set_B_H5_S135", "E9.posterior_3_Set_A_B11_S82", "E9.posterior_4_Set_B_C11_S177", "E9.posterior_3_Set_A_F9_S70",
                      "E9.posterior_3_Set_A_C12_S91", "E9.posterior_3_Set_A_G3_S23", "E9.posterior_3_Set_A_E11_S85", "E9.posterior_3_Set_A_F3_S22",
                      "E9.posterior_3_Set_A_G11_S87")

heart.raw.data <- as.matrix(GetAssayData(heart, slot = "counts"))
DEGs_GoI_pos_neg_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg, verbose = T)
DEGs_GoI_pos_neg_top_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg_top, verbose = T)
DEGs_GoI_pos_neg_bottom_wilcox = FindMarkers(heart.raw.data, cells.1 = Genes_pos, cells.2 = Genes_neg_bottom, verbose = T)


#115 genes
monocle_ordering_genes <- unique(c("Tbx1", 
                            rownames(DEGs_GoI_pos_neg_wilcox)[1:50],
                            rownames(DEGs_GoI_pos_neg_top_wilcox)[1:50],
                            rownames(DEGs_GoI_pos_neg_bottom_wilcox)[1:50],
                            "Tbx5"))

monocle <- setOrderingFilter(monocle, ordering_genes = monocle_ordering_genes)
monocle <- reduceDimension(monocle, method = 'DDRTree', max_components = 2)
monocle <- orderCells(monocle)

pdf("Tbx1_Tbx5_trajectory_115DEGs.pdf")
plot_cell_trajectory(monocle, color_by = 'seurat_clusters')
plot_cell_trajectory(monocle, color_by = 'orig.ident')
plot_cell_trajectory(monocle, color_by = 'Pseudotime')
plot_cell_trajectory(monocle, color_by = 'State')
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Tbx1", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Tbx5", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Gdnf", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Klf16", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Cldn11", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Foxf2", markers_linear = T)
plot_cell_trajectory(monocle, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Aldh1a2", markers_linear = T)
dev.off()

monocle_states <- cbind(sampleNames(monocle), unlist(monocle$State))

pdf("Tbx1_Tbx5_trajectory_states.pdf")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", cells.highlight = monocle_states[monocle_states[,2] %in% "1",][,1], sizes.highlight = 2, cols.highlight = "red") + ggtitle(label = "State 1")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", cells.highlight = monocle_states[monocle_states[,2] %in% "2",][,1], sizes.highlight = 2, cols.highlight = "bisque4") + ggtitle(label = "State 2")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", cells.highlight = monocle_states[monocle_states[,2] %in% "3",][,1], sizes.highlight = 2, cols.highlight = "darkolivegreen4") + ggtitle(label = "State 3")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", cells.highlight = monocle_states[monocle_states[,2] %in% "4",][,1], sizes.highlight = 2, cols.highlight = "blue") + ggtitle(label = "State 4")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", cells.highlight = monocle_states[monocle_states[,2] %in% "5",][,1], sizes.highlight = 2, cols.highlight = "darkorchid3") + ggtitle(label = "State 5")
dev.off()

category <- colnames(heart_Tbx1_Tbx5_cluster012)
category[category %in% c(monocle_states[monocle_states[,2] %in% "1",][,1])] <- "State 1"
category[category %in% c(monocle_states[monocle_states[,2] %in% "2",][,1])] <- "State 2"
category[category %in% c(monocle_states[monocle_states[,2] %in% "3",][,1])] <- "State 3"
category[category %in% c(monocle_states[monocle_states[,2] %in% "4",][,1])] <- "State 4"
category[category %in% c(monocle_states[monocle_states[,2] %in% "5",][,1])] <- "State 5"

heart_Tbx1_Tbx5_cluster012 <- AddMetaData(
  object = heart_Tbx1_Tbx5_cluster012,
  metadata = category,
  col.name = "Trajectory"
)

pdf("Tbx1_Tbx5_trajectory_states_in_original_clustering.pdf")
DimPlot(object = heart_Tbx1_Tbx5_cluster012, reduction="tsne", group.by = "Trajectory", pt.size = 2, label = T)
dev.off()

pdf("Tbx1_Tbx5_trajectory_genes_in_pseudotime.pdf")
plot_genes_in_pseudotime(monocle[c("Tbx1","Tbx5"),])
plot_genes_in_pseudotime(monocle[c("Tbx1","Tbx5", "Gdnf", "Klf16", "Cldn11", "Foxf2", "Aldh1a2"),])
dev.off()

#BEAM takes as input a CellDataSet that's been ordered with orderCells and the name of a branch point in the trajectory. 
#It returns a table of significance scores for each gene. Genes that score significant are said to be branch-dependent in their expression. 

BEAM_res <- BEAM(monocle, branch_point = 1, cores = 3)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
BEAM_res_sub <- subset(BEAM_res[monocle_ordering_genes,], qval < 5e-2)
write.table(BEAM_res_sub, "Tbx1_Tbx5_trajectory_DEGSs_branchpoint1.txt", quote = F, sep = "\t")

BEAM_res_2 <- BEAM(monocle, branch_point = 2, cores = 3)
BEAM_res_2 <- BEAM_res_2[order(BEAM_res_2$qval),]
BEAM_res_2 <- BEAM_res_2[,c("gene_short_name", "pval", "qval")]
BEAM_res_2_sub <- subset(BEAM_res_2[monocle_ordering_genes,], qval < 5e-2)
write.table(BEAM_res_2_sub, "Tbx1_Tbx5_trajectory_DEGSs_branchpoint2.txt", quote = F, sep = "\t")


#You can visualize changes for all the genes that are significantly branch dependent using a special type of heatmap. 
#This heatmap shows changes in both lineages at the same time. It also requires that you choose a branch point to inspect. 
#Columns are points in pseudotime, rows are genes, and the beginning of pseudotime is in the middle of the heatmap. 
#As you read from the middle of the heatmap to the right, you are following one lineage through pseudotime. 
#As you read left, the other. The genes are clustered hierarchically, so you can visualize modules of genes that have similar lineage-dependent expression patterns. 

plot_genes_branched_heatmap(monocle[row.names(BEAM_res_sub),], branch_point = 1, num_clusters = 3, cores = 3, use_gene_short_name = T, show_rownames = T)
plot_genes_branched_heatmap(monocle[row.names(BEAM_res_2_sub),], branch_point = 2, num_clusters = 3, cores = 3, use_gene_short_name = T, show_rownames = T)



#### E8.5 state 5 and 4 to E9.5 posterior Tbx5+ state 4 cells ####
cells_state4 <- monocle_states[monocle_states[,2] %in% "4",][,1] #82
cells_state5 <- monocle_states[monocle_states[,2] %in% "5",][,1] #93
cells_E8_state4_state5 <- unique(c(cells_state4[as.vector(grep('^E8_', cells_state4))], cells_state5[as.vector(grep('^E8_',cells_state5))])) #54

heart_E9.posterior <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% "E9.posterior") #174
cells_E9post_Tbx5 <- WhichCells(object = heart_E9.posterior, expression = (Tbx5>0)) #60
cells_E9post_Tbx5_state4 <- intersect(cells_state4,cells_E9post_Tbx5) #29

cells_of_interest <- unique(c(cells_E8_state4_state5,cells_E9post_Tbx5_state4)) #83

pdf("E8_state5_and_4_to_E9posterior_Tbx5pos_state4_cells.pdf")
DimPlot(object = heart, reduction="tsne", cells.highlight =  cells_E8_state4_state5, pt.size = 1, sizes.highlight = 1, cols.highlight = c("blue")) + 
  ggtitle(label = "E8.5 state 5 and 4 cells") + NoLegend()

DimPlot(object = heart, reduction="tsne", cells.highlight =  cells_E9post_Tbx5_state4, pt.size = 1, sizes.highlight = 1, cols.highlight = c("darkgreen")) + 
  ggtitle(label = "E9.5 posterior Tbx5+ state 4 cells") + NoLegend()

DimPlot(object = heart, reduction="tsne", cells.highlight =  cells_of_interest, pt.size = 1, sizes.highlight = 1, cols.highlight = c("red")) + 
  ggtitle(label = "E8.5 state 5 and 4 to E9.5 posterior Tbx5+ state 4 cells") + NoLegend()
dev.off()

#subset
heart_subset <- SubsetData(heart, cells=cells_of_interest)

monocle_subset <- newimport(heart_subset)
var_genes <- heart_subset[["RNA"]]@var.features
ordering_genes <- var_genes 

monocle_subset <- setOrderingFilter(monocle_subset, ordering_genes)  
monocle_subset <- reduceDimension(monocle_subset, reduction_method = "DDRTree", max_components = 2)
monocle_subset <- orderCells(monocle_subset)

monocle_subset <- setOrderingFilter(monocle_subset, ordering_genes = monocle_ordering_genes)
monocle_subset <- reduceDimension(monocle_subset, method = 'DDRTree', max_components = 2)
monocle_subset <- orderCells(monocle_subset)

pdf("Tbx1_Tbx5_trajectory_115DEGs_subset.pdf")
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters')
plot_cell_trajectory(monocle_subset, color_by = 'orig.ident')
plot_cell_trajectory(monocle_subset, color_by = 'Pseudotime')
plot_cell_trajectory(monocle_subset, color_by = 'State')
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Tbx1", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Tbx5", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Gdnf", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Klf16", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Cldn11", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Foxf2", markers_linear = T)
plot_cell_trajectory(monocle_subset, color_by = 'seurat_clusters', show_state_number = F, cell_size = 2, markers = "Aldh1a2", markers_linear = T)
dev.off()


### tSNE plot of the E9.5 data (anterior and posterior) ###
heart_E9 <- SubsetData(heart, cells = heart@meta.data$orig.ident %in% c("E9.posterior","E9.anterior")) #307 cells
rownames(heart_E9@meta.data)
heart_E9_clustering <- CreateSeuratObject(data, min.cells = 3, min.features = 200)
heart_E9_clustering <- SubsetData(heart_E9_clustering, cells = rownames(heart_E9_clustering@meta.data) %in% rownames(heart_E9@meta.data))
heart_E9_clustering <- NormalizeData(object = heart_E9_clustering)
heart_E9_clustering <- FindVariableFeatures(object = heart_E9_clustering, selection.method = "vst", nfeatures = 500)
all.genes_E9 <- rownames(x=heart_E9_clustering)
heart_E9_clustering <- ScaleData(object = heart_E9_clustering, features = all.genes_E9)
heart_E9_clustering <- RunPCA(object = heart_E9_clustering, features = VariableFeatures(object = heart_E9_clustering))
heart_E9_clustering <- JackStraw(object = heart_E9_clustering, num.replicate = 100)
heart_E9_clustering <- ScoreJackStraw(object = heart_E9_clustering, dims = 1:20)
JackStrawPlot(object = heart_E9_clustering, dims = 1:20)
heart_E9_clustering <- RunUMAP(object = heart_E9_clustering, reduction="pca", dims = 1:9)
heart_E9_clustering <- FindNeighbors(object = heart_E9_clustering, reduction = "pca", dims=1:9)
heart_E9_clustering <- FindClusters(heart_E9_clustering, resolution = 0.5)
heart_E9_clustering <- RunUMAP(object = heart_E9_clustering, dims = 1:9)
heart_E9_clustering <- RunTSNE(object = heart_E9_clustering, dims = 1:9)

pdf("Subclusering_of_E9_anterior_and_posterior.pdf")
DimPlot(object = heart_E9_clustering, reduction="tsne", label = TRUE, pt.size=2) #4 cluster
DimPlot(object = heart_E9_clustering, reduction="tsne", label = TRUE, pt.size=2, group.by = "orig.ident") #4 cluster
dev.off()

#cell browser output
heart_E9_clustering@misc$markers <- FindAllMarkers(heart_E9_clustering)
saveRDS(heart_E9_clustering, "E9_cells.rds")



###################### Marker genes #############################


CPM_genes <- c("Nrg1","Mef2c","Six1","Fgf10","Six2","Kdr","Kazald1","Isl1","Tbx1","Arg1","Tcf21","Hoxa1","Hoxb5","Aplnr","Prdm1","Six2","Meox1","Vamp5","Cd40")
SHF_genes <- c("Hopx","Mef2c","Prox1","Six1","Fgf10","Six2","Tlx1","Hoxb1","Osr1","Nkx2-5","Isl1","Tbx1","Arg1","Tcf21","Hoxa1","Hoxb5","Aplnr","Prdm1","Six2","Meox1","Vamp5","Cd40","Dach1")
HT_genes <- c("Hopx","Nrg1","Mef2c","Prox1","Pitx2","Shox2","Osr1","Gata6","Nkx2-5","Isl1","Tbx1","Hoxa1","Six2","Ppm1h","Bmp7","Bmp4","Acta2","Pitx2")
BM_genes <- c("Mef2c","Prox1","Six1","Ebf1","Ebf2","Edn1","Hgf","Frzb","Chrdl1","Myf5","Pitx2","Lrrn1","Six2","Lhx2","Six4","Tlx1","Isl1","Tbx1","Tcf21","Six2","Ppm1h","Acta2","Myog","MyoD","Pitx2","Eya2","En2")
CT_genes <- c("Scx","Sox6","Sox9","Mroh2a","Rbfa","Ocel1","Grina","Ppm1h","Hpn","Slc6a7","Dcaf17","Kcnj11","Ptprg","Dll4","Qprt","Car8","Egr1","Tcf4")

markers <- (unique(c(CPM_genes, SHF_genes, HT_genes, BM_genes, CT_genes)))

pdf("All_marker_genes_Seurat_cluster.pdf", height = 12, width = 24)
DotPlot(heart, features=markers) + RotatedAxis()
dev.off()

pdf("All_markers_ridge_plots.pdf")
RidgePlot(heart, features=c("Nrg1",   "Mef2c" ,  "Six1" ,   "Fgf10"  , "Six2" ,   "Kdr"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Kazald1", "Isl1",    "Tbx1",    "Arg1"   , "Tcf21" ,  "Hoxa1"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Hoxb5"   ,"Aplnr"  , "Prdm1",   "Meox1"  , "Vamp5","Cd40"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Hopx"  ,  "Prox1" ,  "Tlx1"  ,  "Hoxb1"  , "Osr1"   , "Nkx2-5"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Dach1"  , "Pitx2"  , "Shox2"  , "Gata6" ,  "Ppm1h"   ,"Bmp7"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Bmp4"  ,  "Acta2"  , "Ebf1"  ,  "Ebf2", "Edn1"  ,  "Hgf"), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Frzb"   , "Chrdl1" , "Myf5" ,   "Lrrn1" ,  "Lhx2" ,   "Six4"  ), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Myog"  ,  "MyoD"  ,  "Eya2"  ,  "En2" ,    "Scx"   ,  "Sox6" ), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Sox9"   , "Mroh2a"  ,"Rbfa" ,"Ocel1" ,  "Grina" ,  "Hpn"  ), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Slc6a7"  ,"Dcaf17"  ,"Kcnj11",  "Ptprg"  , "Dll4" ,   "Qprt" ), ncol = 2, nrow=3)
RidgePlot(heart, features=c("Car8"    ,"Egr1",    "Tcf4"    ), ncol = 2, nrow=3)
dev.off()


############# Velocity analysis #######################

### export for velocity analysis ###

write.csv(Cells(heart), file = "cellID_obs.csv")
write.csv(Embeddings(heart, reduction = "tsne"), file = "cell_embeddings.csv")
write.csv(heart@meta.data$seurat_clusters, file = "clusters.csv")
saveRDS(Embeddings(heart, reduction = "tsne"), file = "cell_embeddings.rds")

####################################

#http://pklab.med.harvard.edu/velocyto/notebooks/R/chromaffin2.nb.html

#load lib
library(velocyto.R)

#load loom files
ldat <- read.loom.matrices("onefilepercell_E8_1_Set_A_A10_S72_and_others_E6ZKY.loom")

#reduce the cell names to the corresponding labels
ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub(".bam","",gsub(".*:","",colnames(x)))
  x
})

#Spliced expression magnitude distribution across genes:
png("../velocity/read_counts_per_gene.png")
hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
dev.off()

#cell colors
dplot <- DimPlot(object = heart, reduction="tsne", label = TRUE, pt.size=3)
pbuild <- ggplot2::ggplot_build(dplot)
pdata <- pbuild$data[[1]]
cell.colors = as.array(pdata$colour)
rownames(cell.colors) <- rownames(dplot$data)

write.csv(cell.colors, "seurat_colors.csv")


emb <- Embeddings(heart, reduction = "tsne")

# exonic read (spliced) expression matrix
emat <- ldat$spliced;
# intronic read (unspliced) expression matrix
nmat <- ldat$unspliced
# spanning read (intron+exon) expression matrix
smat <- ldat$spanning;
# filter expression matrices based on some minimum max-cluster averages
emat <- filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 5)
nmat <- filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
smat <- filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
# look at the resulting gene set ([1] 9998)
length(intersect(rownames(emat),rownames(nmat)))

# and if we use spanning reads (smat)
length(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

#Using min/max quantile fit, in which case gene-specific offsets do not require spanning read (smat) fit. Here the fit is based on the top/bottom 5% of cells (by spliced expression magnitude).
fit.quantile <- 0.05;
rvel.qf <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,kCells = 5,fit.quantile = fit.quantile)

#visualize the velocities by projecting observed and extrapolated cells onto the first 5 PC
png("../velocity/velocity_estimates_gene-relative_model_5PCs.png")
pdf("../velocity/velocity_estimates_gene-relative_model_5PCs.pdf")
pca.velocity.plot(rvel.qf, nPcs=5, plot.cols=2, cell.colors=ac(cell.colors, alpha = 0.7), cex=1.2, pcount=0.1, pc.multipliers=c(1,-1,-1,-1,-1))
dev.off()

#Fitting of individual genes can be visualized using “show.gene” option. To save time, we’ll pass previously-calculated velocity (rvel.qf) to save calculation time:
png("../velocity/velocity_estimates_gene-relative_model_5PCs.Tbx1.png")

png("../velocity/velocity_estimates_gene-relative_model_5PCs.tbx1.png", width = 1600, height = 600)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Tbx1',cell.emb=emb,cell.colors=cell.colors)
dev.off()

pdf("../velocity/velocity_estimates_gene-relative_model_5PCs.GenesOfInterest.pdf", width = 16, height = 6)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Tbx1',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Tbx5',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Aldh1a2',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Fgf10',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Mef2c',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Hoxb1',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Nkx2-5',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Nrg1',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Gli1',cell.emb=emb,cell.colors=cell.colors)
gene.relative.velocity.estimates(emat,nmat, kCells = 5,fit.quantile = fit.quantile,old.fit=rvel.qf,show.gene='Cldn11',cell.emb=emb,cell.colors=cell.colors)
dev.off()


#Alternatively, wen ca use spanning reads (smat) to fit the gene offsets. This will result in more accurate offset estimates, but for much fewer genes (spanning reads are rare). 
rvel <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 5, fit.quantile=fit.quantile, diagonal.quantiles = TRUE)
rvel_k10 <- gene.relative.velocity.estimates(emat,nmat,smat=smat, kCells = 10, fit.quantile=fit.quantile, diagonal.quantiles = TRUE)

#visualize the velocity in PCA space
png("../velocity/velocity_estimates_gene-relative_model_5PCs.smat.png")
pdf("../velocity/velocity_estimates_gene-relative_model_5PCs.smat.pdf")
pca.velocity.plot(rvel,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))
dev.off()

#calculate the most basic version of velocity estimates, using relative gamma fit, without cell kNN smoothing (i.e. actual single-cell velocity):
rvel1 <- gene.relative.velocity.estimates(emat,nmat,deltaT=1,deltaT2 = 1,kCells = 1, fit.quantile=fit.quantile)
png("../velocity/velocity_estimates_gene-relative_model_5PCs.no_kNN.png")
#pdf("../velocity/velocity_estimates_gene-relative_model_5PCs.no_kNN.pdf")
pca.velocity.plot(rvel1,nPcs=5,plot.cols=2,cell.colors=ac(cell.colors,alpha=0.7),cex=1.2,pcount=0.1,pc.multipliers=c(1,-1,1,1,1))
dev.off()

#Visualization on an existing embedding (we use t-SNE embedding from our Seurat analysis)
vel <- rvel; arrow.scale=3; cell.alpha=0.4; cell.cex=1; fig.height=4; fig.width=4.5;
png("../velocity/velocity_estimates.Seurat_custer.embessing.png")
pdf("../velocity/velocity_estimates.eurat_custer.embessing.pdf")
show.velocity.on.embedding.cor(emb,vel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
dev.off()

png("../velocity/velocity_estimates.Seurat_custer.embessing.n10.png")
#pdf("../velocity/velocity_estimates.eurat_custer.embessing.n10.pdf")
show.velocity.on.embedding.cor(emb,vel,n=10,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,arrow.lwd=1)
dev.off()

#Alternatively, the same function can be used to calculate a velocity vector field
png("../velocity/velocity_estimates.Seurat_custer.vector_field.png")
#pdf("../velocity/velocity_estimates.eurat_custer.vector_field.pdf")
show.velocity.on.embedding.cor(emb,vel,n=100,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=20,arrow.lwd=2)
dev.off()

png("../velocity/velocity_estimates.Seurat_custer.vector_field.n10.png")
pdf("../velocity/velocity_estimates.eurat_custer.vector_field.n10.pdf")
show.velocity.on.embedding.cor(emb,vel,n=10,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,arrow.scale=arrow.scale,show.grid.flow=TRUE,min.grid.cell.mass=0.5,grid.n=20,arrow.lwd=2)
dev.off()



#Velocity estimate based on gene structure
#ip.mm10 <- readRDS(url("http://pklab.med.harvard.edu/velocyto/chromaffin/ip.mm10.rds"))
#read velocyto.py HDF5 output
#gene.info <- read.gene.mapping.info("../velocity/dump/onefilepercell_E8_1_Set_A_A10_S72_and_others_E6ZKY.hdf5",internal.priming.info=ip.mm10,min.exon.count=10);


#Cell trajectory modeling by directed diffusion on embedding
#The main parameters are set up by sigma (which limits the range of how far a cell can jump in terms of distance) and n (how many nearest neighbors are being considered for jumps).
#For instance, relaxing (increasing) sigma, in particular will eventually lead to sympathoblast cells “jumping” the gap into the into the chromaffin differentiation part.
pdf("cell_trajectory_modeling_sigma25_n30.pdf")
x <- show.velocity.on.embedding.eu(emb,rvel,n=30,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,nPcs=30,sigma=2.5,show.trajectories=TRUE,diffusion.steps=400,n.trajectory.clusters=15,ntop.trajectories=1,embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=4) 
pdf("cell_trajectory_modeling_sigma25_n30.pdf")

pdf("cell_trajectory_modeling_sigma35_n20.pdf")
x <- show.velocity.on.embedding.eu(emb,rvel,n=20,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,nPcs=20,sigma=3.5,show.trajectories=TRUE,diffusion.steps=400,n.trajectory.clusters=15,ntop.trajectories=1,embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=4) 
dev.off()

pdf("cell_trajectory_modeling_sigma35_n30.pdf")
x <- show.velocity.on.embedding.eu(emb,rvel,n=30,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,nPcs=30,sigma=3.5,show.trajectories=TRUE,diffusion.steps=400,n.trajectory.clusters=15,ntop.trajectories=1,embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=4) 
dev.off()

pdf("cell_trajectory_modeling_sigma35_n10.pdf")
x <- show.velocity.on.embedding.eu(emb,rvel,n=10,scale='sqrt',cell.colors=ac(cell.colors,alpha=cell.alpha),cex=cell.cex,nPcs=30,sigma=3.5,show.trajectories=TRUE,diffusion.steps=400,n.trajectory.clusters=15,ntop.trajectories=1,embedding.knn=T,control.for.neighborhood.density=TRUE,n.cores=4) 
dev.off()



##### scvelo analysis ####
library(reticulate)

use_condaenv("r-velocity", required = TRUE)
scv <- import("scvelo")


scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params()


adata = scv.read('onefilepercell_E8_1_Set_A_A10_S72_and_others_E6ZKY.loom', cache=True)



#### NEW TRIAL #####
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/scvelo.html

library(Seurat)
library(SeuratDisk)
library(SeuratWrappers)

ldat <- ReadVelocity(file = "onefilepercell_E8_1_Set_A_A10_S72_and_others_E6ZKY.loom")

ldat <- lapply(ldat,function(x) {
  colnames(x) <-  gsub(".bam","",gsub(".*:","",colnames(x)))
  x
})

bm <- as.Seurat(x = ldat)
bm[["RNA"]] <- bm[["spliced"]]
bm <- SCTransform(bm)
bm <- RunPCA(bm)
bm <- RunUMAP(bm, dims = 1:20)
bm <- RunTSNE(bm, dims = 1:20)
bm <- FindNeighbors(bm, dims = 1:20)
bm <- FindClusters(bm)
DefaultAssay(bm) <- "RNA"
SaveH5Seurat(bm, filename = "mouse2.h5Seurat")
Convert("mouse2.h5Seurat", dest = "h5ad")


#### cluster 0
heart_small.c0 <- subset(heart, idents = 0)
rownames(heart_small.c0@meta.data)
