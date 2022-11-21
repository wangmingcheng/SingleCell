library(Seurat)
library(ggplot2)

args <- commandArgs(T)
seuratObj <- readRDS(args[1])


#seuratObj <- RunPCA(seuratObj)
#seuratObj <- RunUMAP(seuratObj, reduction = "pca", dims = 1:30)

## if seuratObj$celltype
p <- DimPlot(seuratObj, label = T, repel = T) + ggtitle("cluster")
p1 <- DimPlot(seuratObj, label = T, repel = T, group.by = "sample") + ggtitle("sample")
p2 <- DimPlot(seuratObj, label = T, repel = T, group.by = "celltype") + ggtitle("celltype")

ggsave("dimplot_cluster.pdf", p)
ggsave("dimplot_cluster.png", p)
ggsave("dimplot_sample.pdf", p1)
ggsave("dimplot_sample.png", p1)
ggsave("dimplot_celltype.pdf", p2)
ggsave("dimplot_celltype.png", p2)
