library(Seurat)
library(ggplot2)
library(tidyverse)
library(future)
#plan("multisession", workers = 3)
#options(future.globals.maxSize = 100 * 1024^3)

args <- commandArgs(T)
seuratObj <- readRDS(args[1])

#ifnb.list <- list()
Object_list <- list()
for (samp  in unique(seuratObj$sample)){
	seu <- subset(seuratObj, subset = sample == samp)
	seu <- SCTransform(seu, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
	Object_list[[samp]] <- seu
}

print(Object_list)
if(0){'
	seuratObj <- SCTransform(seuratObj, vst.flavor = "v2", verbose = FALSE) %>%
	RunPCA(npcs = 30, verbose = FALSE) %>%
    	RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    	FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
    	FindClusters(resolution = 0.7, verbose = FALSE)

	p <- DimPlot(seuratObj, label = T, repel = T) + ggtitle("Unsupervised clustering")
	p1 <- DimPlot(seuratObj, label = T, repel = T, group.by = "sample") + ggtitle("Unsupervised clustering")

	ggsave("dimplot.pdf", p)
	ggsave("dimplot_sample.pdf", p1)
'}
#seuratObj <- SCTransform(seuratObj, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 30, verbose = FALSE)
ifnb.list <- Object_list
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

combined.sct <- RunPCA(combined.sct, verbose = FALSE)
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
combined.sct <- FindNeighbors(combined.sct, reduction = "pca", dims = 1:30)
#combined.sct <- FindClusters(combined.sct, resolution = 0.3)
saveRDS(combined.sct, file = "intergrated.rds")
