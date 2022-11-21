library(Seurat)
library(ggplot2)
library(scCustomize)

args <- commandArgs(T)

seuratObj <- readRDS(args[1])

feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb", "percent.hb")
#VlnPlot(alldata, group.by = "orig.ident", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()
p <- VlnPlot(seuratObj, group.by = "sample", features = feats, pt.size = 0.1, ncol = 3) + NoLegend()

ggsave("beforeFilter_Vlnplot.pdf", width = 14, height = 7, p)
ggsave("beforeFilter_Vlnplot.png", width = 14, height = 7, p)


p1 <- FeatureScatter(seuratObj, "nCount_RNA", "nFeature_RNA", group.by = "sample", pt.size = 0.1)
ggsave("FeatureScatter.pdf", p1)

#p1 <- DimPlot(seuratObj, label = T, repel = T, group.by = "sample", split.by = "sample") + ggtitle("sample")
#p1 <- DimPlot(seuratObj, label = T, repel = T, cells = Cells(subset(seuratObj, subset=sample=="SC230"))) + ggtitle("sample")
#ggsave("SC230_dimplot_sample.pdf", p1)

samples <- unique(seuratObj$sample)
for (sample in samples){
	p2 <- Meta_Highlight_Plot(seurat_object = seuratObj,
				  meta_data_column = "sample",
				  meta_data_highlight = sample,
				  highlight_color = "firebrick",
				  background_color = "lightgray",
				  pt.size = 0.1
	)
	#ggsave("SC230_dimplot_sample.pdf", p2)
	ggsave(paste(sample,"highlight_dimplot.pdf", sep = "_"), p2)
	ggsave(paste(sample,"highlight_dimplot.png", sep = "_"), p2)

}



#pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 00 & percent.mt < 5)
