library(Seurat)

args <- commandArgs(T)

seuratObj <- readRDS(args[1])

##全部样本按照同一标准过滤
seuratObj <- subset(seuratObj, subset = nFeature_RNA >= 200 & percent.mt <= 20 & percent.hb <= 20)

saveRDS(seuratObj, file = "result_merged_filtered.rds")
