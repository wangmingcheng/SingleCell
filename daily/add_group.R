# methods1
## f <- function(string){str_split(string, "_")[[1]][2]}
## obj$group <- unlist(lapply(obj$sample, f))

library(Seurat)
args <- commandArgs(T)
obj <- readRDS(args[1])
#删除不需要的样本
obj <- subset(obj, subset=sample != "SC341")
#为样本分组
obj$group <- ifelse(obj$sample == "SC230" | obj$sample == "SC247" | obj$sample == "SC248", "group1", "group2")

markers <- FindMarkers(obj, ident.1 = "group1", group.by = "group", min.pct = 0.1)

write.csv(markers, file = "group1_VS_group2_DE.csv", quote = FALSE,)

Idents(obj) <- "celltype"

celltypes <- unique(obj$celltype)

for (celltype in celltypes){
	#celltype_obj <- subset(obj, subset=celltype==celltype)
	markers <- FindMarkers(obj, ident.1 = "group1", group.by = 'groups', subset.ident = celltype)
	write.csv(markers, file = paste(celltype, "group1_VS_group2_DE.csv", sep = "_"), quote = FALSE,)
}
