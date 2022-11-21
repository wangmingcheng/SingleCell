library("Seurat")
library("stringr")
library("tidyverse")
library("ggplot2")
library("magrittr")

args <- commandArgs(T)
obj <- readRDS(args[1])

#查看细胞类型
cells <- unique(obj$celltype)

DefaultAssay(obj) <- "RNA"
#Idents(obj) <- "celltype"
#Idents(obj) <- "sample"

#p1 <- VlnPlot(obj, features = c("PLAUR", "GRPR", "S100A9"))
#ggsave("VlnPlot_gene.pdf", p1)

#p2 <- FeaturePlot(obj, features = c("PLAUR", "GRPR", "S100A9"))

#ggsave("Features_gene.pdf", p2)
#f <- function(string){str_split(string, "_")[[1]][2]}
#obj$group <- unlist(lapply(obj$sample, f))

#Idents(obj) <- "group"
markers <- FindAllMarkers(obj, min.pct = 0.1)
#H <- subset(obj, subset = group == "H")
#HL <- subset(obj, subset = group == "NL")
#markers <- FindMarkers(obj, assay = "RNA", ident.1 = "NL", ident.2 = "LS", group.by = "group", min.pct = 0.1)
#markers <- FindMarkers(obj, cells.1 = Cells(H), cells.2 = Cells(NL), min.pct = 0.1)
#markers <- FindMarkers(obj, group.by = "celltype", min.pct = 0.1)
write.csv(markers, file = "all_NL-LS_markers.csv", quote = FALSE,)

