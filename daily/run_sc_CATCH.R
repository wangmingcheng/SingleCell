library("optparse")

option_list <- list(
        make_option(c("-i", "--inFile"),  help="seurat object file in rds format"),
        make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="./"),
        make_option(c("-s", "--species"), help="Human or Mouse, default %default", default="Human", type = "character"),
        make_option(c("-t", "--tissue"), help="corresponds human or mouse tissues, more detail visit https://github.com/ZJUFanLab/scCATCH/wiki", type = "character")

)
opt <- parse_args(OptionParser(option_list = option_list))


library(scCATCH)
library(Seurat)
library(dplyr)

seu <- readRDS(opt$inFile)
Idents(seu) <- "seurat_clusters"

obj <- createscCATCH(data = seu[['RNA']]@data, cluster = as.character(Idents(seu)))
#obj <- findmarkergene(object = obj, species = "Human", marker = cellmatch, tissue = "Small intestine")
obj <- findmarkergene(object = obj, species = opt$species, if_use_custom_marker  = TRUE, marker = cellmatch, tissue = opt$tissue)
obj <- findcelltype(object = obj)

#anno1 <- data.frame(cluster = seu@meta.data$seurat_clusters, celltype = seu@meta.data$celltype)  %>% distinct(cluster, .keep_all = TRUE) %>% arrange(cluster)
anno2 <- data.frame(cluster = obj@celltype$cluster, scCATCH = obj@celltype$cell_type)

#compare_result <- inner_join(anno1, anno2, by = "cluster")

write.table(anno2, file = paste(opt$outdir, "scCATCH_anno.csv", sep = ""), sep = ",", quote = F, col.names = T, row.names = F)

if(FALSE){
	metadata <- as.data.frame(obj@celltype$cell_type, obj@celltype$cluster)
	metadataDf <- metadata[match(seu$seurat_clusters, rownames(metadata)), drop=FALSE]
	rownames(metadataDf) <- colnames(seu)
	seuratObj <- AddMetaData(object=seu, metadata=metadataDf, col.name=colnames(metadata))
	seuratObj <- AddMetaData(object = seu, metadata = obj@m$cell_type, col.name = colnames(metadata))
	write.table(obj@meta[,c("cell", "cluster")], file = "scCATCH_singlecell_anno.csv", sep = ",", quote = F, col.names = F, row.names = F)
	saveRDS(seuratObj, "liver.rds")
}
