library("optparse")

option_list <- list(
        make_option(c("-i", "--inFile"),  help="seurat object file in rds format"),
        make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="./"),
        make_option(c("-t", "--tissue"), help="Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus", type = "character"),
	make_option(c("-g", "--guess"), help="guess tissue type or not, default %default", default=FALSE, type = "logical")
	
)
opt <- parse_args(OptionParser(option_list = option_list))

lapply(c("dplyr", "Seurat", "HGNChelper", "ggplot2"), library, character.only = T)

# load gene set preparation function
source("sc-type/R/gene_sets_prepare.R")
# load cell type annotation function
source("sc-type/R/sctype_score_.R")
source("sc-type/R/auto_detect_tissue_type.R")

seu <- readRDS(opt$inFile)

db_ = "sc-type/ScTypeDB_full.xlsx";

if(isTRUE(opt$guess)){
	tissue_guess <- auto_detect_tissue_type(path_to_db_file = db_, seuratObject = seu, scaled = TRUE, assay = "RNA")  # if saled = TRUE, make sure the data is scaled, as seuratObject[[assay]]@scale.data is used. If you just created a Seurat object, without any scaling and normalization, set scaled = FALSE, seuratObject[[assay]]@counts will be used

	write.table(tissue_guess, file = paste(opt$output, "tissue_guess.csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = F)
}
# DB file
#tissue = "Liver" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list <- gene_sets_prepare(db_, opt$tissue)

# get cell-type by cell matrix
es.max <- sctype_score(scRNAseqData = seu[["RNA"]]@scale.data, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 
# merge by cluster
cL_resutls <- do.call("rbind", lapply(unique(seu@meta.data$seurat_clusters), function(cl){
    es.max.cl <- sort(rowSums(es.max[, rownames(seu@meta.data[seu@meta.data$seurat_clusters == cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seu@meta.data$seurat_clusters == cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

write.table(sctype_scores %>% arrange(cluster) %>% rename(ScType = type), file = paste(opt$output, "sctype.csv", sep = ""), sep = ",", col.names = T, row.names = F, quote = F)

for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seu@meta.data$customclassif[seu@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

p <- DimPlot(seu, reduction = "umap", label = T, repel = TRUE, group.by = 'customclassif') + labs(title = "")
ggsave(paste(opt$output, "dimplot.pdf", sep = ""), p)
ggsave(paste(opt$output, "dimplot.png", sep = ""), p)

if(FALSE){
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)
                                                                                                     
# prepare edges
cL_resutls=cL_resutls[order(cL_resutls$cluster),]; edges = cL_resutls; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; nodes_lvl1$realname = nodes_lvl1$cluster; nodes_lvl1 = as.data.frame(nodes_lvl1); nodes_lvl2 = c(); 
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_resutls$cluster))){
  dt_tmp = cL_resutls[cL_resutls$cluster == unique(cL_resutls$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

mygraph <- graph_from_data_frame(edges, vertices=nodes)

# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")
}
