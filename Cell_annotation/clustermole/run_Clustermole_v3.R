library("optparse")

option_list <- list(
        make_option(c("-i", "--inFile"),  help="seurat object file"),
        make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="./"),
        make_option(c("-s", "--species"), help="hs or mm, hs represent human and mm represent mouse, default %default", default="hs", type = "character"),
        make_option(c("-r", "--runOverlaps"), help="whether or not run Marker gene overlaps step, default %default", default=FALSE, type = "logical")

)
opt <- parse_args(OptionParser(option_list = option_list))

library(Seurat)
library(tidyverse)
library(clustermole)

seu <- readRDS(opt$inFile)

DefaultAssay(seu) <- "RNA"
Idents(seu) <- "seurat_clusters"


write.table(clustermole_markers(), file = paste(opt$outdir, "clustermole_markers.csv", sep = ""), sep = ",", row.names = F, col.names = T, quote = F)

if(isTRUE(opt$runOverlaps)){
	for (cluster in levels(seu$seurat_clusters)){
		markers_df <- FindMarkers(seu, ident.1 = cluster, min.pct = 0.2, only.pos = TRUE, verbose = FALSE)
		markers <- head(rownames(markers_df), 30)
		overlaps_tbl <- clustermole_overlaps(genes = markers, species = opt$species)
		write.table(overlaps_tbl, file = paste(opt$outdir, "cluster", cluster, "_markersOverlap_annotation.csv", sep = ""), sep = ",", row.names = F, col.names = T, quote = F)
	}
}

avg_exp_mat <- AverageExpression(seu, assays = "RNA", slot = "data")[[1]]
avg_exp_mat <- as.matrix(avg_exp_mat)
avg_exp_mat <- log1p(avg_exp_mat)
enrich_tbl <- clustermole_enrichment(expr_mat = avg_exp_mat, species = opt$species)

write.table(enrich_tbl, file = paste(opt$outdir, "markersEnrichment_annotation.csv", sep = ""), sep = ",", row.names = F, col.names = T, quote = F)

spe <- ifelse(opt$species == "hs", "Human", "Mouse")
print(spe)
clusters <- c()
celltypes <- c()
for (clust in unique(enrich_tbl$cluster)){
        d1 <- enrich_tbl %>% filter(species == spe & cluster == clust) %>% slice(1:5) %>% select(c("cluster", "celltype"))
        clusters <- c(clusters, unique(d1$cluster))
        
        celltypes <- c(celltypes, paste(d1$celltype, collapse = "|"))
}

write.table(tibble(cluster = as.numeric(clusters), clustermole = celltypes) %>% arrange(cluster), file = paste(opt$outdir, "clustermole_annotation.csv", sep = ""), sep = "\t", row.names = F, col.names = T, quote = F)
