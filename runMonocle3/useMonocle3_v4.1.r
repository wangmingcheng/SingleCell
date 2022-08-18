library("optparse")

option_list <- list(
    make_option(c("-i", "--inFile"), help="seurat object rds file"),
    make_option(c("-o", "--outPrefix"), help="output prefix, default %default", default="result"),
    make_option(c("-r", "--recluster"), help="whether re-process data and reduce dimensionality use monolce3 pipline, default %default", default=FALSE, type="logical"),
    make_option(c("-b", "--batch"), help="Remove batch effects, estimate the level of background contamination in each batch of cells and subtract it, default is percent.mt and percent.rp; example percent.hb,Phase", type="character"),
    make_option(c("--use_partition"), help="whether plot cell trajactory in different partitions, default %default", default=TRUE, type="logical"),
    make_option(c("--root_type"), help="specifiy a celltype as root"),
    make_option(c("--root_vertex"), help="specifiy a closest_vertex from groupBy_closest_vertex.png as root node, if --root_type is offered, --root_vertex will be ignored", type="integer"),
    
    make_option(c("--min_expression"), help="the minimum expression level to plot in function plot_genes_in_pseudotime., default %default", default=0.5, type="double"),
    make_option(c("--width"), help="image width, default %default", default=8, type="double"),
    make_option(c("--height"), help="image width, default %default", default=6, type="double"),
    make_option(c("--cores"), help="cpus running graph_test(), default %default", default=12, type="integer")
)
opt <- parse_args(OptionParser(option_list=option_list))

##################################################
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(patchwork)

#library(ggthemes)
#library(ggbeeswarm)
#library(tidyverse)

options(ggrepel.max.overlaps = 10000)

plotQuasirandom <- function(data, groupBy){
    df <- data[,  c("pseudotime", groupBy)]
    colnames(df) <- c("pseudotime", "groupBy")

    p <- ggplot(df, aes(x = pseudotime, y = groupBy, colour = groupBy)) +
        geom_quasirandom(groupOnX = FALSE) +
        scale_color_tableau() + 
        theme_classic() +
        xlab("pseudotime") + 
        ylab(groupBy) +
        guides(colour = "none")

    ggsave(paste(opt$outPrefix, "monocle3.quasirandom", paste0("groupBy_", groupBy), "png", sep="."), width=opt$width, height=opt$height, p)
    ggsave(paste(opt$outPrefix, "monocle3.quasirandom", paste0("groupBy_", groupBy), "pdf", sep="."), width=opt$width, height=opt$height, p)
}

seuratObj <- readRDS(opt$inFile)
DefaultAssay(seuratObj) <- "RNA"
celltypes <- unique(seuratObj$celltype)

root <- names(sort(table(seuratObj[["celltype"]]), decreasing = TRUE)[1])
print (paste(root, "as root for analyse", sep = " "))
if (is.null(opt$root_type)){
	opt$root_type <- root
}
create_cds_obj <- function(input, recluster, batch=NULL){
	if(isTRUE(recluster)){
		expression_matrix <- GetAssayData(input, assay = 'RNA', slot = 'counts')
		cell_metadata <- input@meta.data
		gene_annotation <- data.frame(gene_short_name = rownames(seuratObj))
		rownames(gene_annotation) <- rownames(seuratObj)
		cds <- new_cell_data_set(expression_matrix,
			cell_metadata = cell_metadata,
			#gene_metadata = gene_annotation[row.names(expression_matrix), ]
			gene_metadata = gene_annotation
		)
		cds <- preprocess_cds(cds, num_dim = 50)
		p <- plot_pc_variance_explained(cds)
		ggsave(paste(opt$outPrefix, "elbowPlot.pdf", sep = "_"), p)
		ggsave(paste(opt$outPrefix, "elbowPlot.png", sep = "_"), p)
		
		batch_str <- paste("~ percent.mt ", "percent.rp", sep = " + ")
		if(!is.null(batch)){
			batch <- gsub(",", " + ", opt$batch)
			batch_str <- paste("~ percent.mt ", "percent.rp", batch, sep = " + ")
		}
		cds <- align_cds(cds, residual_model_formula_str = batch_str)
		cds <- reduce_dimension(cds, preprocess_method = "PCA")
		p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('cds.umap')
		cds.embed <- cds@int_colData$reducedDims$UMAP
		seurat.embed <- Embeddings(seuratObj, reduction = "umap")
		seurat.embed <- seurat.embed[rownames(cds.embed),]
		cds@int_colData$reducedDims$UMAP <- seurat.embed
		p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('seurat.umap')
		p = p1|p2
		ggsave(paste(opt$outPrefix,"Reduction_Compare.pdf", sep = "_"), plot = p, width = 10, height = 5)
		ggsave(paste(opt$outPrefix,"Reduction_Compare.png", sep = "_"), plot = p, width = 10, height = 5)
	}else{
		#DefaultAssay(seuratObj) <- "RNA"
		input@active.assay <- "RNA"
		cds <- as.cell_data_set(input, assay = "RNA")
		cds <- estimate_size_factors(cds)
	}
	cds
}

cds <- create_cds_obj(seuratObj, opt$recluster, opt$batch)
#save_monocle_objects(cds = cds, directory_path = paste(opt$outPrefix, "monocle3_cds_objects", sep = "_"), comment = 'For test')
#q()

cds <- cluster_cells(cds = cds)
cds <- learn_graph(cds, use_partition = opt$use_partition)

get_earliest_principal_node <- function(cds, cell_phenotype, root_type, root_vertex){
    root_pr_nodes <- ""
    if(!is.null(root_vertex)){
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(root_vertex)]
	}else{
        cell_ids <- which(colData(cds)[, cell_phenotype] == root_type)
        closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
        closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
        root_pr_nodes <-  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    }
    root_pr_nodes
}

colData(cds)$closest_vertex <- as.factor(cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex)

p <- plot_cells(
    cds,
    color_cells_by = "closest_vertex",
    label_groups_by_cluster = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    cell_size = 1,
    group_label_size = 4,
    #graph_label_size = 3
)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_closest_vertex.png", sep = "."), p)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_closest_vertex.pdf", sep = "."), p)

root_pr_nodes = get_earliest_principal_node(cds, cell_phenotype = "celltype", root_type = opt$root_type, root_vertex = opt$root_vertex)
cds <- order_cells(cds, root_pr_nodes = root_pr_nodes)

p <- plot_cells(
	cds,
	color_cells_by = "pseudotime",
	#label_groups_by_cluster = FALSE,
	#label_leaves = TRUE,
	#label_branch_points = TRUE,
	cell_size = 1,
	graph_label_size = 3
)

ggsave(paste(opt$outPrefix, "monocle3.groupBy_pseudotime.png", sep = "."), width = opt$width, height = opt$height, p)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_pseudotime.pdf", sep = "."), width = opt$width, height = opt$height, p)

df <- cbind(data.frame(barcode = colnames(cds)), colData(cds)[, c("sample", "celltype")], data.frame(pseudotime = pseudotime(cds)))
write.table(df, paste(opt$outPrefix, "monocle3.pseudotime.txt", sep = "."), sep = "\t", row.names = F, quote = F)

#plotQuasirandom(df, groupBy="celltype")
#plotQuasirandom(df, groupBy="sample")

p <- plot_cells(
	cds,
	color_cells_by = "celltype",
	label_groups_by_cluster = FALSE,
	label_leaves = TRUE,
	label_branch_points = TRUE,
	cell_size = 1,
	group_label_size = 4,
	graph_label_size = 3
)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_celltype.png", sep = "."), p)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_celltype.pdf", sep = "."), p)

p <- plot_cells(
	cds,
	color_cells_by = "sample",
	label_groups_by_cluster = FALSE,
	label_leaves = TRUE,
	label_branch_points = TRUE,
	cell_size = 1,
	label_cell_groups = FALSE,
	graph_label_size = 3
)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_sample.png", sep = "."), p)
ggsave(paste(opt$outPrefix, "monocle3.groupBy_sample.pdf", sep = "."), p)

## Calculate size factors using built-in function in monocle3
#cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(seuratObj[["RNA"]])

#Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>% pull(gene_short_name) %>% as.character()
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = opt$cores)
pr_graph_test_res_list <- subset(ciliated_cds_pr_test_res, q_value < 0.05) %>% arrange(desc(morans_I))
Pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
write.table(pr_graph_test_res_list, file = paste(opt$outPrefix, "pr_graph_test_res.xls", sep = "_"), sep = "\t", row.names = F, quote = F)

genelist <- row.names(pr_graph_test_res_list %>% slice(1:4))

if(isTRUE(opt$recluster)){
	gene_module_df <- find_gene_modules(cds[Pr_deg_ids, ], resolution = c(10^seq(-6,-1)))
	write.table(gene_module_df, file = paste(opt$outPrefix, "gene_module_df.xls", sep = "_"), sep = "\t", row.names = F, quote = F)
	cell_group_df <- tibble(cell=row.names(colData(cds)), cell_group = colData(cds)$celltype)
	agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
	row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
	p <- pheatmap(agg_mat, scale="column", clustering_method = "ward.D2")
	ggsave(paste(opt$outPrefix, "celltype_module.heatmap.pdf", sep = "_"), p)
	ggsave(paste(opt$outPrefix, "celltype_module.heatmap.png", sep = "_"), p)
}

lineage_cds <- cds[rowData(cds)$gene_short_name %in% genelist, colData(cds)$celltype %in% celltypes]
print(lineage_cds)
p <- plot_genes_in_pseudotime(lineage_cds, color_cells_by = "celltype", min_expr = opt$min_expression)

ggsave(paste(opt$outPrefix, "top4_genes_in_celltype.pdf", sep = "_"), p)
ggsave(paste(opt$outPrefix, "top4_genes_in_celltype.png", sep = "_"), p)

p <- plot_genes_in_pseudotime(lineage_cds, color_cells_by = "pseudotime", min_expr = opt$min_expression)
ggsave(paste(opt$outPrefix, "top4_genes_in_pseudotime.pdf", sep = "_"), p)
ggsave(paste(opt$outPrefix, "top4_genes_in_pseudotime.png", sep = "_"), p)

p <- plot_cells(cds,
	genes = genelist,
	label_cell_groups = FALSE,
	show_trajectory_graph = FALSE,
)
ggsave(paste(opt$outPrefix, "top4_featureplot.pdf", sep = "_"), p)
ggsave(paste(opt$outPrefix, "top4_featureplot.png", sep = "_"), p)

for (celltype in celltypes){
	print(paste("Dealing with", celltype, sep = " "))
	celltype_cds <- paste0(celltype, "_cds", sep = "")
	celltype_cds <- cds[, grepl(celltype, colData(cds)$celltype, ignore.case=TRUE)]
	pr_graph_test_res <- graph_test(celltype_cds, neighbor_graph="knn", cores = opt$cores)
	pr_graph_test_res_list <- subset(pr_graph_test_res, q_value < 0.05) %>% arrange(desc(morans_I))
	write.table(pr_graph_test_res_list, file = paste(opt$outPrefix, celltype, "pr_test_res.xls", sep = "_"), sep = "\t", row.names = F, quote = F)
	
	genelist <- row.names(pr_graph_test_res_list %>% slice(1:4))
	
	p <- plot_cells(cds,
        	genes = genelist,
        	label_cell_groups = FALSE,
        	show_trajectory_graph = FALSE,
	)
	ggsave(paste(opt$outPrefix, celltype, "top4_featureplot.pdf", sep = "_"), p)
	ggsave(paste(opt$outPrefix, celltype, "top4_featureplot.png", sep = "_"), p)
	
	lineage_cds <- cds[rowData(cds)$gene_short_name %in% genelist,
                       colData(cds)$celltype %in% celltypes]	
	
	p <- plot_genes_in_pseudotime(lineage_cds,
        	color_cells_by = "celltype",
		min_expr = opt$min_expression,
	)
	ggsave(paste(opt$outPrefix, celltype, "top4_genes_in_pseudotime.pdf", sep = "_"), p)
	ggsave(paste(opt$outPrefix, celltype, "top4_genes_in_pseudotime.png", sep = "_"), p)
	print("Done")
}

save_monocle_objects(cds = cds, directory_path = paste(opt$outPrefix, "monocle3_cds_objects", sep = "_"), comment = 'Monocle3 single cell trajectories analysis')
