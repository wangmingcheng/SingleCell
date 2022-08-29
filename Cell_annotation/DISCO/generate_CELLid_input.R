library("optparse")

option_list <- list(
        make_option(c("-i", "--inFile"),  help="seurat object file in rds format"),
        make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="./"),
        make_option(c("-a", "--anno"), help="DISCO annotation file, default name is spreadsheet.csv", type = "character"),
        make_option(c("-t", "--tissue"), help="corresponds Human or Mouse tissues, more detail visit https://github.com/ZJUFanLab/scCATCH/wiki", type = "character")

)
opt <- parse_args(OptionParser(option_list = option_list))

library(Seurat)
library(ggplot2)

#if(!is.null(opt$inFile)){
	seu <- readRDS(opt$inFile)

	Idents(seu) <- "seurat_clusters"

	cluster_average = AverageExpression(seu)
	# This will generate average expression for each cluster
	cluster_average = round(cluster_average$RNA, 2)
	write.table(cluster_average, file = paste(opt$outdir, "CELLiD_input.txt", sep = ""), quote = F, col.names = F, row.names = T, sep="\t")
	# Then, you can upload this file to https://www.immunesinglecell.org/cellpredictor predict cell types
#}

# After you get the results, you can add the predicted cell type to seurat object as follows:
if(!is.null(opt$anno)){
	#predicted.ct = read.csv("spreadsheet.csv")
	predicted.ct = read.csv(opt$anno)
	seu$primary.predict = predicted.ct[as.numeric(seu$seurat_clusters), 2]
	seu$secondary.predict = predicted.ct[as.numeric(seu$seurat_clusters), 3]
	p <- DimPlot(seu, group.by = "primary.predict", label = T) + labs(title = "")

	ggsave(paste(opt$outdir, "Disco.pdf", sep = ""), p)
	ggsave(paste(opt$outdir, "Disco.png", sep = ""), p)
}
