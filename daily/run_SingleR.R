library("optparse")

option_list <- list(
        make_option(c("-i", "--inFile"), help="QC filtered seurat object"),
        make_option(c("-o", "--outPrefix"), help="output prefix, default %default", default="result"),
        make_option(c("-r", "--refGetter"), help="reference database getter function name, like HumanPrimaryCellAtlasData, ref_Mouse_all, ref_Mouse_imm etc."),
        make_option(c( "--refParam"), help="parameter in reference database getter function"),
        make_option(c( "--refLabel"), help="label attribute in reference data used for SingleR prediction, default %default", default="label.main"),
        make_option(c( "--refFilterBy"), help="which label attribute used to filter reference data involved in SingleR prediction, default %default", default="label.main"),
        make_option(c( "--refFilterLabels"), help="which labels in reference data will be used in SingleR prediction"),
        make_option(c("--deMethod"), help="classic or wilcox or t, default %default", default="wilcox"),
        make_option(c("--deN"), help="top differential expression gene number to use for exploring cell type, default %default", default=10),
        make_option(c("--assay"), help="assay used for seurat object of test data, default %default", default="RNA"),
        make_option(c("--assayTypeTest"), help="assay data type used for seurat object of test data, counts or logcounts, default %default", default="counts"),
        make_option(c("--clusters"), help="whether or not annotation by clusters, default %default", default=TRUE, type="logical")
)
opt <- parse_args(OptionParser(option_list=option_list))

##################################################
library(Seurat)
library(tidyverse)
library(SingleR)
library(scRNAseq)
library(scater)
library(BiocParallel)

#ref <- get(load("/rainbow/jianj/bin/singCell/useSingleR/reference/HumanPrimaryCellAtlasData.RData"))
ref <- celldex::HumanPrimaryCellAtlasData()

# load seurat object
seuratObj <- readRDS(opt$inFile)
DefaultAssay(seuratObj) <- opt$assay

cluster <- ifelse(isTRUE(opt$clusters), as.vector(levels(seuratObj$seurat_clusters)), NULL)
test <- as.SingleCellExperiment(seuratObj, assay=opt$assay)

labels <- NULL
if(!is.null(ref[[opt$refLabel]])){
        labels <- ref[[opt$refLabel]]
}else{
	message(paste0("[ERROR] Cant'find ", opt$refLabel," attribute in reference data! Available colnames of coldata of reference data are ", paste(colnames(colData(ref)), collapse=",") ))
	quit()
}
cluster_dir <- dir.create("singleR/cluster_anno", recursive = TRUE)
clusterprefix <- paste("singleR/cluster_anno", opt$outPrefix, sep = "/")
print(clusterprefix)
Idents(seuratObj) <- "seurat_clusters"
result_cluster <- SingleR(
        test = test, 
        ref = ref, 
        labels = ref$label.main, 
        clusters= Idents(seuratObj),
        de.method=opt$deMethod,
        de.n=opt$deN,
        assay.type.test = opt$assayTypeTest,
	BPPARAM = MulticoreParam(8)
)

seuratObj$SingleR_cluster <- result_cluster[seuratObj$seurat_clusters, "pruned.labels"]

stat <- as.data.frame(table(seuratObj$SingleR_cluster))
colnames(stat) <- c("cellType", "number")
write.table(stat, paste(clusterprefix, "SingleR_cluster_stat.txt", sep = "_"), sep="\t", quote=F, row.names=F)

write.csv(data.frame(cluster = as.numeric(unique(rownames(result_cluster))), SingleR = result_cluster$pruned.labels) %>% arrange(cluster), file = "singleR_cluster_anno.txt", row.names = F, quote = F)
result_cluster <- cbind(data.frame(item=rownames(result_cluster)), result_cluster)
write.table(result_cluster, paste(clusterprefix, "SingleR_cluster.txt", sep = ""), sep="\t", quote=F, row.names=F)
saveRDS(seuratObj, paste(clusterprefix, "SingleR_cluster.rds", sep = ""))

if(0){
result_cell <- SingleR(
	test = test, 
	ref = ref, 
	labels = labels, 
        de.n = opt$deN,
        assay.type.test = opt$assayTypeTest
)

seuratObj$SingleR_cell <- result[seuratObj$seurat_clusters, "pruned.labels"]

stat <- as.data.frame(table(seuratObj$SingleR_cluster))
colnames(stat) <- c("cellType", "number")
write.table(stat, paste(opt$outPrefix, "SingleR_cluster", "stat.txt", sep = "."), sep = "\t", quote = F, row.names = F)

write.csv(data.frame(cluster = as.numeric(unique(rownames(result_cluster))), SingleR = result_cluster$pruned.labels) %>% arrange(cluster), file = "singleR_cluster_anno.txt", row.names = F, quote = F)
result <- cbind(data.frame(item=rownames(result)), result)
write.table(result_cell, paste(opt$outPrefix, "SingleR", opt$method, "txt", sep="."), sep="\t", quote=F, row.names=F)
saveRDS(seuratObj, paste(opt$outPrefix, "SingleR", opt$method, "rds", sep="."))
}
