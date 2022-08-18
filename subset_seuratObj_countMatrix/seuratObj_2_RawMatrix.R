library("Seurat")

args <- commandArgs(T)
seuratObj <- readRDS(args[1])

samples <- unique(seuratObj$sample)

for (samp in samples){
	obj <- subset(seuratObj, subset = sample == samp)
	write.table(t(as.matrix(GetAssayData(object = obj, slot = "counts"))), paste(samp, 'counts.csv', sep = "_"), sep = ',', row.names = T, col.names = T, quote = F)
}
