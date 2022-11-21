library("Seurat")
library("tidyverse")
time_start <- Sys.time()

options(scipen=200)
args <- commandArgs(T)

seuratObj <- readRDS(args[1])
samples <- unique(seuratObj$sample)

for (samp in samples){
	obj <- subset(seuratObj, subset = sample == samp)
	qc.metrics <- as_tibble(obj[[]], rownames = "Cell.Barcode")
	obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
	p <- qc.metrics %>% arrange(percent.mt) %>% ggplot(aes(nCount_RNA, nFeature_RNA, colour = percent.mt)) +
		geom_point(size = 1) + 
		scale_color_gradientn(colors = c("black", "blue", "green2", "red", "yellow")) +
		scale_x_log10() +
		scale_y_log10() +
		theme_bw() +
		theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
		geom_smooth(method='lm')
	ggsave(paste(samp, "nCount_nFeature.percent.mt.png", sep = "_"), p)
	ggsave(paste(samp, "nCount_nFeature.percent.mt.pdf", sep = "_"), p)
}
exc_time <- difftime(Sys.time(), time_start, units = 'mins')
print(paste0('程序执行时间：', round(exc_time,2), 'mins'))
