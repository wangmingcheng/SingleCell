library("Seurat")
library("ggplot2")

args <- commandArgs(T)
#sample <- c("Normal54", "Normal55", "Injury50", "Injury51")
#data <- readRDS("result.no_filter.rds")
data <- readRDS(args[1])
samples <- unique(data$sample)

save_plot <- function(result, fun, features=6){
	pdf(paste(result, ".pdf", sep = ""), width = features * 3, height = 6)
	a <- dev.cur()
	png(paste(result, ".png", sep = ""), width = features * 3, height = 6, unit = "in", res = 300)
	dev.control("enable")
	print(fun)
	dev.copy(which = a)
	dev.off()
	dev.off()
}

for (i in samples){
	s <- subset(data, sample == i)
	Idents(s) <- i
	#ii <- trimws(i, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")	
	#feats <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
	feats <- c("nFeature_RNA", "nCount_RNA")
	#p <- VlnPlot(s, features = feats, ncol = 3)	
	p <- VlnPlot(s, features = feats, ncol = 2)	
	save_plot(i, p, length(feats))
}
