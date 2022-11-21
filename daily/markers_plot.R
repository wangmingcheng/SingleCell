library("Seurat")
library("stringr")
library("tidyverse")
library("ggplot2")
library("magrittr")

args <- commandArgs(T)
obj <- readRDS(args[1])

#查看细胞类型
cells <- unique(obj$celltype)

top1_markers <- c()
DefaultAssay(obj) <- "RNA"

for (cell in cells){
	print (cell)
	#celltype <- paste0(cell,"_markers.csv",sep="")
	markers <- FindMarkers(obj, ident.1 = cell, group.by = "celltype", min.pct = 0.1)
	#markers <- FindMarkers(obj, group.by = "celltype", min.pct = 0.1)
	write.csv(markers, file = paste0(cell, "_markers.csv", sep = ""), quote = FALSE,)
	
	top4_markers <- markers %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% slice(1:4)
	write.csv(top4_markers, file = paste0(cell, "_top4_markers.csv", sep = ""), quote = FALSE,)
	
	p <- FeaturePlot(obj, features = rownames(top4_markers), cols = c("lightgrey", "red"))
	ggsave(paste0(cell, "_top4_markers", "_FeaturePlot", ".png"), p)
	ggsave(paste0(cell, "_top4_markers", "_FeaturePlot", ".pdf"), p)
	top1_markers <- c(top1_markers, rownames(top4_markers)[1])

	#markers <- intersect(markers, as.data.frame(rownames(obj[['RNA']])))
	#top8_markers <- markers %>% filter(p_val_adj < 0.05) %>% arrange(desc(avg_log2FC)) %>% slice(1:8)
	#p1 <- DoHeatmap(obj, features =  rownames(top8_markers), group.by = "celltype", size = 3 ) + NoLegend()
	#p1 <- DoHeatmap(obj, features =  rownames(top8_markers), group.by = "celltype", size = 3 )
	#ggsave(paste0(cell, "_top8_markers", "_heatmap", ".png"), p1)
	#ggsave(paste0(cell, "_top8_markers", "_heatmap", ".pdf"), p1)
	
	#p2 <- VlnPlot(obj, features = rownames(top8_markers), group.by = "celltype")
	#ggsave(paste0(cell, "_top8_markers", "_vlnplot", ".png"), p2)
	#ggsave(paste0(cell, "_top8_markers", "_vlnplot", ".pdf"), p2)
	
	#write.csv(top8_markers, file = paste0(cell, "_top8_markers.csv", sep = ""), quote = FALSE,)

}

print (top1_markers)

p <- FeaturePlot(obj, features = top1_markers, cols = c("lightgrey", "red"))
ggsave("top1_markers_FeaturePlot.png", p)
ggsave("top1_markers_FeaturePlot.pdf", p)
write.csv(top1_markers, file = "top1_markers.csv")

#p2 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %>% ggplot(aes(sample, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.3, position="fill") + scale_y_continuous(labels = scales::percent)
p2 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %>% mutate(new = if (sample == "Macrophage"){ "CS3" } else if (sample == "XHM_113"){ "CS2" } else if (sample == "YXQ_21.3.25"){ "CS1" } else { "none" }) %T>% write_tsv(file = "cell_proportion.tsv") %>% ggplot(aes(new, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.3, position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "", y = "", title = "")

ggsave("bar.pdf", p2)
ggsave("bar.png", p2)

#p3 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(celltype) %>% count() %>% ggplot(aes(x = "", y = n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.5) + coord_polar(theta = 'y') + labs(x = '', y = '', title = '') + theme(axis.text = element_blank()) + theme(axis.ticks = element_blank()) + theme(legend.position = "none") + geom_text(aes(y = n/2 + c(0, cumsum(n)[-length(n)]), x = sum(n)/150))

#ggsave("pie.pdf", p3)
#ggsave("pie.png", p3)


#查看样本
#unique(obj$sample)

#unique(obj$seurat_clusters)
