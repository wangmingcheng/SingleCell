library(Seurat)
library(tidyverse)

args <- commandArgs(T)
obj <- readRDS(args[1])

#删除不需要的样本
obj <- subset(obj, subset=sample != "SC341")
#为样本分组
obj$group <- ifelse(obj$sample == "SC230" | obj$sample == "SC247" | obj$sample == "SC248", "group1", "group2")

p1 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %T>% write_tsv(file = "cell_proportion.tsv") %>% ggplot(aes(sample, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.3, position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "", y = "", title = "")

p2 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %T>% write_tsv(file = "cell_proportion.tsv") %>% ggplot(aes(sample, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.2, position = "dodge") + scale_y_continuous(labels = scales::percent) + labs(x = "", y = "", title = "")


p3 <- tibble(group = obj@meta.data$group, celltype = obj@meta.data$celltype) %>% group_by(group, celltype) %>% count() %T>% write_tsv(file = "group_proportion.tsv") %>% ggplot(aes(group, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.3, position = "fill") + scale_y_continuous(labels = scales::percent) + labs(x = "", y = "", title = "")


p4 <- tibble(group = obj@meta.data$group, celltype = obj@meta.data$celltype) %>% group_by(group, celltype) %>% count() %T>% write_tsv(file = "group_proportion.tsv") %>% ggplot(aes(group, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.3, position = "dodge") + scale_y_continuous(labels = scales::percent) + labs(x = "", y = "", title = "")

#p2 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %T>% write_tsv(file = "cell_proportion.tsv") %>% ggplot(aes(sample, n, fill = celltype)) + geom_bar(stat = 'identity', width = 0.5, position = "dodge") + scale_y_continuous() + labs(x = "", y = "", title = "")

#p3 <- tibble(sample = obj@meta.data$sample, celltype = obj@meta.data$celltype) %>% group_by(sample, celltype) %>% count() %T>% write_tsv(file = "cell_proportion.tsv") %>% ggplot(aes(celltype, n, fill = sample)) + geom_bar(stat = 'identity', width = 0.5, position = "dodge") + scale_y_continuous() + labs(x = "", y = "", title = "")

ggsave("sample_celltype_percent1.pdf", height = 7, width =14, p1)
ggsave("sample_celltype_percent1.png", height = 7, width =14, p1)

ggsave("sample_celltype_percent2.pdf", height = 7, width =14, p2)
ggsave("sample_celltype_percent2.png", height = 7, width =14, p2)

ggsave("group_celltype_percent1.pdf", height = 7, width = 7, p3)
ggsave("group_celltype_percent1.png", height = 7, width = 7, p3)

ggsave("group_celltype_percent2.pdf", height = 7, width = 7, p4)
ggsave("group_celltype_percent2.png", height = 7, width = 7, p4)
