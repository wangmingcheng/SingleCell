library(Seurat)
library(tidyverse)

args <- commandArgs(T)
obj <- readRDS(args[1])

#删除不需要的样本
obj <- subset(obj, subset=sample != "SC341")
#为样本分组
obj$group <- ifelse(obj$sample == "SC230" | obj$sample == "SC247" | obj$sample == "SC248", "group1", "group2")

groups <- unique(obj$group)

result <- list()
for (g in groups){
	group_obj <- subset(obj, subset=group==g)
	averExpr <- AverageExpression(group_obj, group.by = "celltype")
	RNA_expr <- as.data.frame(averExpr$RNA)
	RNA_expr$gene <- rownames(RNA_expr)
	RNA_expr <- RNA_expr %>% pivot_longer(cols = colnames(RNA_expr)[-length(colnames(RNA_expr))], names_to = "celltype", values_to = g)
	write.csv(RNA_expr, file = paste(g, "celltype_averExpr.csv", sep = "_"), quote = F, row.names = F)
	result[[g]] <- RNA_expr
}

df <- result[[1]]

for (g in result[[-1]]){
	df <- full_join(df, g, by = c("gene", "celltype"))
	#df <- full_join(df, g)
	#print(colnames(g))
}
write.csv(df, file = paste("group", "celltype_averExpr.csv", sep = "_"), quote = F, row.names = F)
