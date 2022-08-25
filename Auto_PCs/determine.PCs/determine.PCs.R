if(FALSE){
	"主成分累积贡献大于95%;
	PC本身对方差贡献小于3%;
	两个连续PCs之间差异小于0.1%;"
}
library(ggplot2)

args <- commandArgs(T)

obj <- readRDS(args[1])
pct <- obj [["pca"]]@stdev / sum( obj [["pca"]]@stdev ) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 95 & pct < 3)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs

plot_df <- data.frame(pct = pct,   cumu = cumu,   rank = 1:length(pct))


# Elbow plot to visualize 
p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + 
  geom_vline(xintercept = 90, color = "grey") + 
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
ggsave("Elbow_plot.pdf", p)
ggsave("Elbow_plot.png", p)
