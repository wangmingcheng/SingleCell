mat <- as.matrix(mat@assays$RNA@counts)

gene_attr <- data.frame(mean = rowMeans(mat), 
                        detection_rate = rowMeans(mat > 0), 
                        var = apply(mat, 1, var))


gene_attr$log_mean <- log10(gene_attr$mean)
gene_attr$log_var <- log10(gene_attr$var)
rownames(gene_attr) <- rownames(mat)


p <- ggplot(gene_attr, aes(log_mean, log_var)) +
	geom_point(alpha = 0.3, shape = 16) +
	geom_density_2d(size = 0.3) +
	geom_abline(intercept = 0, slope = 1, color = "red")

ggsave("log_mean-log_var_relation.pdf", p)
ggsave("log_mean-log_var_relation.png", p)

x = seq(from = -3, to = 2, length.out = 1000)
poisson_model <- data.frame(log_mean = x, detection_rate = 1 - dpois(0, lambda = 10^x))
p1 <- ggplot(gene_attr, aes(log_mean, detection_rate)) +
	geom_point(alpha = 0.3, shape = 16) +
	geom_line(data = poisson_model, color = "red") +
	theme_gray(base_size = 8)
ggsave("log_mean-detection_rate_relation.pdf", p1)
ggsave("log_mean-detection_rate_relation.png", p1)

cell_attr <- data.frame(n_umi = colSums(mat),
                        n_gene = colSums(mat > 0))
p2 <- ggplot(cell_attr, aes(n_umi, n_gene)) + geom_point(alpha = 0.3, shape = 16) + geom_density_2d(size = 0.3)
ggsave("n_umi-n_gene_relation.pdf", p2)
ggsave("n_umi-n_gene_relation.png", p2)
