library("Seurat")
data <- read.csv("GSM4430459_MS.sample1.clean.data.txt", headeri = T, row.names = 1)
data <- as.matrix(data)
seu <- CreateSeuratObject(counts = data, project = "test", min.cells = 3, min.features = 0)
