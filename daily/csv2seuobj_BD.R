library("optparse")

option_list <- list(
	make_option(c("-i", "--indir"),  help="matrix data directory file"),
	make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="."),
	make_option(c("-n", "--name"),   help="seuratObject output name, rds file, default %default", default="result_merged.rds", type = "character"),
	make_option(c("-t", "--type"), help="C4 or 10X, default %default", default="C4", type = "character")

)
opt <- parse_args(OptionParser(option_list=option_list))

#1 把所有样本的单细胞csv格式的表达矩阵合并然后生成1个rds文件
#2 对合并后的对象进行标准化、归一化、降维和聚类处理
#3 对象未进行任何质控处理，如线粒体百分比等
#4 可以统计过滤前的指标用于和质控后比较，也可以以umap图的形式直观查看批次效应
if(FALSE){"
	indir
	├── GSM4430459_MS.sample1.clean.data.txt.gz
	├── GSM4430460_MS.sample2.clean.data.txt.gz
	├── GSM4430461_MS.sample3.clean.data.txt.gz
	├── GSM4430462_MS.sample4.clean.data.txt.gz
	├── GSM4430463_MS.sample5.clean.data.txt.gz
	├── GSM4430464_MS.sample6.clean.data.txt.gz
	├── GSM4430465_MS.sample7.clean.data.txt.gz
	├── GSM4430466_MS.sample8.clean.data.txt.gz
	├── GSM4430467_MS.sample9.clean.data.txt.gz
	├── GSM4430468_MS.sample10.clean.data.txt.gz
	├── GSM4430469_MS.sample11.clean.data.txt.gz
	├── GSM4430470_MS.sample12.clean.data.txt.gz
	├── GSM4430471_MS.sample13.clean.data.txt.gz
	├── GSM4430472_MS.sample14.clean.data.txt.gz
	├── GSM4430473_MS.sample15.clean.data.txt.gz
	├── GSM4430474_MS.sample16.clean.data.txt.gz
	└── GSM4430475_MS.sample17.clean.data.txt.gz
"}

#str substitution method
#stringr package str_replace_all()
#do package Replace()

library("Seurat")
library("stringr")
library("R.utils")

data_dir <- getAbsolutePath(opt$indir)

myfile <- list.files(opt$indir)
#gzfile <- myfile[grep(myfile, pattern =".gz$")]
gzfile <- myfile[grep(myfile, pattern =".csv$")]

Object_list <- list()
for (sample in gzfile) {
	print (sample)
	filedir <- str_c(data_dir, sample, sep = "/")
	#attention the context before first dot regarded as sample name
	#name <- unlist(strsplit(sample, "[.]"))[[1]]
	name <- unlist(strsplit(sample, "_"))[[1]]
	print(name)
	#data <- read.csv(gzfile(filedir), header = T, row.names = 1)
	#BD 平台的单细胞矩阵行是细胞名，列是基因名，所以t()转置
	data <- t(read.csv(filedir, header = T, row.names = 1, comment.char = "#", check.names = F))
	print(dim(data))
	seuratObj <- CreateSeuratObject(counts = data, project = name, min.cells = 0, min.features = 0)
	seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = "^MT-")
	seuratObj[["percent.rb"]] <- PercentageFeatureSet(seuratObj, pattern = "^RP[SL]")
	seuratObj[["percent.hb"]] <- PercentageFeatureSet(seuratObj, pattern = "^HB[^(P)]")
	seuratObj[["sample"]] <- name
	Object_list[[sample]] <- seuratObj
}

merged_seuratObj <- merge(Object_list[[1]], y = Object_list[-1])

merged_seuratObj <- NormalizeData(merged_seuratObj, normalization.method = "LogNormalize", scale.factor = 10000)
merged_seuratObj <- FindVariableFeatures(merged_seuratObj, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(merged_seuratObj)
merged_seuratObj <- ScaleData(merged_seuratObj, features = all.genes)

merged_seuratObj <- RunPCA(merged_seuratObj, verbose = FALSE)
#merged_seuratObj <- RunTSNE(merged_seuratObj, reduction = "pca", dims = 1:30, verbose = FALSE)
merged_seuratObj <- RunUMAP(merged_seuratObj, reduction = "pca", dims = 1:30, verbose = FALSE)
merged_seuratObj <- FindNeighbors(merged_seuratObj, reduction = "pca", dims = 1:30)

saveRDS(merged_seuratObj, paste(opt$outdir, opt$name, sep = "/"))
