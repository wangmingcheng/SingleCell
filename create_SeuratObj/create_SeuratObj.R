library("optparse")

option_list <- list(
	make_option(c("-i", "--indir"),  help="C4 or 10X data directory file"),
	make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="."),
	make_option(c("-n", "--name"),   help="seuratObject output name, rds file, default %default", default="result_merged.rds", type = "character"),
	make_option(c("-t", "--type"), help="C4 or 10X, default %default", default="C4", type = "character")

)
opt <- parse_args(OptionParser(option_list=option_list))

if(FALSE){"
	indir
	├── C
	│   ├── barcodes.tsv.gz
	│   ├── features.tsv.gz
	│   └── matrix.mtx.gz
	├── D
	│   ├── barcodes.tsv.gz
	│   ├── features.tsv.gz
	│   └── matrix.mtx.gz
	└── E
	    ├── barcodes.tsv.gz
	    ├── features.tsv.gz
	    └── matrix.mtx.gz
"}

#str substitution method
#stringr package str_replace_all()
#do package Replace()

library("Seurat")
library("stringr")
library("R.utils")

data_dir <- getAbsolutePath(opt$indir)

Object_list <- list()
index <- 1
for (sample in list.files(opt$indir)) {
	filedir <- str_c(data_dir, sample, sep = "/")
	print (filedir)
	if (opt$type == "C4"){
		#barcodes名字中的字符"_"替换成字符":",之后在名字末尾加"-1"(第1个样本加"-1",第2个样本加"-2",以此类推)
		barcodes <- read.table(paste(filedir, "barcodes.tsv.gz", sep = "/"), header = FALSE)
		for (i in 1:nrow(barcodes)) {
                	barcodes[i,1] = paste(gsub("_", ":", barcodes[i,1]), index, sep = "-")
		}
		write.table(barcodes, file=gzfile(paste(filedir, "barcodes.tsv.gz", sep = "/")), quote = FALSE, row.names = FALSE, col.names = FALSE)
	
		#features中字符"_"替换成"-"
		features <- read.table(paste(filedir, "features.tsv.gz", sep = "/"), header = FALSE)
		for (j in 1:nrow(features)) {
                	barcodes[j,1] = gsub("_", "-", features[j,1])
        	}
		write.table(barcodes, file=gzfile(paste(filedir, "features.tsv.gz", sep = "/")), quote = FALSE, row.names = FALSE, col.names = FALSE)
		index <- index + 1
		counts <- Read10X(filedir, gene.column = 1)
	}else if(opt$type == "10X"){
		counts <- Read10X(filedir, gene.column = 2)
	}else{
		stop("You shall choose C4 or 10X")
	}

	seuratObj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3, min.features = 0)
	seuratObj[["sample"]] <- sample
	Object_list[[sample]] <- seuratObj
}

merged_seuratObj <- merge(Object_list[[1]], y = Object_list[-1])
saveRDS(merged_seuratObj, paste(opt$outdir, opt$name, sep = "/"))
