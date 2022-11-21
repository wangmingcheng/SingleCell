library("optparse")

option_list <- list(
	make_option(c("-i", "--indir"),  help="matrix data directory file"),
	make_option(c("-o", "--outdir"), help="seuratObject output directory, default %default", default="."),
	make_option(c("-n", "--name"),   help="seuratObject output name, rds file, default %default", default="result_merged.rds", type = "character"),
	make_option(c("-t", "--type"), help="C4 or 10X, default %default", default="C4", type = "character")

)
opt <- parse_args(OptionParser(option_list=option_list))

#Rownames are gene, colnames are cell
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
gzfile <- myfile[grep(myfile, pattern =".gz$")]

Object_list <- list()
for (sample in gzfile) {
	filedir <- str_c(data_dir, sample, sep = "/")
	#attention the context before first dot regarded as sample name
	name <- unlist(strsplit(sample, "[.]"))[[1]]
	print(name)
	data <- read.csv(gzfile(filedir), header = T, row.names = 1)
	seuratObj <- CreateSeuratObject(counts = data, project = name, min.cells = 3, min.features = 0)
	seuratObj[["sample"]] <- name
	Object_list[[sample]] <- seuratObj
}

merged_seuratObj <- merge(Object_list[[1]], y = Object_list[-1])
saveRDS(merged_seuratObj, paste(opt$outdir, opt$name, sep = "/"))
