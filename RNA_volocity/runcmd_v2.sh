Rscript Seuratobj_deal.R result.renamed.rds
生成4个文件counts.mtx，metadata.csv，gene_names.csv，pca.csv

#需要上一步生成的4个文件，loom.file.txt文件，输出文件dynamic.h5ad
python Anndata_v2.py counts.mtx metadata.csv gene_names.csv pca.csv loom.file.txt dynamic.h5ad

#
python run_RNA_velocity.py dynamic.h5ad
