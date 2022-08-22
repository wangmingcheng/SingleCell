import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import loompy
import sys

args = sys.argv
# load sparse matrix:
X = io.mmread(args[1])
print (args[1])

# create anndata object
adata = anndata.AnnData(X=X.transpose().tocsr(), dtype=np.float32)

# load cell metadata:
cell_meta = pd.read_csv(args[2])

# load gene names:
with open(args[3], 'r') as f:
    gene_names = f.read().splitlines()

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv(args[4])
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T

# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=['sample'], frameon=False, save=True)

# save dataset as anndata format
#adata.write('my_data.h5ad')

#adata = sc.read_h5ad('my_data.h5ad')
#merge loom
files = []
inf = open(args[5], "r")
for line in inf:
        if not line.strip().startswith("#"):
                files.append(line.strip())
inf.close()

loompy.combine(files, "merged.loom", key="Accession")

#load loom files for spliced/unspliced matrices for each sample:
ldata = scv.read('merged.loom', cache=True)
# rename barcodes in order to merge:
barcodes = [bc.split(':')[1] for bc in ldata.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_10' for bc in barcodes]
ldata.obs.index = barcodes

# make variable names unique
ldata.var_names_make_unique()
# merge matrices into the original adata object
adata = scv.utils.merge(adata, ldata)

adata.write(args[6], compression = 'gzip')
