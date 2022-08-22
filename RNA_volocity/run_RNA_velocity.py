import scanpy as sc
import scvelo as scv
import pandas as pd
import numpy as np
import os
import anndata
import cellrank as cr
import matplotlib.pyplot as plt
import sys

args = sys.argv
# read h5ad file
adata= scv.read(args[1])

scv.pl.proportions(adata, groupby='celltype', save = "proportions.pdf")
scv.pl.proportions(adata, groupby='celltype', save = "proportions.png")
#scv.pl.proportions(adata)
plt.savefig("proportions.pdf")
sc.pl.umap(adata, color="celltype", save="_celltype.pdf")
print(adata)
# Running RNA Velocity

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization

scv.pp.filter_and_normalize(adata,min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=15, n_neighbors=30)
#scv.tl.velocity(adata, mode = "stochastic")
scv.tl.recover_dynamics(adata, n_jobs=10)
scv.tl.velocity(adata, mode = "dynamical")
scv.tl.velocity_graph(adata)

#scv.pl.velocity_embedding(adata, basis='X_umap', arrow_size=5)
#ident_colours = ["#F8766D", "#A3A500", "#00BF7D", "#00B0F6", "#E76BF3"]
scv.pl.velocity_embedding(adata, basis='umap', color = "celltype",  arrow_length=3, arrow_size=2, dpi=480, save = "embedding.pdf")
scv.pl.velocity_embedding(adata, basis='umap', color = "celltype",  arrow_length=3, arrow_size=2, dpi=480, save = "embedding.png")

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color = 'celltype',  save = "embedding_streaming.pdf")
scv.pl.velocity_embedding_stream(adata, basis = 'umap', color = 'celltype',  save = "embedding_streaming.png")

scv.pl.velocity_embedding_grid(adata, basis = 'umap', color = 'celltype',  save = "embedding_grid.pdf")
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color = 'celltype',  save = "embedding_grid.png")

#plt.savefig("velocity_embedding_stream.png")
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8)
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8, save="UMAP_stream.png", figsize=(7,5), legend_fontsize = 9, show=False, title='')
#scv.pl.velocity_embedding_stream(adata, basis='X_umap',color = "celltype", palette = ident_colours, size = 20,alpha =0.8, save="UMAP_stream.pdf", figsize=(7,5), legend_fontsize = 9, show=False, title='')

#cr.tl.terminal_states(adata, cluster_key="seurat_clusters", weight_connectivities=0.2)
#cr.pl.terminal_states(adata, save = "terminal_states.png")

scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)
os.makedirs("tables", exist_ok=True)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
write = pd.ExcelWriter("tables/rank_velocity_genes.xls")
df.to_excel(write)
write.save()

df.to_html("tables/rank_velocity_genes.html", border = 1)
df.to_csv("tables/rank_velocity_genes.csv")
print (df.head())

for cluster in df.columns:
	kwargs = dict(frameon=False, size=10, linewidth=1.5, add_outline = cluster, color="celltype")
	scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, figsize = (15, 9), save = cluster + "_scatter.pdf", **kwargs)
	scv.pl.scatter(adata, df[cluster][:5], ylabel=cluster, figsize = (15, 9), save = cluster + "_scatter.png", **kwargs)
	scv.pl.velocity(adata, df[cluster][:4], save = cluster + "_velocity.pdf", ncols=2)
	scv.pl.velocity(adata, df[cluster][:4], save = cluster + "_velocity.png", ncols=2)
        #gene = df['B_cell'][:3]
	#print (gene)
	#scv.pl.scatter(adata, gene, save="cluste_scatter.pdf", **kwargs)

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="length_velocity_confidence.pdf")
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95], save="length_velocity_confidence.png")


scv.tl.velocity_pseudotime(adata)
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot')

adata.uns['neighbors']['distances'] = adata.obsp['distances']
adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']

scv.tl.paga(adata, groups='celltype')
#df = scv.get_df(adata, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')

print ("Fine")
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save="paga.pdf")
scv.pl.paga(adata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5, save="paga.png")

'''
df = adata_subset.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
	pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
	pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
	pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

	plt.savefig("hist.pdf")
	plt.savefig("hist.png")

scv.get_df(adata_subset, 'fit*', dropna=True).to_csv("adata_subset.csv")
'''

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save = "latent_time_scatter.pdf")
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80, save = "latent_time_scatter.png")

