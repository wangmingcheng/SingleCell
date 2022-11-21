#!/usr/bin/env Rscript
setwd('~/analysis')
###########################
id='ROSMAP_BBB_mouse.integration.CCA'
##########################
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
wr <- colorRampPalette(colors = c( "white", "red"))(100)
rwb <- colorRampPalette(colors = c("blue", "white", "red"))(100)
ryb=colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)
prgn=colorRampPalette(rev(brewer.pal(n = 10, name ="PRGn")))(100)
wb <- colorRampPalette(colors = c( "white", "blue"))(100)
#####################################
ids=read.table('~/Dropbox (MIT)/work/dataset/human2mouse.geneSymbol.txt',sep='\t')
## mouse
mouse=readRDS('GSE98816.mouse.brain_vascular.seurat.rds')
meta=mouse@meta.data
table(meta$celltype)
mouse=subset(mouse,cells=rownames(meta[meta$celltype %in% c('Endothelial','Fibroblast','Pericyte','SmoothMuscleCell'),]))
genes=rownames(mouse)
newgenes=c()
for (g in genes){
  if (g %in% ids$V2){
    newg=as.character(ids[ids$V2==g,1])
  }else{newg=toupper(g)}
  newgenes=c(newgenes,newg)
}
rownames(mouse@assays$RNA@counts)=newgenes
rownames(mouse@assays$RNA@data)=newgenes
rownames(mouse@assays$RNA@scale.data)=newgenes
mouse$orig.ident=rep('mouse',nrow(mouse@meta.data))
mouse$species=rep('mouse',nrow(mouse@meta.data))
mouse$cellsubtype=rep('undefined',nrow(mouse@meta.data))
meta=mouse@meta.data
count=mouse@assays$RNA@counts
mouse= CreateSeuratObject(counts = count, project = "mouse", meta.data = meta)
mouse= NormalizeData(mouse, normalization.method = "LogNormalize", scale.factor = 10000)
mouse=FindVariableFeatures(mouse, selection.method = "vst", nfeatures = 2000)
mouse
##
## ROSMAP
rosmap=readRDS('ROSMAP.VascularCells.seurat.rds')
rosmap$orig.ident=rep('Frozen',nrow(rosmap@meta.data))
rosmap$species=rep('human',nrow(rosmap@meta.data))
new=gsub('Endo','Endothelial',rosmap$celltype)
new=gsub('Fib','Fibroblast',new)
new=gsub('SMC','SmoothMuscleCell',new)
new=gsub('Per','Pericyte',new)
rosmap$celltype=new
rosmap$cellsubtype=rosmap$subtype
## FRESH
bbb=readRDS('brain.BBB.human.vascular.rds')
bbb
head(bbb@meta.data)
bbb$orig.ident=rep('Fresh',nrow(bbb@meta.data))
bbb$species=rep('human',nrow(bbb@meta.data))
new=gsub('aEndo','Endothelial',bbb$cellsubtype)
new=gsub('capEndo','Endothelial',new)
new=gsub('vEndo','Endothelial',new)
new=gsub('Fib1','Fibroblast',new)
new=gsub('Fib2','Fibroblast',new)
new=gsub('Fib3','Fibroblast',new)
new=gsub('aSMC','SmoothMuscleCell',new)
new=gsub('vSMC','SmoothMuscleCell',new)
new=gsub('Per1','Pericyte',new)
new=gsub('Per2','Pericyte',new)
table(new)
bbb$celltype=new
bbb@meta.data=bbb@meta.data[,c(1:9,17:19,20)]
####3
reference.list <- c(rosmap,bbb,mouse)
brain.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)

brain <- IntegrateData(anchorset = brain.anchors, dims = 1:30)

library(ggplot2)
library(cowplot)
library(patchwork)
# switch to integrated assay. The variable features of this assay are automatically
# set during IntegrateData
DefaultAssay(brain) <- "integrated"

# Run the standard workflow for visualization and clustering
brain <- ScaleData(brain, verbose = FALSE)
brain <- RunPCA(brain, npcs = 30, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)
brain <- FindNeighbors(brain, dims = 1:30)
brain <- FindClusters(brain, resolution = 0.5)

DimPlot(brain, reduction = "umap", group.by = "orig.ident")
DimPlot(brain, reduction = "umap", group.by = "celltype", label = TRUE, repel = TRUE) + NoLegend()
DimPlot(brain, reduction = "umap", group.by = "cellsubtype", label = TRUE, repel = TRUE) + NoLegend()

brain=subset(brain,idents = c(0:8,10))

pdf(file=paste(id,'.umap.pdf',sep=''),height = 6,width = 7)
DimPlot(brain, reduction = "umap", group.by = "orig.ident",cols=col.all[1:length(table(brain$orig.ident))])
DimPlot(brain, reduction = "umap",label = T)
DimPlot(brain, reduction = "umap", group.by = "celltype", label = TRUE,repel = TRUE,cols = col.all[1:length(table(brain$celltype))])
DimPlot(brain, reduction = "umap", group.by = "celltype", label = F,repel = TRUE,cols = col.all[1:length(table(brain$celltype))])
DimPlot(brain, reduction = "umap", group.by = "cellsubtype", label = TRUE,repel = TRUE)
DimPlot(brain, reduction = "umap", group.by = "cellsubtype", label = F,repel = TRUE)
DimPlot(brain, reduction = "umap", group.by = "cellsubtype", label = TRUE,repel = TRUE,cols = col.all[1:length(table(brain$cellsubtype))])
dev.off()

pdf(file=paste(id,'.umap_cluster.pdf',sep=''),height = 6,width = 7)
DimPlot(brain, reduction = "umap",label = T)
DimPlot(brain, reduction = "umap",label = F)
dev.off()

saveRDS(brain,file=paste(id,'.rds',sep=''))
write.table(brain@meta.data,file=paste(id,'.metadata.txt',sep=''),quote = F,sep='\t')
########### 

############## DEGs between human and mouse ########3
DefaultAssay(brain) <- "RNA"

brain.sel=brain


brain.sel$celltype.ident=paste(brain.sel$celltype,brain.sel$orig.ident,sep='_')
Idents(brain.sel)=brain.sel$celltype.ident

ctps=names(table(brain.sel$celltype))
ids=names(table(brain.sel$orig.ident))

ids
## Fresh vs Frozen
res=c()
alldegs=c()
for (c in ctps){
  id1=paste(c,ids[1],sep='_')
  id2=paste(c,ids[2],sep='_')
  degs=FindMarkers(brain.sel,ident.1=id1,ident.2 = id2)
  degs$gene=rownames(degs)
  degs$celltype=rep(c,nrow(degs))
  res=rbind(res,head(degs))
  alldegs=rbind(alldegs,degs)
}
table(alldegs$celltype)
write.table(alldegs,file='Fresh_Frozen.DEGs.txt',sep = '\t',quote = F)
## Fresh vs Mouse
res=c()
alldegs=c()
for (c in ctps){
  id1=paste(c,ids[1],sep='_')
  id2=paste(c,ids[3],sep='_')
  degs=FindMarkers(brain.sel,ident.1=id1,ident.2 = id2)
  degs$gene=rownames(degs)
  degs$celltype=rep(c,nrow(degs))
  res=rbind(res,head(degs))
  alldegs=rbind(alldegs,degs)
}
table(alldegs$celltype)
write.table(alldegs,file='Fresh_Mouse.DEGs.txt',sep = '\t',quote = F)

## Frozen vs Mouse
res=c()
alldegs=c()
for (c in ctps){
  id1=paste(c,ids[2],sep='_')
  id2=paste(c,ids[3],sep='_')
  degs=FindMarkers(brain.sel,ident.1=id1,ident.2 = id2,logfc.threshold = log(2))
  degs$gene=rownames(degs)
  degs$celltype=rep(c,nrow(degs))
  res=rbind(res,head(degs))
  alldegs=rbind(alldegs,degs)
}
table(alldegs$celltype)
write.table(alldegs,file='Frozen_Mouse.DEGs.txt',sep = '\t',quote = F)

##
