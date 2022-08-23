library("Seurat")

ChooseClusterResolutionDownsample <- function(input.srobj, res.low = .1, res.high = 1.2, by = 0.02, n.pcs = 30, bias = "over", figdir) {

	######## step 1: save the input seurat object as a new temporary object, 
	########         dont want to overwrite or change the original one with all of the parameter scans
    
	srobj.tmp <- input.srobj 

    	######## step 2: calculate the FindClusters over a large range of resolutions
	print("Performing parameter scan over multiple resolutions...")

    	#set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=3)
    	set.res <- seq(from = res.low, to = res.high, by = by)
    	srobj.tmp <- FindClusters(srobj.tmp, resolution = set.res)
    
    	######## step 3: output plot of how the resolution changes the number of clusters you get
    	snn <- ifelse(DefaultAssay(srobj.tmp) == "SCT", "SCT_snn_res.", "RNA_snn_res.")
	n.clusters = vector(mode = "numeric", length = length(set.res))
    	names(n.clusters) = set.res
    	for(i in 1:length(n.clusters)){
        	n.clusters[i] = length(table(as.vector(srobj.tmp@meta.data[, paste(snn, names(n.clusters)[i], sep = "")])))
    	}
	
    	######## step 4: calculate the silhouette width for each resolution
    	print("Computing a silhouette width for each cell, for each resolution...")
    	require(cluster)

    	dist.temp <- cor(t(srobj.tmp@reductions$pca@cell.embeddings[, 1:n.pcs]), method = "pearson")
    	random.cells.choose <- sample(1:nrow(dist.temp), round(nrow(dist.temp)/10, digits = 0))
    	dist.temp.downsample <- dist.temp[random.cells.choose, random.cells.choose]
    	sil.all.matrix = matrix(data=NA, nrow = nrow(dist.temp.downsample), ncol = 0)

    	#print(names(srobj.tmp[[]]))

    	for(i in 1:length(set.res)){
        	clusters.temp = as.numeric(as.vector(srobj.tmp@meta.data[random.cells.choose, paste(snn, set.res[i], sep = "")]))
        	if(length(table(clusters.temp)) > 1){
            		sil.out = silhouette(clusters.temp, as.dist(1 - as.matrix(dist.temp.downsample)))
            		sil.all.matrix = cbind(sil.all.matrix, sil.out[, 3])
        	}
        	if(length(table(clusters.temp)) == 1){
            		sil.all.matrix = cbind(sil.all.matrix, rep(0, length(clusters.temp)))
        	}
        	print(paste("          ", round(100*i/length(set.res)), "% done with silhouette calculation", sep = ""))
    
    	}

    	######## step 5: calculate summary metric to compare the silhouette distributions,
    	########         average has worked well so far... could get fancier

    	print("Identifying a best resolution to maximize silhouette width")
    	sil.average = colMeans(sil.all.matrix)
    	names(sil.average) = set.res


    	######## step 6: automate choosing resolution that maximizes the silhouette 
    	hist.out = hist(sil.average, length(sil.average)/1.2,  plot = FALSE)
    
    	#  take the ones that fall into the top bin, 
    	#  and the max OR MIN of those  ******* can change this to under vs over cluster
    	if(bias == "over"){
        	resolution.choice = as.numeric(max(names(sil.average[which(sil.average > hist.out$breaks[length(hist.out$breaks) - 1])])))
    	}
    	if(bias == "under"){
        	resolution.choice = as.numeric(min(names(sil.average[which(sil.average > hist.out$breaks[length(hist.out$breaks) - 1])])))
    	}
    
    	# get the silhouette of the best resolution: 
    	silhouette.best = as.numeric(sil.average[paste(resolution.choice)])

    	print(paste("Best Resolution Choice: ", resolution.choice, ", with average silhouette score of: ",
        	round(silhouette.best, digits=3), ", giving ", as.numeric(n.clusters[paste(resolution.choice)]), " clusters", sep = ""))

    	######### step 7: output plot and data 
        
        print(paste0("Ouptutting summary statistics and returning seurat object... ",
              "This will create a pdf in your output directory,",
              " and will return your input seurat object ammended with the best choice",
              " for clusters (found as Best.Clusters in the meta.data matrix, and set to your new ident)..."))

        pdf(paste(figdir, "/", "resolution_benchmark.pdf", sep = ""),
        width = 10, height = 4, useDingbats = FALSE)
        par(mfrow = c(1,3))
        # Resolution vs # of Clusters
        plot(set.res, n.clusters, col = "black", pch = 20,
             type = "p", xlab = "Resolution", ylab = "Clusters",
             main = "Resolution vs. Clusters")
        # Resolution vs Average Silhouette
        plot(set.res, sil.average, col = "black", pch = 20,
	     type = "p", xlab = "Resolution", ylab = "Average Silhouette",
             main = "Resolution vs. Average Silhouette")
        abline(h = hist.out$breaks[length(hist.out$breaks) - 1], col = "firebrick3", lty = 2)
        abline(v = resolution.choice, col = "dodgerblue2", lty = 2)

        # N Clusters vs Average Silhouette
        plot(n.clusters, sil.average, col = "black", pch = 20,
             type = "p", xlab = "Clusters", ylab = "Average Silhouette",
             main = "Clusters vs. Average Silhouette")
        abline(h = hist.out$breaks[length(hist.out$breaks) - 1], col = "firebrick3", lty = 2)
        abline(v = as.numeric(n.clusters[paste(resolution.choice)]), col = "dodgerblue2", lty = 2)
        dev.off()
    
    	######## step 8: return the original seurat object, with the metadata containing a 
    	########         concatenated vector with the clusters defined by the best choice here,
    	########         as well as the ident set to this new vector
    
    	Best.Clusters = srobj.tmp@meta.data[, paste(snn, resolution.choice, sep = "")]
    
    	input.srobj$seurat_clusters <- Best.Clusters
    	Idents(input.srobj) <- input.srobj$seurat_clusters
    	input.srobj@misc$resolution.choice <- resolution.choice
    	return(input.srobj)
}

obj <- readRDS("clustree/result.findClusters.rds")
obj <- ChooseClusterResolutionDownsample(obj, bias = "over", figdir = "result")
saveRDS(obj, "optimal_resolution.rds")
