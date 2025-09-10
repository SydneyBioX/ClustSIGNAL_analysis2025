#### simulation generation functions ####

sim_patch <- function(coords) {
    # generate Patched spatial pattern
    # coords = data frame with two columns, for x and y coordinates 
    
    celltypes <- apply(coords, 1, function(crd) {
        x <- crd[1]
        y <- crd[2]
        # partitioning space by xy-coordinates and allocating labels
        if (x <= 100 & y > 200 & y <= 300) {
            return("CellType1")
        } else if (x > 100 & x <= 200 & y > 200 & y <= 300) {
            return("CellType2")
        } else if (x > 200 & x <= 300 & y > 200 & y <= 300) {
            return("CellType3")
        } else if (x <= 100 & y > 100 & y <= 200) {
            return("CellType4")
        } else if (x > 100 & x <= 200 & y > 100 & y <= 200) {
            return("CellType5")
        } else if (x > 200 & x <= 300 & y > 100 & y <= 200) {
            return("CellType6")
        } else if (x <= 100 & y <= 100) {
            return("CellType7")
        } else if (x > 100 & x <= 200 & y <= 100) {
            return("CellType8")
        } else if (x > 200 & x <= 300 & y <= 100) {
            return("CellType9")
        } 
    })
    return(celltypes)
}

sim_gradient <- function(coords) {
    # generate Gradient spatial pattern
    # coords = data frame with two columns, for x and y coordinates 
    
    celltypes <- apply(coords, 1, function(crd) {
        x <- crd[1]
        y <- crd[2]
        # partitioning space by xy-coordinates and allocating labels
        if (x + y <= 90) {
            return("CellType1")
        } else if (x + y > 90 & x + y <= 110) {
            return("CellType2")
        } else if (x + y > 110 & x + y <= 150) {
            return("CellType3")
        } else if (x + y > 150 & x + y <= 170) {
            return("CellType4")
        } else if (x + y > 170 & x + y <= 230) {
            return("CellType5") 
        } else if (x + y > 230 & x + y <= 240) {
            return("CellType6")
        } else if (x + y > 240 & x + y <= 280) {
            return("CellType7") 
        } else if (x + y > 280 & x + y <= 310) {
            return("CellType8")
        } else if (x + y > 310 & x + y <= 350) {
            return("CellType9")  
        } else if (x + y > 350 & x + y <= 380) {
            return("CellType10")  
        } else if (x + y > 380 & x + y <= 390) {
            return("CellType11")
        } else if (x + y > 390 & x + y <= 450) {
            return("CellType12")
        } else if (x + y > 450 & x + y <= 480) {
            return("CellType13")
        } else if (x + y > 480) {
            return("CellType14") 
        }
    })
    return(celltypes)                 
}

sim_complex <- function(coords, sp_center) {
    # generate Gradient spatial pattern
    # coords = data frame with two columns, for x and y coordinates 
    # sp_center = numeric vector specifying coordinates of center of 
    #   xy-coordinate space
    
    celltypes <- apply(coords, 1, function(crd) {
        x <- crd[1]
        y <- crd[2]
        distance <- sqrt((x - sp_center[1])^2 + (y - sp_center[2])^2)
        # partitioning space by xy-coordinates and/or distance from center 
        # and allocating labels
        if (abs(x + y - 300) < 10 ) {
            return("CellType1") 
        } else if (abs(x - y) < 10) {
            return("CellType2") 
        } else if (distance < 70 & x > y) {
            return("CellType3")
        } else if (distance < 70 & x < y) {
            return("CellType4")
        } else if (distance < 100 & x > y) {
            return("CellType5")
        } else if (distance < 100 & x < y) {
            return("CellType6")
        } else if (distance < 120 & x > y) {
            return("CellType7")
        } else if (distance < 120 & x < y) {
            return("CellType8")
        } else {
            return(sample(c("CellType9", "CellType10", "CellType11", "CellType12"), 1))
        }
    })
    return(celltypes)
}

sim_uniform <- function(coords, n_celltypes) {
    # generate uniform spatial pattern
    # coords = data frame with two columns, for x and y coordinates
    # n_celltypes = numeric value specifying the number of labels to allocate
    
    row_ind <- c(seq_len(nrow(coords)))
    # determining cell pool size for each label, i.e., the number of cells 
    # allotted to a label
    size <- nrow(coords) %/% n_celltypes # get the quotient for sample size
    celltypes <- data.frame(matrix(nrow = 0, ncol = 2, 
                                   dimnames = list(NULL, c("index", "cell"))))
    for (j in seq_len(n_celltypes)) {
        # randomly select specific number of cells and allocate a label
        coord_ind <- sample(row_ind, size, replace = FALSE)
        celltypes <-rbind(celltypes, cbind(index = coord_ind, 
                                           cell = paste0("CellType", j)))
        row_ind <- row_ind[!(row_ind %in% coord_ind)]
    }
    celltypes <-rbind(celltypes, cbind(index = row_ind, 
                                       cell = paste0("CellType", n_celltypes)))
    celltypes$index <- as.integer(celltypes$index)
    celltypes <- celltypes[order(celltypes$index), ]
    return(celltypes$cell)
}


#### simulation benchmarking functions ####

smoothing100 <- function(spe, nnCells, threads = 1) {
    # function for complete smoothing and clustering
    # spe - spe object containing log normalised counts
    # nnCells - matrix containing ids of nearest neighbours of each cell
    # threads - number of cores
    
    #### smoothing
    G_mat <- as(logcounts(spe), "sparseMatrix") # log normalised counts
    wt_val <- 1 / ncol(nnCells) # scaled weight
    # melt matrix and recreate as sparse matrix - for computational efficiency
    wts_df <- cbind(reshape2::melt(t(nnCells))[, 2:3], 
                    rep(wt_val, ncol(nnCells) * nrow(nnCells)))
    colnames(wts_df) <- c("c", "ci", "wci")
    wts_df$c <- as.integer(match(wts_df$c, colnames(spe)))
    wts_df$ci <- as.integer(match(wts_df$ci, colnames(spe)))
    W_mat <- sparseMatrix(i = wts_df$ci, j = wts_df$c,
                          x = wts_df$wci)
    # perform complete smoothing
    smoothMat <- G_mat %*% W_mat
    colnames(smoothMat) <- colnames(spe)
    # add smoothed data to spe object
    assay(spe, "smoothed.100") <- smoothMat
    
    #### clustering
    # dimension reduction of smoothed data
    spe <- runPCA(spe, assay.type = "smoothed.100", name = "PCA.smooth100")
    mat <- reducedDim(spe, "PCA.smooth100")
    # applying default values used in ClustSIGNAL for clustering
    clustVal <- min(as.integer(ncol(spe) / 5), 5000)
    reClust <- clusterRows(mat, TwoStepParam(
        first = KmeansParam(centers = clustVal, iter.max = 30),
        second = NNGraphParam(k = 10, cluster.fun = "louvain", num.threads = threads)))
    spe$reClust100 <- factor(reClust)
    return(spe)
}

smoothing0 <- function(spe, k){
    # function for no smoothing and clustering
    # spe - spe object containing PCA low embedding
    # k - number of clusters to identify
    
    mat <- reducedDim(spe, "PCA")
    nClust <- clusterRows(mat, KmeansParam(centers = k, iter.max = 30))
    spe$reClust0 <- factor(nClust)
    return(spe)
}


#### clustering method functions ####

mbkmeans_clust <- function(spe, dimRed = "PCA", cluster_summary){
    # performing mbkmeans clustering and subclustering
    # spe - spe object containing log normalised counts
    # dimRed - matrix of low embedding
    # cluster_summary - data frame of cluster labels and number of subclusters 
    #   within each cluster
    
    #### clustering
    k <- dim(cluster_summary)[1] # number of clusters
    res <- mbkmeans(spe, clusters = k, reduceMethod = dimRed) # run mbkmeans
    # add cluster labels to spe object as "initial clusters"
    spe$initCluster <- res$Clusters 
    print("Mbkmeans clustering performed.")
    
    #### subclustering
    subclusters_list <- list()
    for (cl in sort(unique(spe$initCluster))) {
        speX <- spe[, spe$initCluster == cl] # subset spe object to cluster cells
        # get number of subclusters to identify in the cluster
        k_num <- cluster_summary[cluster_summary$clusters == cl, 
                                 "num_subclusters"]
        res_sub <- mbkmeans(speX, clusters = k_num, reduceMethod = NA,
                            whichAssay = "logcounts") # run mbkmeans
        subclusters_list[[cl]] <- setNames(paste0(cl, ".", res_sub$Cluster),
                                           colnames(speX))
    }
    # add subcluster labels to spe object as "initial subclusters"
    spe$initSubcluster <- unlist(subclusters_list)[colnames(spe)]
    print("Mbkmeans subclustering performed.")
    return(spe)
}

walktrap_clust <- function(spe, dimRed = "PCA") {
    # performing walktrap clustering and subclustering
    # spe - spe object containing log normalised counts
    # dimRed - matrix of low embedding
    
    #### clustering
    # build neighbourhood graph
    data_graph <- buildSNNGraph(spe, use.dimred = dimRed, k = 10)
    # perform walktrap community detection and add cluster labels to spe object 
    # as "initial clusters"
    spe$initCluster <- cluster_walktrap(data_graph)$membership
    print("Walktrap clustering performed.")
    
    #### subclustering
    subclusters_list <- list()
    for (cl in sort(unique(spe$initCluster))) {
        speX <- spe[, spe$initCluster == cl] # subset spe object to cluster cells
        # build neighbourhood graph
        data_graphX <- buildSNNGraph(speX, use.dimred = NULL, k = 10, 
                                     assay.type = "logcounts")
        # perform walktrap community detection
        clusters_out <- cluster_walktrap(data_graphX)$membership
        subclusters_list[[cl]] <- setNames(paste0(cl, ".", clusters_out),
                                           colnames(speX))
    }
    # add subcluster labels to spe object as "initial subclusters"
    spe$initSubcluster <- unlist(subclusters_list)[colnames(spe)]
    print("Walktrap subclustering performed.")
    return(spe)
}

dbscan_clust <- function(spe) {
    # performing walktrap clustering and subclustering
    # spe - spe object containing low embedding data
    
    #### clustering
    mat <- reducedDim(spe, "PCA")
    # perform dbscan clustering and add cluster labels to spe object 
    # as "initial clusters"
    spe$initCluster <- clusterRows(mat, DbscanParam()) |> as.numeric()
    # replace NA with 999 to avoid errors during subclustering
    spe$initCluster <- replace_na(spe$initCluster, 999)
    print("Dbscan clustering performed.")
    
    #### subclustering
    subclusters_list <- list()
    for (cl in sort(unique(spe$initCluster))) {
        speX <- spe[, spe$initCluster == cl] # subset spe object to cluster cells
        speX <- runPCA(speX) # dimension reduction
        matX <- reducedDim(speX, "PCA")
        # perform dbscan clustering
        clusters_out <- clusterRows(matX, DbscanParam()) |> as.numeric()
        # replace NA with 999 to avoid errors during downstream analysis
        clusters_out <- replace_na(clusters_out, 999)
        subclusters_list[[cl]] <- setNames(paste0(cl, ".", clusters_out),
                                           colnames(speX))
    }
    # add subcluster labels to spe object as "initial subclusters"
    spe$initSubcluster <- unlist(subclusters_list)[colnames(spe)]
    print("Dbscan subclustering performed.")
    return(spe)
}