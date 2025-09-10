#### wrapper functions for method packages ####
runBANKSY <- function(spe, annots_label, sample_label, SEED, k_geom = c(15, 30), 
                      batch = FALSE, batch_by = "None", lambda = 0.2, res = 1, 
                      npcs = 20, use_agf = TRUE, compute_agf = TRUE){
    
    #### parameters
    # spe - spe object containing raw counts
    # annots_label - metadata column name containing annotation labels
    # sample_label - metadata column name containing sample IDs
    # batch - whether to perform batch correction
    # batch_by - metadata column name containing batch groups
    
    #### Banksy parameters
    # compute_agf - TRUE computes both weighted neighborhood mean (H_0) and
    #   the azimuthal Gabor filter (H_1).
    # lambda - mixing parameter, ranges from 0-1. Smaller lambda for cell-typing 
    #   mode (recommended value is 0.2)
    # k_geom - numeric value or vector defining neighbourhood sizes of H_0 and 
    #   H_1, respectively. Recommended values c(15, 30).
    # res - Leiden clustering resolution. Higher value gives more clusters.
    # SEED - to set.seed() for reproducibility.
    # npcs - number of PCA dimensions to calculate. Default 20.
    
    if (batch == TRUE) {
        show(paste("Multisample run with batch correction. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        colnames(spe) <- paste0(colnames(spe), "_", spe[[sample_label]])
        
        # Staggering spatial coordinates
        locs <- spatialCoords(spe)
        locs <- cbind(locs, sample = factor(spe[[sample_label]]))
        locs_dt <- data.table(locs)
        colnames(locs_dt) <- c("sdimx", "sdimy", "group")
        locs_dt[, sdimx := sdimx - min(sdimx), by = group]
        global_max <- max(locs_dt$sdimx) * 1.5
        locs_dt[, sdimx := sdimx + group * global_max]
        locs <- as.matrix(locs_dt[, 1:2])
        rownames(locs) <- colnames(spe)
        spatialCoords(spe) <- locs
        show(paste("Spatial coordinates of samples staggered. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Following normalization approach in vignette
        seu <- as.Seurat(spe, data = NULL)
        # normalizing data
        scale_factor <- median(colSums(assay(spe, "counts")))
        seu <- NormalizeData(seu, scale.factor = scale_factor,
                             normalization.method = "RC")
        # Adding data to spe object
        assay(spe, "normcounts") <- GetAssayData(seu)
        show(paste("Seurat normalisation complete. Time",
                   format(Sys.time(),'%H:%M:%S')))
        
        # Running BANKSY
        show(paste("BANKSY run started. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        spe <- computeBanksy(spe, assay_name = "normcounts",
                             compute_agf = compute_agf, k_geom = k_geom)
        spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda,
                            npcs = npcs, seed = SEED)
        # Harmony batch correction
        PCA_label <- paste0("PCA_M", as.numeric(use_agf), "_lam", lambda)
        set.seed(SEED)
        harmony_embedding <- RunHarmony(data_mat = reducedDim(spe, PCA_label),
                                        meta_data = colData(spe),
                                        vars_use = batch_by,
                                        verbose = FALSE)
        reducedDim(spe, "PCA_harmony") <- harmony_embedding
        show(paste("Batch correction completed. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        # Banksy clustering
        spe <- clusterBanksy(spe, dimred = "PCA_harmony", use_agf = use_agf,
                             lambda = lambda, resolution = res, seed = SEED)
        show(paste("BANKSY clustering completed. Time", 
                   format(Sys.time(),'%H:%M:%S')))
    }
    else {
        # Following normalization approach in vignette
        # separating samples into individual spe objects
        sample_names <- unique(spe[[sample_label]])
        spe_list <- lapply(sample_names, function(x) spe[, spe[[sample_label]] == x])
        # Seurat - normalizing data
        seu_list <- lapply(spe_list, function(x) {
            x_seu <- as.Seurat(x, data = NULL)
            x_seu <- NormalizeData(x_seu, scale.factor = 5000,
                                   normalization.method = "RC")
            return(x_seu)})
        # Adding data to spe object
        spe_list <- Map(function(spe, seu) {
            assay(spe, "normcounts") <- GetAssayData(seu)
            spe},
            spe_list, seu_list)
        show(paste("Seurat normalisation complete. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Running BANKSY
        show(paste("BANKSY run started. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        spe_list <- lapply(spe_list, computeBanksy, assay_name = "normcounts",
                           compute_agf = compute_agf, k_geom = k_geom)
        # merging samples for downstream steps
        spe <- do.call(cbind, spe_list)
        spe <- runBanksyPCA(spe, use_agf = use_agf, lambda = lambda,
                            group = sample_label, seed = SEED)
        spe <- clusterBanksy(spe, use_agf = use_agf, lambda = lambda,
                             resolution = res, seed = SEED)
        show(paste("BANKSY clustering completed. Time", 
                   format(Sys.time(),'%H:%M:%S')))
    }
    return(spe)
}


runBASS <- function(cntm, xym, C, R) {
    
    #### BASS parameters
    # cntm - list of counts matrices
    # xym - list of xy-coordinate matrices
    # C - number of expected cell types
    # R - number of expected domains
    
    # Set up BASS object
    BASS <- createBASSObject(cntm, xym, C, R, beta_method = "SW", 
                             init_method = "kmeans")
    # Data pre-processing
    BASS <- BASS.preprocess(BASS, doLogNormalize = TRUE, doPCA = TRUE, 
                            scaleFeature = TRUE, nPC = 20)
    # Run BASS algorithm
    BASS <- BASS.run(BASS)
    
    # post-process posterior samples
    BASS <- BASS.postprocess(BASS)
    return (BASS)
}


runSpatialPCA <- function(spe, samples, sparkv = 'sparkx', ncores = 1) {
    
    #### parameters
    # spe - spe object with one or more samples
    # samples - metadata column containing sample IDs
    
    #### SpatialPCA parameters
    # sparkv - 'sparkx' for large datasets,
    # ncores - number of cpu cores to use
    
    sample_IDs <- unique(spe[[samples]]) # get sample IDs
    # apply SpatialPCA to each sample
    clust_list <- lapply(sample_IDs, function(i) {
        speX <- spe[, spe[[samples]] == i] # subset spe object
        cnts <- counts(speX) # raw counts
        locs <- spatialCoords(speX) # spatial locations
        
        # Create SpatialPCA object
        spca <- CreateSpatialPCAObject(counts = cnts, location = locs, 
                                       sparkversion = sparkv, 
                                       numCores_spark = ncores)
        show(paste("Created SpatialPCA object for", i, ". Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Build SpatialPCA kernel
        spca <- SpatialPCA_buildKernel(spca, kerneltype = 'gaussian',
                                       bandwidthtype = 'Silverman', 
                                       sparseKernel = TRUE,
                                       sparseKernel_tol = 1e-5,
                                       sparseKernel_ncore = ncores)
        show(paste("Built SpatialPCA kernel. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Estimate loadings
        spca <- SpatialPCA_EstimateLoading(spca)
        show(paste("Estimated SpatialPCA loadings. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Generate SpatialPCA embedding
        spca <- SpatialPCA_SpatialPCs(spca)
        show(paste("Calculated embeddings. Time", 
                   format(Sys.time(),'%H:%M:%S')))
        
        # Clustering SpatialPCA embedding using ClustSIGNAL's default parameters
        clustVal <- min(as.integer(ncol(speX) / 5), 5000)
        clusters <- clusterRows(t(spca@SpatialPCs), TwoStepParam(
            first = KmeansParam(centers = clustVal, iter.max = 30),
            second = NNGraphParam(k = 10, num.threads = ncores, 
                                  cluster.fun = "louvain")))
        names(clusters) <- colnames(spca@normalized_expr) # adding cell IDs
        
        # Accounting for unlabelled cells
        unlab_cells <- setdiff(colnames(speX), colnames(spca@normalized_expr))
        unlab_clusts <- rep("Unlabelled", length(unlab_cells))
        names(unlab_clusts) <- unlab_cells
        
        # Adding unlabelled cells to final result
        clusters <- append(clusters, unlab_clusts) 
        
        # Reordering cells according to their index in spe object
        clusters <- clusters[colnames(speX)] 
        show(paste("Clustering performed on", i, ". Time", 
                   format(Sys.time(),'%H:%M:%S')))
        return(clusters)
    })
    return(clust_list)
}