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


runBASS <- function(cntm, xym, C, R, batch = T) {
    
    #### BASS parameters
    # cntm - list of counts matrices
    # xym - list of xy-coordinate matrices
    # C - number of expected cell types
    # R - number of expected domains
    
    # Set up BASS object
    BASS <- createBASSObject(cntm, xym, C, R, beta_method = "SW", 
                             init_method = "kmeans")
    # Data pre-processing
    BASS <- BASS.preprocess(BASS, , doBatchCorrect = batch)
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


#### Python method wrappers via reticulate ####

.resolveSampleLabel <- function(sample_label = NULL, samples = NULL) {
    if (!is.null(sample_label) && !is.null(samples) && sample_label != samples) {
        stop("`sample_label` and `samples` must match when both are provided.")
    }
    label <- if (!is.null(sample_label)) sample_label else samples
    if (is.null(label)) {
        stop("Provide a sample column via `sample_label` or `samples`.")
    }
    return(label)
}


.expandBySample <- function(value, sample_names, arg_name) {
    if (is.list(value) && !is.null(names(value))) {
        missing_samples <- setdiff(sample_names, names(value))
        if (length(missing_samples) > 0) {
            stop(paste0("`", arg_name, "` is missing entries for: ",
                        paste(missing_samples, collapse = ", ")))
        }
        out <- unname(value[sample_names])
    }
    else if (!is.null(names(value))) {
        missing_samples <- setdiff(sample_names, names(value))
        if (length(missing_samples) > 0) {
            stop(paste0("`", arg_name, "` is missing entries for: ",
                        paste(missing_samples, collapse = ", ")))
        }
        out <- unname(as.list(value[sample_names]))
    }
    else if (length(value) == 1) {
        out <- rep(list(value), length(sample_names))
    }
    else if (length(value) == length(sample_names)) {
        out <- as.list(value)
    }
    else {
        stop(paste0("`", arg_name, "` must have length 1 or ",
                    length(sample_names), "."))
    }
    return(out)
}


.preparePythonBenchmarkInputs <- function(spe, sample_label, assay_name = "counts") {
    sample_vec <- as.character(spe[[sample_label]])
    sample_names <- unique(sample_vec)
    gene_names <- rownames(spe)
    if (is.null(gene_names)) {
        gene_names <- paste0("gene_", seq_len(nrow(spe)))
    }
    
    spe_list <- lapply(sample_names, function(x) spe[, sample_vec == x])
    counts_list <- lapply(spe_list, function(x) assay(x, assay_name))
    coords_list <- lapply(spe_list, spatialCoords)
    barcode_list <- lapply(spe_list, function(x) {
        ids <- colnames(x)
        if (is.null(ids)) {
            ids <- paste0("cell_", seq_len(ncol(x)))
        }
        as.character(ids)
    })
    
    names(counts_list) <- sample_names
    names(coords_list) <- sample_names
    names(barcode_list) <- sample_names
    
    return(list(
        sample_names = sample_names,
        counts = unname(counts_list),
        coords = unname(coords_list),
        barcodes = unname(barcode_list),
        genes = as.character(gene_names)
    ))
}


.getBenchmarkPythonModule <- function(py_file = NULL, python_path = NULL) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package `reticulate` is required for GraphST/SpaGCN/STAGATE wrappers.")
    }
    if (!is.null(python_path)) {
        reticulate::use_python(python_path, required = TRUE)
    }
    if (is.null(py_file)) {
        py_file <- file.path("Benchmarking_scripts", "benchmarking_methods_py.py")
    }
    py_file <- normalizePath(py_file, mustWork = TRUE)
    module_name <- tools::file_path_sans_ext(basename(py_file))
    reticulate::import_from_path(module_name,
                                 path = dirname(py_file),
                                 convert = TRUE)
}


.pyMatrixToR <- function(x, n_rows = NULL) {
    if (is.null(x)) {
        return(NULL)
    }
    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }
    else if (!is.matrix(x) && !is.array(x)) {
        x <- unlist(x, use.names = FALSE)
    }
    if (is.null(dim(x))) {
        if (is.null(n_rows)) {
            x <- matrix(x, nrow = 1)
        }
        else {
            n_cols <- length(x) / n_rows
            x <- matrix(x, nrow = n_rows, ncol = n_cols, byrow = TRUE)
        }
    }
    return(x)
}


.formatPythonBenchmarkResult <- function(py_res, cell_order) {
    sample_ids <- as.character(unlist(py_res$sample_ids, use.names = FALSE))
    barcodes <- py_res$barcodes
    
    clusters <- setNames(lapply(sample_ids, function(sid) {
        vals <- as.character(unlist(py_res$clusters[[sid]], use.names = FALSE))
        ids <- as.character(unlist(barcodes[[sid]], use.names = FALSE))
        stats::setNames(vals, ids)
    }), sample_ids)
    
    cluster_vector <- unlist(clusters, use.names = FALSE)
    names(cluster_vector) <- unlist(lapply(clusters, names), use.names = FALSE)
    cluster_vector <- factor(cluster_vector[cell_order])
    
    embeddings <- setNames(lapply(sample_ids, function(sid) {
        emb <- py_res$embeddings[[sid]]
        ids <- as.character(unlist(barcodes[[sid]], use.names = FALSE))
        emb <- .pyMatrixToR(emb, n_rows = length(ids))
        if (is.null(emb)) {
            return(NULL)
        }
        rownames(emb) <- ids
        colnames(emb) <- paste0("dim_", seq_len(ncol(emb)))
        emb
    }), sample_ids)
    
    cluster_columns <- setNames(lapply(sample_ids, function(sid) {
        cols <- py_res$cluster_columns[[sid]]
        ids <- as.character(unlist(barcodes[[sid]], use.names = FALSE))
        if (is.null(cols)) {
            return(NULL)
        }
        cols <- lapply(cols, function(x) as.character(unlist(x, use.names = FALSE)))
        out <- as.data.frame(cols, stringsAsFactors = FALSE)
        rownames(out) <- ids
        out
    }), sample_ids)
    
    out <- list(
        method = py_res$method,
        clusters = clusters,
        cluster_vector = cluster_vector,
        embeddings = embeddings,
        cluster_columns = cluster_columns,
        params = py_res$params,
        raw = py_res
    )
    class(out) <- c("benchmarkPythonResult", class(out))
    return(out)
}


runGraphST <- function(spe, sample_label = NULL, samples = NULL, n_clusters, SEED,
                       assay_name = "counts", cluster_method = "mclust",
                       refinement = TRUE, radius = 50, datatype = "Slide",
                       epochs = 600, device = "cpu", start = 0.1, end = 3.0,
                       increment = 0.01, py_file = NULL, python_path = NULL) {
    
    sample_label <- .resolveSampleLabel(sample_label, samples)
    inputs <- .preparePythonBenchmarkInputs(spe, sample_label, assay_name)
    n_clusters <- .expandBySample(n_clusters, inputs$sample_names, "n_clusters")
    datatype <- .expandBySample(datatype, inputs$sample_names, "datatype")
    py_mod <- .getBenchmarkPythonModule(py_file, python_path)
    
    show(paste("GraphST run started. Time", format(Sys.time(), '%H:%M:%S')))
    py_res <- py_mod$run_graphst(
        count_matrices = inputs$counts,
        coord_matrices = inputs$coords,
        barcodes = inputs$barcodes,
        genes = inputs$genes,
        sample_ids = inputs$sample_names,
        n_clusters = n_clusters,
        seed = as.integer(SEED),
        device = device,
        cluster_method = cluster_method,
        refinement = refinement,
        radius = as.integer(radius),
        datatype = datatype,
        epochs = as.integer(epochs),
        start = start,
        end = end,
        increment = increment
    )
    show(paste("GraphST clustering completed. Time", format(Sys.time(), '%H:%M:%S')))
    .formatPythonBenchmarkResult(py_res, colnames(spe))
}


runSpaGCN <- function(spe, sample_label = NULL, samples = NULL, n_clusters, SEED,
                      assay_name = "counts", images = NULL,
                      use_histology = FALSE, refine = FALSE,
                      refine_shape = "square", alpha = 1, beta = 49, p = 0.5,
                      l_value = NULL, res = NULL, init_spa = TRUE,
                      init = "louvain", tol = 5e-3, lr = 0.05,
                      max_epochs = 200, search_res_start = 0.7,
                      search_res_step = 0.1, search_res_tol = 5e-3,
                      search_res_epochs = 20, l_search_start = 0.01,
                      l_search_end = 1000, l_search_tol = 0.01,
                      l_search_max_run = 100, min_cells = 3,
                      py_file = NULL, python_path = NULL) {
    
    sample_label <- .resolveSampleLabel(sample_label, samples)
    inputs <- .preparePythonBenchmarkInputs(spe, sample_label, assay_name)
    n_clusters <- .expandBySample(n_clusters, inputs$sample_names, "n_clusters")
    refine_shape <- .expandBySample(refine_shape, inputs$sample_names, "refine_shape")
    if (!is.null(images)) {
        images <- .expandBySample(images, inputs$sample_names, "images")
    }
    if (!is.null(l_value)) {
        l_value <- .expandBySample(l_value, inputs$sample_names, "l_value")
    }
    if (!is.null(res)) {
        res <- .expandBySample(res, inputs$sample_names, "res")
    }
    py_mod <- .getBenchmarkPythonModule(py_file, python_path)
    
    show(paste("SpaGCN run started. Time", format(Sys.time(), '%H:%M:%S')))
    py_res <- py_mod$run_spagcn(
        count_matrices = inputs$counts,
        coord_matrices = inputs$coords,
        barcodes = inputs$barcodes,
        genes = inputs$genes,
        sample_ids = inputs$sample_names,
        n_clusters = n_clusters,
        seed = as.integer(SEED),
        images = images,
        use_histology = use_histology,
        refine = refine,
        refine_shape = refine_shape,
        alpha = alpha,
        beta = beta,
        p = p,
        l_value = l_value,
        res = res,
        init_spa = init_spa,
        init = init,
        tol = tol,
        lr = lr,
        max_epochs = as.integer(max_epochs),
        search_res_start = search_res_start,
        search_res_step = search_res_step,
        search_res_tol = search_res_tol,
        search_res_epochs = as.integer(search_res_epochs),
        l_search_start = l_search_start,
        l_search_end = l_search_end,
        l_search_tol = l_search_tol,
        l_search_max_run = as.integer(l_search_max_run),
        min_cells = as.integer(min_cells)
    )
    show(paste("SpaGCN clustering completed. Time", format(Sys.time(), '%H:%M:%S')))
    .formatPythonBenchmarkResult(py_res, colnames(spe))
}


runSTAGATE <- function(spe, sample_label = NULL, samples = NULL, n_clusters, SEED,
                       assay_name = "counts", device = "cpu",
                       rad_cutoff = 150, cluster_method = "mclust",
                       n_top_genes = 3000, resolution = 1,
                       use_rep = "STAGATE", target_sum = 1e4,
                       n_epochs = NULL, py_file = NULL,
                       python_path = NULL) {
    
    sample_label <- .resolveSampleLabel(sample_label, samples)
    inputs <- .preparePythonBenchmarkInputs(spe, sample_label, assay_name)
    n_clusters <- .expandBySample(n_clusters, inputs$sample_names, "n_clusters")
    rad_cutoff <- .expandBySample(rad_cutoff, inputs$sample_names, "rad_cutoff")
    py_mod <- .getBenchmarkPythonModule(py_file, python_path)
    
    show(paste("STAGATE run started. Time", format(Sys.time(), '%H:%M:%S')))
    py_res <- py_mod$run_stagate(
        count_matrices = inputs$counts,
        coord_matrices = inputs$coords,
        barcodes = inputs$barcodes,
        genes = inputs$genes,
        sample_ids = inputs$sample_names,
        n_clusters = n_clusters,
        seed = as.integer(SEED),
        device = device,
        rad_cutoff = rad_cutoff,
        cluster_method = cluster_method,
        n_top_genes = as.integer(n_top_genes),
        resolution = resolution,
        use_rep = use_rep,
        target_sum = target_sum,
        n_epochs = n_epochs
    )
    show(paste("STAGATE clustering completed. Time", format(Sys.time(), '%H:%M:%S')))
    .formatPythonBenchmarkResult(py_res, colnames(spe))
}
