#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0 || is.na(x[[1]])) {
        return(y)
    }
    x
}


getScriptPath <- function() {
    file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
    if (length(file_arg) == 0) {
        stop("Unable to determine script path. Run this file with `Rscript`.")
    }
    normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE)
}


attachPackages <- function(pkgs) {
    for (pkg in pkgs) {
        suppressPackageStartupMessages(
            library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        )
    }
}


missingPackages <- function(pkgs) {
    pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
}


missingPythonModules <- function(modules, python_path = NULL) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        return(modules)
    }
    if (!is.null(python_path) && nzchar(python_path)) {
        Sys.setenv(RETICULATE_PYTHON = python_path)
        reticulate::use_python(python_path, required = TRUE)
    }
    modules[!vapply(modules, reticulate::py_module_available, FUN.VALUE = logical(1))]
}


downloadInputObject <- function(destfile, url) {
    dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
    utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
    destfile
}


loadInputObject <- function(input_path = NULL,
                            input_url = "https://raw.githubusercontent.com/SydneyBioX/clustSIGNAL/main/data/mEmbryo2.RData") {
    cache_dir <- file.path(tempdir(), "benchmark_smoke_test")
    local_file <- input_path %||% file.path(cache_dir, basename(input_url))
    if (!file.exists(local_file)) {
        message("Downloading input object to ", local_file)
        downloadInputObject(local_file, input_url)
    }
    
    e <- new.env(parent = emptyenv())
    object_names <- load(local_file, envir = e)
    if (length(object_names) == 0) {
        stop("No objects were found in ", local_file)
    }
    
    spe_name <- object_names[[1]]
    for (nm in object_names) {
        obj <- get(nm, envir = e)
        if (inherits(obj, "SpatialExperiment")) {
            spe_name <- nm
            break
        }
    }
    obj <- get(spe_name, envir = e)
    message("Loaded object `", spe_name, "` from ", local_file)
    if (inherits(obj, "SpatialExperiment")) {
        return(obj)
    }
    
    object_names <- ls(e)
    if ("me_expr" %in% object_names && "me_data" %in% object_names) {
        counts_mat <- e[["me_expr"]]
        meta_df <- e[["me_data"]]
        if (!inherits(counts_mat, "Matrix") && !is.matrix(counts_mat)) {
            stop("`me_expr` is not a matrix-like counts object.")
        }
        if (!is.data.frame(meta_df)) {
            stop("`me_data` is not a data.frame.")
        }
        if (!all(c("X", "Y") %in% colnames(meta_df))) {
            stop("`me_data` must contain `X` and `Y` coordinate columns.")
        }
        if (!identical(rownames(meta_df), colnames(counts_mat))) {
            stop("`me_data` rownames do not match count matrix column names.")
        }
        
        message("Constructing SpatialExperiment from `me_expr` + `me_data`.")
        return(SpatialExperiment::SpatialExperiment(
            assays = list(counts = counts_mat),
            colData = S4Vectors::DataFrame(meta_df),
            spatialCoords = as.matrix(meta_df[, c("X", "Y"), drop = FALSE])
        ))
    }
    
    obj
}


detectSampleColumn <- function(spe) {
    candidates <- c(
        "embryo", "sample_id", "Run_Tissue_name", "samples", "sample",
        "section", "slice", "batch", "orig.ident"
    )
    present <- intersect(candidates, colnames(SummarizedExperiment::colData(spe)))
    if (length(present) > 0) {
        return(present[[1]])
    }
    
    variable_cols <- colnames(SummarizedExperiment::colData(spe))[vapply(
        as.data.frame(SummarizedExperiment::colData(spe)),
        function(x) length(unique(as.character(x))) > 1,
        FUN.VALUE = logical(1)
    )]
    if (length(variable_cols) > 0) {
        return(variable_cols[[1]])
    }
    
    stop("Could not detect a suitable sample column in colData(spe).")
}


detectAnnotationColumn <- function(spe) {
    candidates <- c(
        "celltype_mapped_refined", "celltype", "cell_type", "annotation",
        "annot", "labels", "cluster", "celltype_mapped"
    )
    present <- intersect(candidates, colnames(SummarizedExperiment::colData(spe)))
    if (length(present) > 0) {
        return(present[[1]])
    }
    
    non_constant_cols <- colnames(SummarizedExperiment::colData(spe))[vapply(
        as.data.frame(SummarizedExperiment::colData(spe)),
        function(x) length(unique(as.character(x))) > 1,
        FUN.VALUE = logical(1)
    )]
    if (length(non_constant_cols) > 0) {
        return(non_constant_cols[[1]])
    }
    
    stop("Could not detect a suitable annotation column in colData(spe).")
}


chooseSubset <- function(spe, sample_col, assay_name, max_cells = 120L,
                         max_genes = 300L, seed = 1L) {
    sample_vec <- as.character(spe[[sample_col]])
    sample_tab <- sort(table(sample_vec), decreasing = TRUE)
    sample_id <- names(sample_tab)[[1]]
    spe_sub <- spe[, sample_vec == sample_id]
    
    set.seed(seed)
    n_keep_cells <- min(max_cells, ncol(spe_sub))
    cell_idx <- sort(sample(seq_len(ncol(spe_sub)), size = n_keep_cells))
    spe_sub <- spe_sub[, cell_idx]
    
    total_counts <- Matrix::rowSums(SummarizedExperiment::assay(spe_sub, assay_name))
    positive_genes <- which(total_counts > 0)
    if (length(positive_genes) == 0) {
        stop("No non-zero genes found in the selected subset.")
    }
    gene_order <- positive_genes[order(total_counts[positive_genes], decreasing = TRUE)]
    n_keep_genes <- min(max_genes, length(gene_order))
    spe_sub <- spe_sub[gene_order[seq_len(n_keep_genes)], ]
    
    message(
        "Using sample `", sample_id, "` with ", ncol(spe_sub), " cells and ",
        nrow(spe_sub), " genes for smoke tests."
    )
    spe_sub
}


assertNamedClusters <- function(x, expected_ids) {
    if (length(x) != length(expected_ids)) {
        stop("Expected ", length(expected_ids), " clusters, got ", length(x), ".")
    }
    if (is.null(names(x))) {
        stop("Cluster vector is missing names.")
    }
    if (!identical(sort(names(x)), sort(expected_ids))) {
        stop("Cluster vector names do not match the input cell IDs.")
    }
}


validateBANKSY <- function(res, spe_small) {
    if (!inherits(res, "SpatialExperiment")) {
        stop("BANKSY did not return a SpatialExperiment.")
    }
    clust_cols <- grep("^clust_", colnames(SummarizedExperiment::colData(res)), value = TRUE)
    if (length(clust_cols) == 0) {
        stop("BANKSY output is missing clustering columns.")
    }
    invisible(TRUE)
}


validateBASS <- function(res, spe_small) {
    if (!isS4(res) || !"results" %in% methods::slotNames(res)) {
        stop("BASS output does not contain a `results` slot.")
    }
    bass_clusters <- methods::slot(res, "results")$c
    if (is.null(bass_clusters)) {
        stop("BASS output is missing `results$c`.")
    }
    if (length(unlist(bass_clusters)) != ncol(spe_small)) {
        stop("BASS cluster assignments do not match the input cell count.")
    }
    invisible(TRUE)
}


validateSpatialPCA <- function(res, spe_small, sample_col) {
    if (!is.list(res) || length(res) != 1) {
        stop("SpatialPCA smoke test expected a single-sample list output.")
    }
    assertNamedClusters(res[[1]], colnames(spe_small))
    invisible(TRUE)
}


validatePythonWrapper <- function(res, spe_small) {
    if (!is.list(res) || is.null(res$cluster_vector)) {
        stop("Python wrapper output is missing `cluster_vector`.")
    }
    if (length(res$cluster_vector) != ncol(spe_small)) {
        stop("`cluster_vector` length does not match the input cell count.")
    }
    if (is.null(names(res$cluster_vector))) {
        stop("`cluster_vector` is missing cell names.")
    }
    if (!identical(sort(names(res$cluster_vector)), sort(colnames(spe_small)))) {
        stop("`cluster_vector` names do not match the input cell IDs.")
    }
    invisible(TRUE)
}


runSmokeTest <- function(name, expr_fn, validator,
                         required_r = character(),
                         attach_r = character(),
                         required_py = character(),
                         python_path = NULL) {
    result <- list(name = name, status = "PASS", message = NA_character_)
    
    missing_r <- missingPackages(required_r)
    if (length(missing_r) > 0) {
        result$status <- "SKIP"
        result$message <- paste("Missing R packages:", paste(missing_r, collapse = ", "))
        return(result)
    }
    
    if (length(attach_r) > 0) {
        attachPackages(attach_r)
    }
    
    if (length(required_py) > 0) {
        missing_py <- missingPythonModules(required_py, python_path = python_path)
        if (length(missing_py) > 0) {
            result$status <- "SKIP"
            result$message <- paste(
                "Missing Python modules:", paste(missing_py, collapse = ", ")
            )
            return(result)
        }
    }
    
    tryCatch({
        out <- expr_fn()
        validator(out)
        result$message <- "ok"
        result
    }, error = function(e) {
        result$status <- "FAIL"
        result$message <- conditionMessage(e)
        result
    })
}


script_path <- getScriptPath()
repo_root <- dirname(dirname(script_path))
functions_file <- file.path(repo_root, "Benchmarking_scripts", "benchmarking_functions.R")

attachPackages(c("Matrix", "SummarizedExperiment", "SpatialExperiment"))
source(functions_file)

input_path <- Sys.getenv("BENCHMARK_INPUT", unset = "")
input_url <- Sys.getenv(
    "BENCHMARK_INPUT_URL",
    unset = "https://raw.githubusercontent.com/SydneyBioX/clustSIGNAL/main/data/mEmbryo2.RData"
)
python_path <- Sys.getenv("BENCHMARK_PYTHON", unset = "")
python_path <- if (nzchar(python_path)) python_path else NULL
if (!is.null(python_path)) {
    Sys.setenv(RETICULATE_PYTHON = python_path)
}

spe <- loadInputObject(
    input_path = if (nzchar(input_path)) input_path else NULL,
    input_url = input_url
)

if (!inherits(spe, "SpatialExperiment")) {
    stop("Input object is not a SpatialExperiment.")
}
if (ncol(spatialCoords(spe)) < 2) {
    stop("Input object does not contain 2D spatial coordinates.")
}

sample_col <- detectSampleColumn(spe)
annotation_col <- detectAnnotationColumn(spe)
assay_name <- if ("counts" %in% SummarizedExperiment::assayNames(spe)) {
    "counts"
} else {
    SummarizedExperiment::assayNames(spe)[[1]]
}

spe_small <- chooseSubset(
    spe = spe,
    sample_col = sample_col,
    assay_name = assay_name,
    max_cells = 120L,
    max_genes = 300L,
    seed = 1L
)

n_labels <- length(unique(as.character(spe_small[[annotation_col]])))
n_clusters <- max(2L, min(4L, n_labels))
n_domains <- max(2L, min(3L, n_clusters))

cntm <- list(SummarizedExperiment::assay(spe_small, assay_name))
xym <- list(spatialCoords(spe_small))

method_names <- c(
    "runBANKSY", "runBASS", "runSpatialPCA", "runGraphST", "runSpaGCN",
    "runSTAGATE"
)
selected_method <- Sys.getenv("BENCHMARK_METHOD", unset = "")
isolate_methods <- identical(Sys.getenv("BENCHMARK_ISOLATE", unset = ""), "1")

if (isolate_methods && !nzchar(selected_method)) {
    output_dir <- Sys.getenv("BENCHMARK_OUTPUT_DIR", unset = tempdir())
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    isolated_dir <- file.path(output_dir, "isolated")
    dir.create(isolated_dir, recursive = TRUE, showWarnings = FALSE)

    child_results <- lapply(method_names, function(method_name) {
        method_dir <- file.path(isolated_dir, method_name)
        dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
        log_file <- file.path(method_dir, "smoke_test.log")
        child_env <- c(
            "BENCHMARK_ISOLATE=0",
            paste0("BENCHMARK_METHOD=", method_name),
            paste0("BENCHMARK_OUTPUT_DIR=", method_dir)
        )
        for (env_name in c(
            "BENCHMARK_INPUT", "BENCHMARK_INPUT_URL", "BENCHMARK_PYTHON",
            "R_MAKEVARS_USER", "MPLCONFIGDIR"
        )) {
            env_value <- Sys.getenv(env_name, unset = "")
            if (nzchar(env_value)) {
                child_env <- c(child_env, paste0(env_name, "=", env_value))
            }
        }

        status <- system2(
            "Rscript",
            args = shQuote(script_path),
            env = child_env,
            stdout = log_file,
            stderr = log_file
        )
        result_file <- file.path(method_dir, "smoke_test_results.rds")
        if (file.exists(result_file)) {
            result <- readRDS(result_file)$summary[1, , drop = FALSE]
            if (!identical(status, 0L) && result$status == "PASS") {
                result$status <- "FAIL"
                result$message <- paste("Subprocess exited with status", status)
            }
            return(result)
        }

        data.frame(
            method = method_name,
            status = "FAIL",
            message = paste(
                "Subprocess exited before saving results with status", status,
                "- see", log_file
            ),
            stringsAsFactors = FALSE
        )
    })

    summary_df <- do.call(rbind, child_results)
    utils::write.csv(
        summary_df,
        file = file.path(output_dir, "smoke_test_summary.csv"),
        row.names = FALSE
    )
    saveRDS(
        list(
            summary = summary_df,
            results = child_results,
            sample_col = sample_col,
            annotation_col = annotation_col,
            assay_name = assay_name,
            n_cells = ncol(spe_small),
            n_genes = nrow(spe_small),
            n_clusters = n_clusters,
            n_domains = n_domains
        ),
        file = file.path(output_dir, "smoke_test_results.rds")
    )
    print(summary_df, row.names = FALSE)
    message("Saved smoke test results to ", normalizePath(output_dir))
    if (any(summary_df$status == "FAIL")) {
        stop("One or more benchmark smoke tests failed.")
    }
    quit(status = 0L)
}

runSelectedSmokeTest <- function(name, ...) {
    if (nzchar(selected_method) && name != selected_method) {
        return(NULL)
    }
    runSmokeTest(name = name, ...)
}

results <- list(
    runSelectedSmokeTest(
        name = "runBANKSY",
        expr_fn = function() runBANKSY(
            spe_small,
            annots_label = annotation_col,
            sample_label = sample_col,
            SEED = 1L,
            k_geom = c(8, 16),
            lambda = 0.2,
            res = 0.4,
            npcs = 10,
            use_agf = TRUE,
            compute_agf = TRUE
        ),
        validator = function(res) validateBANKSY(res, spe_small),
        required_r = c("Banksy", "Seurat"),
        attach_r = c("Banksy", "Seurat")
    ),
    runSelectedSmokeTest(
        name = "runBASS",
        expr_fn = function() runBASS(cntm = cntm, xym = xym, C = n_clusters, R = n_domains),
        validator = function(res) validateBASS(res, spe_small),
        required_r = c("BASS"),
        attach_r = c("BASS")
    ),
    runSelectedSmokeTest(
        name = "runSpatialPCA",
        expr_fn = function() runSpatialPCA(
            spe_small,
            samples = sample_col,
            sparkv = "sparkx",
            ncores = 1
        ),
        validator = function(res) validateSpatialPCA(res, spe_small, sample_col),
        required_r = c("SpatialPCA", "bluster"),
        attach_r = c("SpatialPCA", "bluster")
    ),
    runSelectedSmokeTest(
        name = "runGraphST",
        expr_fn = function() runGraphST(
            spe_small,
            sample_label = sample_col,
            n_clusters = n_clusters,
            SEED = 1L,
            assay_name = assay_name,
            cluster_method = "leiden",
            refinement = FALSE,
            radius = 20,
            epochs = 10,
            device = "cpu",
            python_path = python_path
        ),
        validator = function(res) validatePythonWrapper(res, spe_small),
        required_r = c("reticulate"),
        required_py = c("anndata", "numpy", "pandas", "scipy", "torch", "GraphST")
    ),
    runSelectedSmokeTest(
        name = "runSpaGCN",
        expr_fn = function() runSpaGCN(
            spe_small,
            sample_label = sample_col,
            n_clusters = n_clusters,
            SEED = 1L,
            assay_name = assay_name,
            use_histology = FALSE,
            refine = FALSE,
            init = "kmeans",
            l_value = 0.1,
            res = 0.4,
            max_epochs = 20,
            search_res_epochs = 5,
            min_cells = 1,
            python_path = python_path
        ),
        validator = function(res) validatePythonWrapper(res, spe_small),
        required_r = c("reticulate"),
        required_py = c("anndata", "numpy", "pandas", "scanpy", "scipy", "SpaGCN")
    ),
    runSelectedSmokeTest(
        name = "runSTAGATE",
        expr_fn = function() runSTAGATE(
            spe_small,
            sample_label = sample_col,
            n_clusters = n_clusters,
            SEED = 1L,
            assay_name = assay_name,
            device = "cpu",
            rad_cutoff = 80,
            cluster_method = "leiden",
            n_top_genes = min(200L, nrow(spe_small)),
            resolution = 0.4,
            n_epochs = 10,
            python_path = python_path
        ),
        validator = function(res) validatePythonWrapper(res, spe_small),
        required_r = c("reticulate"),
        required_py = c("anndata", "numpy", "pandas", "scanpy", "scipy", "torch", "STAGATE_pyG")
    )
)
results <- Filter(Negate(is.null), results)

summary_df <- do.call(
    rbind,
    lapply(results, function(x) {
        data.frame(
            method = x$name,
            status = x$status,
            message = x$message,
            stringsAsFactors = FALSE
        )
    })
)

print(summary_df, row.names = FALSE)

output_dir <- Sys.getenv("BENCHMARK_OUTPUT_DIR", unset = "")
if (nzchar(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(
        summary_df,
        file = file.path(output_dir, "smoke_test_summary.csv"),
        row.names = FALSE
    )
    saveRDS(
        list(
            summary = summary_df,
            results = results,
            sample_col = sample_col,
            annotation_col = annotation_col,
            assay_name = assay_name,
            n_cells = ncol(spe_small),
            n_genes = nrow(spe_small),
            n_clusters = n_clusters,
            n_domains = n_domains
        ),
        file = file.path(output_dir, "smoke_test_results.rds")
    )
    message("Saved smoke test results to ", normalizePath(output_dir))
}

if (any(summary_df$status == "FAIL")) {
    stop("One or more benchmark smoke tests failed.")
}
