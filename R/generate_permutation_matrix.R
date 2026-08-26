#' Generate a permutation matrix of module scores
#'
#' This function computes a permutation matrix by generating random gene sets
#' and scoring them using Seurat's `AddModuleScore()`. Results are aggregated
#' as median scores per cluster. Supports both sequential and parallel
#' computation, with reproducible random gene-set generation.
#'
#' @param se A Seurat object containing gene expression and metadata.
#' @param gnum Number of genes to sample for each permutation (default: 50).
#' @param permutation_nr Number of permutations to perform (default: 1000).
#' @param cluster_anno The name of the column in the metadata with cluster
#'   identities (default: `"seurat_clusters"`).
#' @param workers Optional. Number of parallel workers. If `NULL` or
#'   `workers <= 1`, runs sequentially.
#' @param ctrl_nr Number of control genes used by `AddModuleScore()`
#'   for each gene in the input gene set (default: 100).
#' @param nbin_nr Number of expression-level bins used by
#'   `AddModuleScore()` (default: 24).
#'
#' @return A data frame where rows represent permutations and columns
#'   represent clusters. Each entry is the median module score of a randomly
#'   sampled gene set for that cluster.
#' @export
#'
#' @examples
#' # Run sequentially
#' perm.mat <- generate_permutation_matrix(
#'   se,
#'   gnum = 50,
#'   permutation_nr = 100
#' )
#'
#' # Run in parallel
#' perm.mat <- generate_permutation_matrix(
#'   se,
#'   gnum = 50,
#'   permutation_nr = 100,
#'   workers = 4
#' )
generate_permutation_matrix <- function(
    se,
    gnum = 50,
    permutation_nr = 1000,
    cluster_anno = "seurat_clusters",
    workers = NULL,
    ctrl_nr = 100,
    nbin_nr = 24
) {
  
  parallel <- !is.null(workers) && workers > 1
  
  if (!(cluster_anno %in% colnames(se@meta.data))) {
    stop(
      paste(
        "The specified cluster column",
        cluster_anno,
        "is not found in the Seurat object metadata."
      )
    )
  }
  
  cluster_names <- sort_cluster_names(
    unique(se@meta.data[[cluster_anno]])
  )
  
  Perm_names <- paste0("Perm_nr", seq_len(permutation_nr))
  
  gene_sets <- generate_random_gene_sets(
    genes = rownames(se),
    gnum = gnum,
    permutation_nr = permutation_nr
  )
  
  execution_time <- system.time({
    
    if (parallel) {
      options(future.globals.maxSize = 4 * 1024^3)
      plan(multisession, workers = workers)
      
      perm_results <- future.apply::future_lapply(
        gene_sets,
        calculate_permutation_score,
        se = se,
        cluster_anno = cluster_anno,
        ctrl_nr = ctrl_nr,
        nbin_nr = nbin_nr,
        future.seed = TRUE
      )
      
    } else {
      perm_results <- lapply(
        gene_sets,
        calculate_permutation_score,
        se = se,
        cluster_anno = cluster_anno,
        ctrl_nr = ctrl_nr,
        nbin_nr = nbin_nr
      )
    }
    
    perm.mat <- do.call(rbind, perm_results)
    rownames(perm.mat) <- Perm_names
    colnames(perm.mat) <- cluster_names
  })
  
  print(execution_time)
  
  as.data.frame(perm.mat)
}


#' Generate random gene sets for permutation analysis
#'
#' Generates a list of randomly sampled gene sets to be used for permutation
#' analysis. A separate seed is used for each permutation to make the
#' generated gene sets reproducible.
#'
#' @param genes A character vector containing the genes from which to sample.
#' @param gnum Number of genes to sample for each permutation.
#' @param permutation_nr Number of random gene sets to generate.
#'
#' @return A list containing one randomly sampled gene set per permutation.
generate_random_gene_sets <- function(
    genes,
    gnum,
    permutation_nr
) {
  
  lapply(seq_len(permutation_nr), function(i) {
    set.seed(i)
    sample(genes, gnum)
  })
}


#' Calculate median module scores for a random gene set
#'
#' Calculates Seurat module scores for a randomly sampled gene set and
#' summarises the scores as the median module score for each cluster.
#'
#' @param se A Seurat object containing gene expression and metadata.
#' @param genes A character vector containing the genes to score.
#' @param cluster_anno The name of the metadata column containing cluster
#'   identities.
#' @param ctrl_nr Number of control genes used by `AddModuleScore()`
#'   for each gene in the input gene set.
#' @param nbin_nr Number of expression-level bins used by
#'   `AddModuleScore()`.
#'
#' @return A numeric vector containing the median module score for each
#'   cluster.
calculate_permutation_score <- function(
    se,
    genes,
    cluster_anno,
    ctrl_nr,
    nbin_nr
) {
  
  se_temp <- AddModuleScore(
    se,
    list(genes),
    ctrl = ctrl_nr,
    nbin = nbin_nr,
    name = "NULL_gene_set"
  )
  
  df <- se_temp@meta.data[, c(cluster_anno, "NULL_gene_set1")]
  colnames(df)[1] <- "cluster"
  
  df %>%
    group_by(cluster) %>%
    summarise(
      median_score = median(NULL_gene_set1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    `[[`("median_score")
}


#' Sort cluster names
#'
#' Sorts cluster identifiers by placing numeric cluster identifiers first
#' in numerical order, followed by non-numeric cluster identifiers in
#' alphabetical order.
#'
#' @param cluster_names A vector containing cluster identifiers.
#'
#' @return A vector of sorted cluster identifiers.
sort_cluster_names <- function(cluster_names) {
  
  suppressWarnings({
    numeric_ids <- cluster_names[
      !is.na(as.numeric(as.character(cluster_names)))
    ]
    
    non_numeric_ids <- cluster_names[
      is.na(as.numeric(as.character(cluster_names)))
    ]
  })
  
  c(
    sort(as.numeric(as.character(numeric_ids))),
    sort(as.character(non_numeric_ids))
  )
}
