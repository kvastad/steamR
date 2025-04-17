#' Generate a permutation matrix of module scores
#'
#' This function computes a permutation matrix by generating random gene sets and scoring them
#' using Seurat’s `AddModuleScore()`. Results are aggregated as median scores per cluster.
#' Supports both sequential and parallel computation for reproducibility and speed.
#'
#' @param se A Seurat object containing gene expression and metadata.
#' @param gnum Number of genes to sample for each permutation (default: 50).
#' @param permutation_nr Number of permutations to perform (default: 1000).
#' @param cluster_col The name of the column in the metadata with cluster identities (default: `"seurat_clusters"`).
#' @param workers Optional. Number of parallel workers. If `NULL` or `workers <= 1`, runs sequentially.
#'
#' @return A data frame where rows represent permutations and columns represent clusters.
#'         Each entry is the median module score of a randomly sampled gene set for that cluster.
#' @export
#'
#' @examples
#' # Run sequentially
#' perm.mat <- generate_permutation_matrix(se, gnum = 50, permutation_nr = 100)
#'
#' # Run in parallel
#' perm.mat <- generate_permutation_matrix(se, gnum = 50, permutation_nr = 100, workers = 4)
generate_permutation_matrix <- function(se, 
                                        gnum = 50, 
                                        permutation_nr = 1000, 
                                        cluster_col = "seurat_clusters", 
                                        workers = NULL) {
  parallel <- !is.null(workers) && workers > 1
  
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  cluster_names <- unique(se@meta.data[[cluster_col]])
  
  suppressWarnings({
    numeric_ids <- cluster_names[!is.na(as.numeric(as.character(cluster_names)))]
    non_numeric_ids <- cluster_names[is.na(as.numeric(as.character(cluster_names)))]
  })
  
  sorted_names <- c(sort(as.numeric(as.character(numeric_ids))),
                    sort(as.character(non_numeric_ids)))
  cluster_names <- sorted_names
  Perm_names <- paste0("Perm_nr", seq_len(permutation_nr))
  ref_genes_se <- rownames(se)
  
  gene_sets <- lapply(seq_len(permutation_nr), function(i) {
    set.seed(i)
    sample(ref_genes_se, gnum)
  })
  
  execution_time <- system.time({
    score_fun <- function(i) {
      genes <- gene_sets[[i]]
      se_temp <- AddModuleScore(se, list(genes), ctrl = 100, name = "NULL_gene_set")
      df <- se_temp@meta.data[, c(cluster_col, "NULL_gene_set1")]
      colnames(df)[1] <- "cluster"
      df %>%
        group_by(cluster) %>%
        summarise(median_score = median(NULL_gene_set1, na.rm = TRUE), .groups = "drop") %>%
        `[[`("median_score")
    }
    
    if (parallel) {
      options(future.globals.maxSize = 4 * 1024^3)
      plan(multisession, workers = workers)
      perm_results <- future.apply::future_lapply(seq_len(permutation_nr), score_fun, future.seed = TRUE)
    } else {
      perm_results <- lapply(seq_len(permutation_nr), score_fun)
    }
    
    perm.mat <- do.call(rbind, perm_results)
    rownames(perm.mat) <- Perm_names
    colnames(perm.mat) <- cluster_names
  })
  
  print(execution_time)
  as.data.frame(perm.mat)
}
