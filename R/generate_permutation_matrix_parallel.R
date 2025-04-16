#' Generate a permutation matrix
#'
#' @param se A Seurat object
#' @param gnum Number of genes to sample for each permutation
#' @param permutation_nr Number of permutations
#' @param workers Number of parallel workers
#' @param cluster_col The name of the column in the metadata that contains cluster information
#'
#' @returns A data frame containing the permutation matrix
#' @export
#'
#' @examples
#' perm.mat <- generate_permutation_matrix_parallel(se, gnum = 50, permutation_nr = 100, workers = 4)
generate_permutation_matrix_parallel <- function(se, gnum = 50, permutation_nr = 10000, workers = 10, cluster_col = "seurat_clusters") {
  
  options(future.globals.maxSize = 4 * 1024^3)
  plan(multisession, workers = workers)
  
  # Check cluster column
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  execution_time <- system.time({
    
    ref_genes_se <- rownames(se)
    cluster_names <- unique(se@meta.data[[cluster_col]])
    Perm_names <- paste0("Perm_nr", seq_len(permutation_nr))
    
    # Pre-generate gene sets for each permutation
    set.seed(1)
    gene_sets <- replicate(permutation_nr, sample(ref_genes_se, gnum), simplify = FALSE)
    
    # Run AddModuleScore on each gene set
    perm_results <- future_lapply(seq_len(permutation_nr), function(i) {
      genes <- gene_sets[[i]]
      
      se_temp <- AddModuleScore(se, list(genes), ctrl = 100, name = "NULL_gene_set")
      se_metadata_NULL <- se_temp@meta.data[, c(cluster_col, "NULL_gene_set1")]
      
      colnames(se_metadata_NULL)[1] <- "cluster"
      
      se_metadata_by_cluster_NULL <- se_metadata_NULL %>%
        group_by(cluster) %>%
        summarise(median_score = median(NULL_gene_set1, na.rm = TRUE), .groups = "drop")
      
      return(se_metadata_by_cluster_NULL$median_score)
    }, future.seed = TRUE)
    
    perm.mat <- do.call(rbind, perm_results)
    rownames(perm.mat) <- Perm_names
    colnames(perm.mat) <- cluster_names
    perm.mat <- as.data.frame(perm.mat)
  })
  
  print(execution_time)
  return(perm.mat)
}
