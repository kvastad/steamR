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
  
  # Check if the specified cluster column exists
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  execution_time <- system.time({
    
    ref_genes_se <- rownames(se)
    genes.in.NULL.set <- gnum
    
    # Use the specified cluster column
    cluster_names <- unique(se@meta.data[[cluster_col]])
    Perm_names <- paste0("Perm_nr", seq_len(permutation_nr))
    perm.mat <- matrix(data = NA, nrow = permutation_nr, ncol = length(cluster_names), dimnames = list(Perm_names, cluster_names))
    
    perm_results <- future_lapply(seq_len(permutation_nr), function(i) {
      NULL_gene_set_ <- sample(ref_genes_se, genes.in.NULL.set)
      genes <- unlist(list(NULL_gene_set_))
      
      se_temp <- AddModuleScore(se, list(genes), ctrl = 100, name = "NULL_gene_set")
      se_metadata_NULL <- se_temp@meta.data[, c(cluster_col, "NULL_gene_set1")]
      
      # Rename the cluster column for consistent processing
      colnames(se_metadata_NULL)[1] <- "cluster"
      
      se_metadata_by_cluster_NULL <- se_metadata_NULL %>%
        group_by(cluster) %>%
        summarise(median_score = median(NULL_gene_set1, na.rm = TRUE))
      
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
