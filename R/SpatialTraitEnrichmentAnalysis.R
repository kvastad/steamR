#' Spatial Trait Enrichment Analysis
#'
#' Performs permutation-based enrichment analysis of a trait-specific gene score 
#' (e.g. from OpenTargets) across spatial or single-cell clusters, comparing 
#' observed medians to a null distribution.
#'
#' @param se A Seurat object containing gene scores and clustering metadata.
#' @param perm.mat.label.data A data frame of null median scores for the full gene list,
#'        where columns correspond to cluster names.
#' @param perm.mat.window50.data A data frame of null median scores for ranked gene sets
#'        (e.g., sliding windows), with columns corresponding to cluster names.
#' @param window_rank_list_abr_label A list or vector of window rank identifiers 
#'        used to reference the rank-specific gene sets.
#' @param gene_list A character string matching the name of the gene score column in `se@meta.data`
#'        to test for enrichment (e.g., "OpenTargets_SCZ_Genetic_1").
#' @param cluster_anno Column name in `se@meta.data` specifying the clustering to use 
#'        (e.g., "seurat_clusters" or "supercluster_term").
#' @param imputation Strategy for handling p-value calculation when no permutations exceed the observed value.
#'        Options are:
#'        - "all": Always add 1 to numerator and denominator (default)
#'        - "none": No imputation, can result in p=0
#'        - "dynamic": Only impute when no permutations exceed observed value
#' @param log_file Path to the log file for imputed p-values
#'
#' @returns A data frame of unadjusted p-values for each cluster.
#'
#' @export
#'
#' @examples
#' pval_mat <- SpatialTraitEnrichmentAnalysis(
#'   se = se,
#'   perm.mat.label.data = perm.mat.genetic.data,
#'   perm.mat.window50.data = perm.mat.window50.data,
#'   window_rank_list_abr_label = window50_rank_list_SCZ_Genetic,
#'   gene_list = "OpenTargets_SCZ_Genetic_1",
#'   cluster_anno = "supercluster_term",
#'   imputation = "dynamic"
#' )
SpatialTraitEnrichmentAnalysis <- function(
    se,
    perm.mat.label.data,
    perm.mat.window50.data,
    window_rank_list_abr_label,
    gene_list,
    cluster_anno = "seurat_clusters",
    imputation = "all",
    log_file = NULL
) {
  # Validate imputation parameter
  if (!imputation %in% c("all", "none", "dynamic")) {
    stop("imputation must be one of: 'all', 'none', 'dynamic'")
  }
  
  # Check that clustering column exists
  if (!(cluster_anno %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_anno, "is not found in the Seurat object metadata."))
  }
  
  # Check gene list column exists
  if (!(gene_list %in% colnames(se@meta.data))) {
    stop(paste("The specified gene list column", gene_list, "is not found in the metadata."))
  }
  
  cluster_numbers <- sort(unique(se@meta.data[[cluster_anno]]))
  cluster_names <- paste0("cluster_", cluster_numbers)
  
  # Initialize results matrix
  permutation_nr <- dim(perm.mat.label.data)[1]
  nr_of_tests <- length(cluster_numbers)
  p_val_mat <- matrix(NA, nrow = 1, ncol = length(cluster_names),
                      dimnames = list("p.val", cluster_names))
  
  # Open a connection to the log file if it's provided
  if (!is.null(log_file)) {
    log_con <- file(log_file, open = "wt")
  }
  
  # Process each cluster
  for (index in seq_along(cluster_numbers)) {
    i <- cluster_numbers[index]
    cluster_name <- as.character(i)
    
    if (!cluster_name %in% colnames(perm.mat.label.data)) {
      warning(paste("Cluster", cluster_name, "not found in permutation matrix. Skipping."))
      next
    }
    
    # Get cells for this cluster
    cells_i <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] == i]
    
    # Calculate median for the gene list in this cluster
    actual_median <- median(se@meta.data[cells_i, gene_list], na.rm = TRUE)
    
    # Get matched columns for this cluster
    matched_cols <- grep(gene_list, colnames(se@meta.data), value = TRUE)
    
    if (length(matched_cols) == 0) {
      warning(paste("No matching columns found for:", gene_list, "in cluster:", i))
      next
    }
    
    # Calculate medians for matched columns
    cluster_medians <- sapply(matched_cols, function(col) {
      median(se@meta.data[cells_i, col], na.rm = TRUE)
    })
    
    # Process windows
    for (j in seq_along(window_rank_list_abr_label)) {
      OT_label_window <- cluster_medians[j]
      
      perm_mat_window_cluster_bigger <- subset(perm.mat.window50.data,
                                               perm.mat.window50.data[[cluster_name]] >= OT_label_window)
      
      # Calculate p-value with specified imputation strategy
      n_bigger <- nrow(perm_mat_window_cluster_bigger)
      if (n_bigger == 0) {
        warning_msg <- sprintf("No median scores in the null distribution were equal to or larger than the queried median score for cluster %s. Consider using imputation 'all' or 'dynamic'.", cluster_name)
        if (!is.null(log_file)) {
          writeLines(warning_msg, log_con)
        } else {
          warning(warning_msg)
        }
        cluster_i_w_p_val <- 1 / (permutation_nr + 1)
      } else {
        cluster_i_w_p_val <- n_bigger / permutation_nr
      }
      
      perm_mat_label_bigger <- subset(perm.mat.label.data,
                                      perm.mat.label.data[[cluster_name]] >= actual_median)
      
      # Calculate p-value with specified imputation strategy
      n_bigger <- nrow(perm_mat_label_bigger)
      if (n_bigger == 0) {
        warning_msg <- sprintf("No median scores in the null distribution were equal to or larger than the queried median score for cluster %s. Consider using imputation 'all' or 'dynamic'.", cluster_name)
        if (!is.null(log_file)) {
          writeLines(warning_msg, log_con)
        } else {
          warning(warning_msg)
        }
        cluster_i_p_val <- 1 / (permutation_nr + 1)
      } else {
        cluster_i_p_val <- n_bigger / permutation_nr
      }
      
      p_val_mat[1, index] <- cluster_i_p_val
    }
  }
  
  # Close the connection to the log file
  if (!is.null(log_file)) {
    close(log_con)
  }
  
  return(as.data.frame(p_val_mat))
}