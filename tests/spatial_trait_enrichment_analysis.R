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
#' @param cluster_col Column name in `se@meta.data` specifying the clustering to use 
#'        (e.g., "seurat_clusters" or "supercluster_term").
#'
#' @returns A matrix of enrichment results for each cluster. Rows contain:
#' \describe{
#'   \item{-log10(p.val)}{The log-transformed unadjusted p-values}
#'   \item{-log10(p.val.adj)}{The log-transformed Bonferroni-adjusted p-values}
#' }
#' Columns correspond to each cluster analyzed.
#'
#' @export
#'
#' @examples
#' pval_mat <- spatial_trait_enrichment_analysis(
#'   se = se,
#'   perm.mat.label.data = perm.mat.genetic.data,
#'   perm.mat.window50.data = perm.mat.window50.data,
#'   window_rank_list_abr_label = window50_rank_list_SCZ_Genetic,
#'   gene_list = "OpenTargets_SCZ_Genetic_1",
#'   cluster_col = "supercluster_term"
#' )
spatial_trait_enrichment_analysis <- function(
    se,
    perm.mat.label.data,
    perm.mat.window50.data,
    window_rank_list_abr_label,
    gene_list,
    cluster_col = "seurat_clusters"
) {
  se_subset <- list()
  meta_data_subset <- list()
  meta_data_abr_label <- list()
  
  # Check that clustering column exists
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  # Check gene list column exists
  if (!(gene_list %in% colnames(se@meta.data))) {
    stop(paste("The specified gene list column", gene_list, "is not found in the metadata."))
  }
  
  cluster_numbers <- sort(unique(se@meta.data[[cluster_col]]))
  
  for (i in cluster_numbers) {
    subset_name <- paste0("se_", i)
    
    cells_i <- rownames(se@meta.data)[se@meta.data[[cluster_col]] == i]
    se_subset[[subset_name]] <- subset(se, cells = cells_i)
    
    meta_data_subset[[subset_name]] <- se_subset[[subset_name]][[]]
    matched_cols <- grep(gene_list, colnames(meta_data_subset[[subset_name]]), value = TRUE)
    
    if (length(matched_cols) == 0) {
      warning(paste("No matching columns found for:", gene_list, "in cluster:", i))
      next
    }
    
    meta_data_abr_label[[subset_name]] <- subset(meta_data_subset[[subset_name]], select = matched_cols)
  }
  
  permutation_nr <- dim(perm.mat.label.data)[1]
  nr_of_tests <- length(cluster_numbers)
  cluster_names <- paste0("cluster_", cluster_numbers)
  p_val_mat <- matrix(NA, nrow = 2, ncol = length(cluster_names),
                      dimnames = list(c("-log10(p.val)", "-log10(p.val.adj)"), cluster_names))
  
  for (index in seq_along(cluster_numbers)) {
    i <- cluster_numbers[index]
    subset_name <- paste0("se_", i)
    cluster_name <- as.character(i)
    
    if (!cluster_name %in% colnames(perm.mat.label.data)) {
      warning(paste("Cluster", cluster_name, "not found in permutation matrix. Skipping."))
      next
    }
    
    abr_label_structure_i <- data.frame()
    for (k in 1:ncol(meta_data_abr_label[[subset_name]])) {
      median_value <- median(meta_data_abr_label[[subset_name]][, k], na.rm = TRUE)
      abr_label_structure_i <- rbind(abr_label_structure_i, data.frame(Rank = as.character(k), Median = median_value))
    }
    
    for (j in seq_along(window_rank_list_abr_label)) {
      OT_label_window <- abr_label_structure_i[j, "Median"]
      
      perm_mat_window_cluster_bigger <- subset(perm.mat.window50.data,
                                               perm.mat.window50.data[[cluster_name]] >= OT_label_window)
      
      cluster_i_w_p_val <- (nrow(perm_mat_window_cluster_bigger) + 1) / (permutation_nr + 1)
      adjusted_p_val <- p.adjust(cluster_i_w_p_val, method = "bonferroni", n = nr_of_tests)
      
      actual_median <- median(as.numeric(unlist(se_subset[[subset_name]][[gene_list]])), na.rm = TRUE)
      perm_mat_label_bigger <- subset(perm.mat.label.data,
                                      perm.mat.label.data[[cluster_name]] >= actual_median)
      
      cluster_i_p_val <- (nrow(perm_mat_label_bigger) + 1) / (permutation_nr + 1)
      cluster_i_p_val_adj <- p.adjust(cluster_i_p_val, method = "bonferroni", n = nr_of_tests)
      
      p_val_mat[1, index] <- -log10(cluster_i_p_val)
      p_val_mat[2, index] <- -log10(cluster_i_p_val_adj)
    }
  }
  
  return(p_val_mat)
}