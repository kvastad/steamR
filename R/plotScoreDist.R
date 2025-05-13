#' Plot Enrichment Score Distributions for Clusters
#'
#' Generates density plots of null distribution scores for each cluster
#' and overlays the actual enrichment score to visually compare significance.
#' This function is useful for assessing whether observed enrichment scores
#' deviate from a random background.
#'
#' @param se A Seurat object containing enrichment scores in the metadata.
#' @param perm.mat A null distribution permutation matrix of enrichment scores per cluster (global medians). 
#' The permutation matrix is generated using the generatePermutationMatrix function and must be generated 
#' from a gene list of equal size as the gene list score given in 'enrichment_score_col'.
#' @param cluster_anno The column in the Seurat metadata that contains cluster labels (default is "seurat_clusters").
#' @param clusters_to_plot Optional vector of indices to subset and plot specific clusters only.
#' @param enrichment_score_col A string specifying the enrichment score column to use.
#'
#' @returns No return value. Produces a series of ggplot2 density plots, one per cluster.
#' @export
#'
#' @examples
#' plotScoreDist(
#'   se = se,
#'   perm.mat = perm.mat.genetic.data,
#'   cluster_anno = "supercluster_term",
#'   enrichment_score_col = "OpenTargets_SCZ_Genetic_1"
#' )
plotScoreDist <- function(se, 
                          perm.mat, 
                          cluster_anno = "seurat_clusters",
                          clusters_to_plot = NULL,
                          enrichment_score_col) {
  
  theme_set(theme_classic())
  
  cluster_labels <- unique(se@meta.data[[cluster_anno]])
  
  # Separate numeric and non-numeric IDs
  suppressWarnings({
    numeric_ids <- cluster_labels[!is.na(as.numeric(as.character(cluster_labels)))]
    non_numeric_ids <- cluster_labels[is.na(as.numeric(as.character(cluster_labels)))]
  })
  
  # Sort numerically and alphabetically
  sorted_ids <- c(
    sort(as.numeric(as.character(numeric_ids))),
    sort(as.character(non_numeric_ids))
  )
  
  cluster_labels <- sorted_ids
  
  # Subset to specified clusters if provided
  if (!is.null(clusters_to_plot)) {
    cluster_labels <- cluster_labels[clusters_to_plot]
  }
  
  # Find all relevant OT enrichment columns once
  start_with_abr_label <- grep(enrichment_score_col, colnames(se@meta.data), value = TRUE)
  if (length(start_with_abr_label) == 0) {
    stop("No matching columns found for term")
  }
  
  for (cluster_label in cluster_labels) {
    message("\nProcessing cluster: ", cluster_label)
    
    # Get cells for this cluster
    cells_in_cluster <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] == cluster_label]
    if (length(cells_in_cluster) == 0) {
      warning("No cells found for cluster: ", cluster_label)
      next
    }
    
    # Calculate medians for all relevant columns for this cluster
    meta_data_cluster <- se@meta.data[cells_in_cluster, start_with_abr_label, drop = FALSE]
    cluster_medians <- apply(meta_data_cluster, 2, median, na.rm = TRUE)
    
    # Use plain cluster name
    cluster_name <- as.character(cluster_label)
    
    if (!(cluster_name %in% colnames(perm.mat))) {
      warning("Cluster name not found in permutation matrix: ", cluster_name)
      next
    }
    
    # Plot global distribution
    median_score_null <- data.frame(Median_scores = perm.mat[[cluster_name]])
    
    if (!(enrichment_score_col %in% colnames(se@meta.data))) {
      warning("Missing OT score column: ", enrichment_score_col, " in cluster: ", cluster_label)
      next
    }
    
    actual_median <- median(se@meta.data[cells_in_cluster, enrichment_score_col], na.rm = TRUE)
    
    p <- ggplot(median_score_null, aes(x = Median_scores)) +
      geom_density(fill = "grey", alpha = 0.8) +
      geom_vline(xintercept = actual_median, color = "red", size = 1) +
      ggtitle(paste("Cluster:", cluster_label)) +
      xlab("Median Score") +
      ylab("Density")
    
    print(p)
  }
}
