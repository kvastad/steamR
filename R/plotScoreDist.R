#' Plot Enrichment Score Distributions for Clusters
#'
#' Generates density plots of NULL distribution scores for each cluster
#' and overlays the actual enrichment score to visually compare significance.
#' This function is useful for assessing whether observed enrichment scores
#' deviate from a random background.
#'
#' @param se A Seurat object containing enrichment scores in the metadata.
#' @param perm.mat A permutation matrix of enrichment scores per cluster (global medians).
#' @param perm.mat.window50 A permutation matrix based on ranked gene windows (optional but preloaded).
#' @param window_rank_list A list or vector of rank indices used to define the enrichment windows.
#' @param ot_gene_set_label A string label for the gene set (e.g., `"Genetic"` or `"Drugs"`).
#' @param disease_abbr A short string indicating the trait or disease (e.g., `"SCZ"` or `"ALZ"`).
#' @param cluster_col The column in the Seurat metadata that contains cluster labels (default is `"seurat_clusters"`).
#' @param clusters_to_plot Optional vector of indices to subset and plot specific clusters only.
#'
#' @returns No return value. Produces a series of ggplot2 density plots, one per cluster.
#' @export
#'
#' @examples
#' plotScoreDist(
#'   se = se,
#'   perm.mat = perm.mat.genetic.data,
#'   perm.mat.window50 = perm.mat.window50.data,
#'   window_rank_list = window50_rank_list_SCZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "SCZ",
#'   cluster_col = "supercluster_term"
#' )
plotScoreDist <- function(se, 
                          perm.mat, 
                          perm.mat.window50, 
                          window_rank_list, 
                          ot_gene_set_label = "Genetic", 
                          disease_abbr = "ALZ",
                          cluster_col = "seurat_clusters",
                          clusters_to_plot = NULL) {
  
  theme_set(theme_classic())
  
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  cluster_labels <- unique(se@meta.data[[cluster_col]])
  
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
  start_with_abr_label <- grep(term, colnames(se@meta.data), value = TRUE)
  if (length(start_with_abr_label) == 0) {
    stop("No matching columns found for term")
  }
  
  for (cluster_label in cluster_labels) {
    message("\nProcessing cluster: ", cluster_label)
    
    # Get cells for this cluster
    cells_in_cluster <- rownames(se@meta.data)[se@meta.data[[cluster_col]] == cluster_label]
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
    median_score_NULL <- data.frame(Median_scores = perm.mat[[cluster_name]])
    ot_name <- paste0("OpenTargets_", disease_abbr, "_", ot_gene_set_label, "_1")
    
    if (!(ot_name %in% colnames(se@meta.data))) {
      warning("Missing OT score column: ", ot_name, " in cluster: ", cluster_label)
      next
    }
    
    actual_median <- median(se@meta.data[cells_in_cluster, ot_name], na.rm = TRUE)
    
    p <- ggplot(median_score_NULL, aes(x = Median_scores)) +
      geom_density(fill = "grey", alpha = 0.8) +
      geom_vline(xintercept = actual_median, color = "red", size = 1) +
      ggtitle(paste(ot_gene_set_label, disease_abbr, "Cluster:", cluster_label)) +
      xlab("Median Score") +
      ylab("Density")
    
    print(p)
  }
}
