#' Visualize Rank-Based Enrichment Scores with Spatial Feature Plots
#'
#' Generates spatial `SpatialFeaturePlot()` visualizations for enrichment scores computed 
#' over ranked gene windows. Allows plotting either all clusters together or one 
#' at a time.
#'
#' @param se A Seurat object with spatial transcriptomics data and precomputed enrichment scores.
#' @param window_rank_list A list of gene rank windows used for enrichment calculation.
#' @param ot_gene_set_label Character. The gene set label (e.g., "Genetic").
#' @param disease_abbr Character. Abbreviation for the disease/trait (e.g., "SCZ").
#' @param cluster_number Optional. A specific cluster to visualize (used when `all_clusters = FALSE`).
#' @param ranks_per_plot Number of rank plots to show per panel (default is `6`).
#' @param spot_alpha Alpha value for spot transparency in plots (default is `0.5`).
#' @param point_size Point size for `SpatialFeaturePlot()` (default is `1`).
#' @param all_clusters Logical. Whether to visualize all clusters together (default is `FALSE`).
#' @param cluster_anno Name of the metadata column to use for coloring (default: "cluster_anno").
#'
#' @returns No return value. Generates spatial enrichment plots using `SpatialFeaturePlot()` per cluster/rank.
#' @export
#'
#' @examples
#' featurePlotForRanks(
#'   se = se,
#'   window_rank_list = window50_rank_list_SCZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "SCZ",
#'   cluster_number = 5,
#'   ranks_per_plot = 4,
#'   all_clusters = FALSE
#' )
featurePlotForRanks <- function(
    se, 
    window_rank_list, 
    ot_gene_set_label, 
    disease_abbr, 
    cluster_number = NULL, 
    ranks_per_plot = 6,
    spot_alpha = 0.5, 
    point_size = 1,
    all_clusters = FALSE,
    cluster_anno = "cluster_anno"
) {
  
  if (all_clusters) {
    cluster_numbers <- unique(se@meta.data[[cluster_anno]])
    
    # Calculate number of plots needed
    num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
    
    # Process each set of ranks
    for (plot_idx in 1:num_plots) {
      start_rank <- (plot_idx - 1) * ranks_per_plot + 1
      end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
      
      # Create feature names for this set of ranks
      features <- sapply(start_rank:end_rank, function(rank) {
        paste(disease_abbr, ot_gene_set_label, paste("Rank", rank, sep = ""), "1", sep = "_")
      })
      
      print(SpatialFeaturePlot(se, 
                           features = features,
                           pt.size.factor = point_size,
                           ncol = min(4, length(features)),
                           image.alpha = 0) +
            scale_color_gradient(low = "blue", high = "red"))
    }
    
  } else {
    # For single cluster visualization
    cluster_numbers <- cluster_number
    
    # Calculate number of plots needed
      num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
      
    # Process each set of ranks
      for (plot_idx in 1:num_plots) {
        start_rank <- (plot_idx - 1) * ranks_per_plot + 1
        end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
        
      # Create feature names for this set of ranks
      features <- sapply(start_rank:end_rank, function(rank) {
        paste(disease_abbr, ot_gene_set_label, paste("Rank", rank, sep = ""), "1", sep = "_")
      })
      
      # Create a mask for the selected cluster
      cells_use <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] %in% cluster_numbers]
      
      print(SpatialFeaturePlot(se, 
                             features = features,
                           cells = cells_use,
                           pt.size.factor = point_size,
                           ncol = min(4, length(features)),
                           image.alpha = 0) +
            scale_color_gradient(low = "blue", high = "red"))
    }
  }
}
