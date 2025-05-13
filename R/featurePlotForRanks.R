#' Visualize Rank-Based Enrichment Scores with Spatial Feature Plots
#'
#' Generates spatial visualizations for enrichment scores computed over ranked gene windows.
#' Works with both standard Seurat objects and STUtility-generated Seurat objects.
#'
#' @param se A Seurat object with spatial transcriptomics data and precomputed enrichment scores.
#' @param window_rank_list A list of gene rank windows used for enrichment calculation.
#' @param ot_gene_set_label Character. The gene set label (e.g., "Genetic").
#' @param disease_abbr Character. Abbreviation for the disease/trait (e.g., "SCZ").
#' @param cluster_number Optional. A specific cluster to visualize (used when `all_clusters = FALSE`).
#' @param ranks_per_plot Number of rank plots to show per panel (default is `6`).
#' @param spot_alpha Alpha value for spot transparency in plots (default is `0.5`).
#' @param point_size Point size for spatial plots (default is `1`).
#' @param all_clusters Logical. Whether to visualize all clusters together (default is `FALSE`).
#' @param cluster_anno Name of the metadata column to use for coloring (default: "cluster_anno").
#' @param seurat_type Optional. Type of Seurat object ("seurat" or "stutility"). If NULL, will attempt to auto-detect.
#'
#' @returns No return value. Generates spatial enrichment plots per cluster/rank.
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
    cluster_anno = "cluster_anno",
    seurat_type = NULL
) {
  # Auto-detect Seurat object type if not specified
  if (is.null(seurat_type)) {
    if ("Staffli" %in% names(se@tools) || 
        any(grepl("^STUtility", class(se))) ||
        "STUtility" %in% names(se@misc)) {
      seurat_type <- "stutility"
    } else {
      seurat_type <- "seurat"
    }
  } else {
    if (!seurat_type %in% c("seurat", "stutility")) {
      stop("seurat_type must be either 'seurat' or 'stutility'")
    }
  }
  
  # Show warning only for STUtility objects when STUtility package is not available
  if (seurat_type == "stutility" && !requireNamespace("STUtility", quietly = TRUE)) {
    warning("Package 'STUtility' is required for spatial visualization of STUtility generated Seurat objects. Some features may not work as expected.")
  }
  
  # Ensure cluster_anno exists
  if (!(cluster_anno %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_anno, "does not exist in the Seurat object metadata."))
  }
  
  if (all_clusters) {
    # Calculate number of plots needed
    num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
    
    for (plot_idx in 1:num_plots) {
      start_rank <- (plot_idx - 1) * ranks_per_plot + 1
      end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
      
      features <- paste(disease_abbr, ot_gene_set_label, paste("Rank", start_rank:end_rank, sep = ""), "1", sep = "_")
      
      if (seurat_type == "stutility") {
        print(ST.FeaturePlot(se, 
                           features = features,
                           pt.alpha = spot_alpha,
                           pt.size = point_size,
                           grid.ncol = 4,
                           cols = c("black", "darkblue", "cyan", "yellow", "red", "darkred"),
                           value.scale = "all"))
      } else {
        print(SpatialFeaturePlot(se, 
                               features = features,
                               pt.size.factor = point_size,
                               ncol = min(4, length(features)),
                               image.alpha = 0) +
              scale_color_gradient(low = "blue", high = "red"))
      }
    }
    
  } else {
    # For single cluster visualization
    cluster_numbers <- cluster_number
    
    # Calculate number of plots needed
    num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
    
    for (plot_idx in 1:num_plots) {
      start_rank <- (plot_idx - 1) * ranks_per_plot + 1
      end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
      
      features <- paste(disease_abbr, ot_gene_set_label, paste("Rank", start_rank:end_rank, sep = ""), "1", sep = "_")
      
      # Create a mask for the selected cluster
      cells_use <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] %in% cluster_numbers]
      
      if (seurat_type == "stutility") {
        print(ST.FeaturePlot(se, 
                           features = features,
                           cells = cells_use,
                           pt.alpha = spot_alpha,
                           pt.size = point_size,
                           grid.ncol = 4,
                           cols = c("black", "darkblue", "cyan", "yellow", "red", "darkred"),
                           value.scale = "all"))
      } else {
        se_subset <- suppressWarnings(subset(se, cells = cells_use))
        suppressWarnings(print(
          SpatialFeaturePlot(se_subset, 
            features = features,
            pt.size.factor = point_size,
            ncol = min(4, length(features)),
            image.alpha = 0
          ) + scale_color_gradient(low = "blue", high = "red")
        ))
      }
    }
  }
}
