#' Rank Score Distribution Plots for Enrichment Structures
#'
#' Generates line plots of median enrichment scores across ranked windows 
#' for each cluster. This helps visualize how enrichment scores vary 
#' by window rank and how they compare to permutation null distributions.
#'
#' @param window_results Output from WindowRankEnrichmentAnalysis containing pre-calculated median scores and quantiles
#' @param perm.mat.window A null distribution permutation matrix for window-based median scores, used to calculate and visualize the null distribution's min, median and max values. This permutation matrix must be generated from random gene lists of the same size as the window size being tested.
#' @param clusters_to_plot Optional numeric vector indicating specific clusters to plot. If NULL, all are plotted.
#' @param ylim A numeric vector of length two specifying the y-axis limits (default: c(-0.3, 0.3))
#'
#' @returns Line plots per cluster showing the structure of rank-wise median scores.
#' @export
#'
#' @examples
#' window_results <- WindowRankEnrichmentAnalysis(
#'   se = se,
#'   perm.mat.window.data = perm.mat.window.data,
#'   window_rank_list = window_rank_list_ALZ_Drugs,
#'   ot_gene_set_label = "Drugs",
#'   disease_abbr = "ALZ"
#' )
#' 
#' rankPlotForClusters(
#'   window_results = window_results,
#'   perm.mat.window = perm.mat.window.data,
#'   ylim = c(-0.3, 0.3)
#' )
rankPlotForClusters <- function(window_results, 
                               perm.mat.window, 
                               clusters_to_plot = NULL,
                               ylim = c(-0.3, 0.3)) {
    
    # Get unique clusters in their original order (maintained by factor levels)
    unique_clusters <- levels(window_results$cluster)
    if (is.null(unique_clusters)) {
        unique_clusters <- sort(unique(window_results$cluster))
    }
    
    # Filter clusters if specified
  if (!is.null(clusters_to_plot)) {
        cluster_names <- paste0("cluster_", clusters_to_plot)
        unique_clusters <- unique_clusters[unique_clusters %in% cluster_names]
    }
    
    # Create plots in the correct order
    for (cluster in unique_clusters) {
        cluster_data <- window_results[window_results$cluster == cluster, ]
        
        # Sort by window number to ensure correct plotting order
        cluster_data <- cluster_data[order(cluster_data$window), ]
        
        plot(cluster_data$window,
             cluster_data$observed_score,
             main = paste("Cluster:", gsub("cluster_", "", cluster)),
      xlab = "Window Rank",
      ylab = "Median Score",
      type = "o",
      ylim = ylim,
             xaxt = "n")
        
        axis(1, at = seq_along(cluster_data$window), cex.axis = 0.5)
    
        abline(h = cluster_data$q95[1], col = "red", lwd = 2, lty = 2)
        abline(h = cluster_data$q05[1], col = "red", lwd = 2, lty = 2)
        
        if (!is.null(perm.mat.window)) {
            null_dist <- perm.mat.window[[gsub("cluster_", "", cluster)]]
            if (!is.null(null_dist)) {
                abline(h = max(null_dist), col = "darkred", lwd = 2, lty = 2)
                abline(h = min(null_dist), col = "darkred", lwd = 2, lty = 2)
                abline(h = median(null_dist), col = "blue", lwd = 2, lty = 2)
            }
        }
  }
}
