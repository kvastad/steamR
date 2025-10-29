#' Plot Cluster Distributions for MetaCells or MetaSpots
#'
#' This function takes a list containing `cell_types` and `metacells` (Seurat objects)
#' and plots histograms of the number of units (cells or spots) per cluster
#' for each type. It also prints summary information to the console.
#'
#' @param metaspots A list containing:
#'   \describe{
#'     \item{cell_types}{Character vector of cell type names.}
#'     \item{metacells}{List of Seurat objects corresponding to each cell type.}
#'   }
#'
#' @return Invisible NULL. The function is called for its side effect of plotting
#'   and printing cluster summaries.
#'
#' @examples
#' \dontrun{
#' plot_metaclusters(sorted_metaspots)
#' }
#'
#' @export
plot_metaclusters <- function(metaspots) {
  cell_types <- metaspots$cell_types
  metacells <- metaspots$metacells
  
  invisible(lapply(seq_along(metacells), function(i) {
    se_obj <- metacells[[i]]   # Seurat object
    ct_name <- cell_types[i]
    
    # Extract cluster assignments from Seurat metadata
    clusters <- se_obj$seurat_clusters
    
    # Number of clusters
    num_clusters <- length(unique(clusters))
    cat("Number of clusters for", ct_name, ":", num_clusters, "\n")
    
    # Cluster distribution
    cluster_distribution <- table(clusters)
    cat("Cluster distribution:\n")
    print(cluster_distribution)
    
    # Plot histogram
    hist(cluster_distribution,
         main = paste("Cluster distribution for:", ct_name),
         xlab = "Units in MetaCluster",
         ylab = "Number of MetaClusters",
         col = "lightblue",
         border = "black")
    
    # Add mean and median markers
    abline(v = mean(cluster_distribution), col = "red", lwd = 2, lty = 2)
    abline(v = median(cluster_distribution), col = "blue", lwd = 2, lty = 2)
    
    cat(strrep("_", 80), "\n")
  }))
}
