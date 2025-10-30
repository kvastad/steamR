#' Plot Cluster Distributions for MetaCells or MetaSpots
#'
#' This function takes a list containing `cell_types` and `metacells` (Seurat objects)
#' and plots histograms of the number of units (cells or spots) per cluster
#' for each type and optionally adds mean and/or median reference lines. 
#' It also prints summary information to the console.
#'
#'
#' @param metaspots A list with elements:
#'   \describe{
#'     \item{cell_types}{Character or numeric vector of cell type names.}
#'     \item{metacells}{List of Seurat objects corresponding to each cell type.}
#'   }
#' @param show_mean Logical; whether to display a red dashed line for the mean cluster size. Default is \code{TRUE}.
#' @param show_median Logical; whether to display a blue dashed line for the median cluster size. Default is \code{TRUE}.
#'
#' @return Invisibly returns \code{NULL}. Generates plots as a side effect.
#'
#' @examples
#' \dontrun{
#' plot_metaclusters(sorted_metacells, show_mean = TRUE, show_median = FALSE)
#' }
#'
#' @export
plot_metaclusters <- function(metaspots, show_mean = TRUE, show_median = TRUE) {
  cell_types <- metaspots$cell_types
  metacells <- metaspots$metacells
  
  invisible(lapply(seq_along(metacells), function(i) {
    se_obj <- metacells[[i]]
    ct_name <- cell_types[i]
    
    # Extract cluster assignments
    clusters <- se_obj$seurat_clusters
    
    # Number of clusters
    num_clusters <- length(unique(clusters))
    cat("Number of clusters for", ct_name, ":", num_clusters, "\n")
    
    # Cluster distribution
    cluster_distribution <- table(clusters)
    cat("Cluster distribution:\n")
    print(cluster_distribution)
    
    # Detect whether Seurat object is spatial or single-cell
    x_label <- if ("images" %in% slotNames(se_obj) && length(se_obj@images) > 0) {
      "Spots in MetaCluster"
    } else {
      "Cells in MetaCluster"
    }
    
    # Plot histogram
    hist(cluster_distribution,
         main = paste("Cluster distribution for:", ct_name),
         xlab = x_label,
         ylab = "Number of MetaClusters",
         col = "lightblue",
         border = "black")
    
    # Optionally add mean and median lines
    if (show_mean) {
      abline(v = mean(cluster_distribution), col = "red", lwd = 2, lty = 2)
    }
    if (show_median) {
      abline(v = median(cluster_distribution), col = "blue", lwd = 2, lty = 2)
    }
    
    cat(strrep("_", 80), "\n")
  }))
}
