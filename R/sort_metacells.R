#' Sort MetaCells or MetaSpots by Type
#'
#' This function takes a Seurat object and a list of metacells, extracts the unique
#' cluster annotations (cell types or spots), sorts them alphabetically, and
#' reorders the metacells accordingly.
#'
#' @param se A Seurat object containing the `cluster_anno` column in its metadata.
#' @param metacells A list of Seurat objects corresponding to the clusters in `se`.
#' @param cluster_anno Character string specifying the metadata column to use for sorting.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{cell_types}{Character vector of sorted cell type names.}
#'     \item{metacells}{List of metacells reordered to match sorted cell types.}
#'   }
#'
#' @examples
#' \dontrun{
#' sorted <- sort_metacells(se, metacells)
#' sorted$cell_types
#' sorted$metacells
#' }
#'
#' @export
sort_metacells <- function(se, metacells, cluster_anno = "cluster_anno") {
  # Check that the metadata column exists
  if (!cluster_anno %in% colnames(se@meta.data)) {
    stop(paste("Column", cluster_anno, "not found in Seurat metadata"))
  }
  
  # Remaining cell types in the Seurat object
  remaining_cell_types <- unique(se@meta.data[[cluster_anno]])
  
  # If metacells is unnamed, automatically name them from their Seurat metadata
  if (is.null(names(metacells))) {
    names(metacells) <- sapply(metacells, function(x) {
      unique_vals <- unique(x@meta.data[[cluster_anno]])
      if (length(unique_vals) != 1) {
        stop("Each Seurat object in metacells must have exactly one unique cell type in the specified column")
      }
      unique_vals
    })
  }
  
  # Filter metacells to only those remaining
  metacells <- metacells[names(metacells) %in% remaining_cell_types]
  
  # Sort remaining cell types alphabetically
  sorted_indices <- order(names(metacells))
  metacells <- metacells[sorted_indices]
  cell_types <- names(metacells)
  
  list(cell_types = cell_types, metacells = metacells)
}