#' Sort MetaCells or MetaSpots by Type
#'
#' This function takes a Seurat object and a list of metacells, extracts the unique
#' cluster annotations (cell types or spots), sorts them alphabetically, and
#' reorders the metacells accordingly.
#'
#' @param se A Seurat object containing the `cluster_anno` column in its metadata.
#' @param metacells A list of Seurat objects corresponding to the clusters in `se`.
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
sort_metacells <- function(se, metacells) {
  # Extract unique cluster annotations (cell types or spots) and sort them alphabetically
  cell_types <- unique(se$cluster_anno)
  sorted_indices <- order(cell_types)
  cell_types <- cell_types[sorted_indices]
  metacells <- metacells[sorted_indices]
  
  # Return both updated objects
  list(cell_types = cell_types, metacells = metacells)
}
