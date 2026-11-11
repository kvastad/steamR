#' Generate trajectory MetaCell windows
#'
#' Divides the cells in a Seurat object into consecutive windows after
#' ordering them by a specified metadata column. Each window is returned
#' as a separate Seurat object.
#'
#' @param se Seurat object containing the cells to be ordered.
#' @param rank.col Name of the metadata column used to order cells along
#'   the trajectory.
#' @param window.size Number of cells included in each window.
#' @param step.size Number of cells by which the window start advances.
#'   If `NULL`, defaults to `window.size`.
#'
#' @return A named list of Seurat objects, with one element for each
#'   trajectory window. Element names indicate the rank range included
#'   in each window.
#'
#' @export
generate_trajectory_metacells <- function(
    se,
    rank.col = "cell_rank",
    window.size = 10,
    step.size = NULL
) {
  
  # Default step size = window size
  if (is.null(step.size)) step.size <- window.size
  
  # Extract metadata
  meta <- se@meta.data
  
  if (!rank.col %in% colnames(meta)) {
    stop(paste("Metadata column", rank.col, "not found in object"))
  }
  
  # Order cells only by rank.col
  meta <- meta[order(meta[[rank.col]]), , drop = FALSE]
  
  ranked.cells <- rownames(meta)
  
  # Generate window start positions
  starts <- seq(1, length(ranked.cells), by = step.size)
  
  # Create list of Seurat subsets
  window.list <- lapply(starts, function(s) {
    e <- s + window.size - 1
    cells.in.window <- ranked.cells[s:min(e, length(ranked.cells))]
    subset(se, cells = cells.in.window)
  })
  
  # Name each element
  names(window.list) <- paste0(
    "rank_", starts,
    "_to_",
    pmin(starts + window.size - 1, length(ranked.cells))
  )
  
  return(window.list)
}
