#' Order cells by a predefined biological lineage
#'
#' Assigns consecutive ranks to cells according to a manually specified order
#' of cell types. Within each cell type, cells retain their relative ordering
#' based on a numeric ranking variable, such as diffusion pseudotime (DPT).
#' This allows known biological lineage information to be incorporated into
#' downstream analyses of stratified trajectory enrichment patterns.
#'
#' The Seurat object should first be subset to the lineage being analysed.
#' Every cell type present in the object's metadata must be included in
#' `lineage_order`; otherwise, the function stops with an informative error.
#'
#' @param se A Seurat object containing cells from the lineage being
#'   analysed.
#' @param cell_rank_col Character scalar naming a metadata column containing
#'   numeric cell ranks or values used to order cells within each cell type.
#'   Defaults to `"cell_rank"`.
#' @param cell_type_col Character scalar naming a metadata column containing
#'   cell-type annotations. Defaults to `"cell_anno"`.
#' @param lineage_order Character vector specifying the desired order of cell
#'   types, from earliest to latest in the lineage.
#'
#' @return A numeric vector of consecutive ranks, in the same cell order as
#'   `se@meta.data`. Names correspond to the metadata row names.
#'
#' @examples
#' \dontrun{
#' lineage_order <- c("NSC", "INP", "NB", "YGN", "GN")
#'
#' # Subset the Seurat object to the lineage being analysed.
#' lineage_cells <- subset(
#'   se,
#'   subset = cell_anno %in% lineage_order
#' )
#'
#' lineage_cells$cell_rank_ordered <- order_cells_by_lineage(
#'   lineage_cells,
#'   cell_rank_col = "DPT",
#'   cell_type_col = "cell_anno",
#'   lineage_order = lineage_order
#' )
#' }
#'
#' @export
order_cells_by_lineage <- function(se,
                                   cell_rank_col = "cell_rank",
                                   cell_type_col = "cell_anno",
                                   lineage_order) {
  if (!is.character(cell_rank_col) || length(cell_rank_col) != 1L ||
      is.na(cell_rank_col) || !nzchar(cell_rank_col)) {
    stop("`cell_rank_col` must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.character(cell_type_col) || length(cell_type_col) != 1L ||
      is.na(cell_type_col) || !nzchar(cell_type_col)) {
    stop("`cell_type_col` must be a single non-empty character string.",
         call. = FALSE)
  }
  
  if (!is.character(lineage_order) || length(lineage_order) == 0L ||
      anyNA(lineage_order) || any(!nzchar(lineage_order))) {
    stop("`lineage_order` must be a non-empty character vector without missing or empty values.",
         call. = FALSE)
  }
  
  if (anyDuplicated(lineage_order)) {
    stop("`lineage_order` must not contain duplicate cell types.",
         call. = FALSE)
  }
  
  if (!"meta.data" %in% methods::slotNames(se)) {
    stop("`se` must have a `meta.data` slot, as in a Seurat object.",
         call. = FALSE)
  }
  
  meta <- se@meta.data
  
  if (!(cell_rank_col %in% colnames(meta))) {
    stop(sprintf("Metadata column `%s` was not found.", cell_rank_col),
         call. = FALSE)
  }
  
  if (!(cell_type_col %in% colnames(meta))) {
    stop(sprintf("Metadata column `%s` was not found.", cell_type_col),
         call. = FALSE)
  }
  
  cell_ranks <- meta[[cell_rank_col]]
  if (!is.numeric(cell_ranks) || !is.null(dim(cell_ranks))) {
    stop(sprintf("Metadata column `%s` must be numeric.", cell_rank_col),
         call. = FALSE)
  }
  
  cell_types <- as.character(meta[[cell_type_col]])
  
  if (anyNA(cell_ranks) || any(!is.finite(cell_ranks))) {
    stop(sprintf("Metadata column `%s` must not contain missing or non-finite values.",
                 cell_rank_col),
         call. = FALSE)
  }
  
  if (anyNA(cell_types)) {
    stop(sprintf("Metadata column `%s` must not contain missing cell-type annotations.",
                 cell_type_col),
         call. = FALSE)
  }
  
  unlisted_types <- setdiff(unique(cell_types), lineage_order)
  if (length(unlisted_types) > 0L) {
    stop(
      sprintf(
        paste0(
          "The Seurat object contains cell types not listed in `lineage_order`: %s. ",
          "Subset the object to the lineage of interest or include all present ",
          "cell types in `lineage_order`."
        ),
        paste(unlisted_types, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  
  new_ranks <- numeric(length(cell_ranks))
  names(new_ranks) <- rownames(meta)
  counter <- 1L
  
  for (ct in lineage_order) {
    idx <- which(cell_types == ct)
    if (length(idx) == 0L) {
      next
    }
    
    # Stable ordering within each cell type; input order breaks ties.
    sorted_idx <- idx[order(cell_ranks[idx], seq_along(idx))]
    n <- length(sorted_idx)
    
    new_ranks[sorted_idx] <- seq.int(from = counter, length.out = n)
    counter <- counter + n
  }
  
  return(new_ranks)
}
