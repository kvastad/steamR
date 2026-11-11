#' Rank cells by a numeric variable
#'
#' Assigns each cell a rank based on a numeric variable, such as diffusion
#' pseudotime (DPT), stored in a vector. By default, cells are ranked in
#' ascending order, so the cell with the lowest value receives rank 1.
#' Set `reverse = TRUE` to rank in descending order, so the cell with the
#' highest value receives rank 1.
#' 
#' Ties are resolved by input order: when cells have the same value,
#' the cell appearing first in `values` receives the lower rank.
#'
#' @param values Numeric vector of diffusion pseudotime values, with one
#'   value per cell. Missing values are not allowed.
#' @param reverse Logical scalar. If `FALSE` (default), rank values in
#'   ascending order. If `TRUE`, rank values in descending order.
#'
#' @return An integer vector of ranks with the same length and order as
#'   `values`.
#'
#' @examples
#' dpt <- c(cell_a = 0.8, cell_b = 0.2, cell_c = 0.5)
#' rank_cells(dpt)
#' # cell_a cell_b cell_c
#' #      3      1      2
#'
#' rank_cells(dpt, reverse = TRUE)
#' # cell_a cell_b cell_c
#' #      1      3      2
#'
#' @export
rank_cells <- function(values, reverse = FALSE) {
  if (!is.numeric(values) || !is.null(dim(values))) {
    stop("`values` must be a numeric vector.", call. = FALSE)
  }
  
  if (anyNA(values)) {
    stop("`values` must not contain missing values.", call. = FALSE)
  }
  
  if (!is.logical(reverse) || length(reverse) != 1L || is.na(reverse)) {
    stop("`reverse` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  
  if (reverse) {
    return(rank(-values, ties.method = "first"))
  }
  
  rank(values, ties.method = "first")
}
