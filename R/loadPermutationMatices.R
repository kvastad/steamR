#' Load Permutation Matrices from RDS Files
#'
#' Loads a list of permutation matrix `.rds` files and returns them as a named list of data frames or matrices.
#'
#' @param file_paths A character vector of file paths pointing to `.rds` files containing permutation matrices.
#'
#' @returns A named list of permutation matrices, where names correspond to the file base names.
#' @export
#'
#' @examples
#' file_paths <- list.files("data/permutations", pattern = "\\.rds$", full.names = TRUE)
#' permutation_matrices <- loadPermutationMatices(file_paths)
loadPermutationMatices <- function(file_paths) {
  matrices <- lapply(file_paths, readRDS)
  names(matrices) <- basename(file_paths)
  return(matrices)
}