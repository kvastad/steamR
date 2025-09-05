#' Calculate Gene Percentage per Cell
#'
#' Computes the percentage of total expression attributable to genes matching a given pattern,
#' across all cells in a Seurat object.
#'
#' @param se A Seurat object.
#' @param pattern A regular expression pattern used to match gene names (e.g., `"^MT-"` for mitochondrial genes).
#'
#' @returns A numeric vector of gene percentages (0–100) for each cell in the Seurat object.
#' @export
#'
#' @examples
#' # Calculate mitochondrial gene percentage
#' mito_pct <- calculate_gene_percentage(se, pattern = "^MT-")
calculate_gene_percentage <- function(se, pattern) {
  genes <- grep(pattern, rownames(se), value = TRUE)
  (Matrix::colSums(GetAssayData(se, slot = "counts")[genes, ]) /
      Matrix::colSums(GetAssayData(se, slot = "counts"))) * 100
}