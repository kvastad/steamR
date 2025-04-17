#' Generate Sliding Windows from a Ranked Gene List
#'
#' Creates overlapping gene windows from a ranked character vector. This is commonly
#' used to segment ranked gene lists (e.g., by association score) into fixed-size windows
#' for downstream enrichment analysis. Only full windows of the specified size are returned.
#'
#' @param char_vector A character vector of ranked gene symbols.
#' @param window_size Integer. The number of genes in each window.
#' @param step_size Integer. The step size between window starts (i.e., the number of genes to slide).
#' @param disease_abbr Character. A short abbreviation for the disease/trait (e.g., `"ALZ"`).
#' @param ot_gene_set_label Character. A label for the gene set type (e.g., `"Genetic"`).
#'
#' @return A named list of character vectors, where each vector represents a gene window.
#'         Names follow the format `"ALZ_Genetic_Rank1"`, `"ALZ_Genetic_Rank2"`, etc.
#'
#' @export
#'
#' @examples
#' ranked_genes <- c("Gene1", "Gene2", "Gene3", "Gene4", "Gene5", "Gene6")
#' generate_sliding_windows(ranked_genes, window_size = 3, step_size = 2, disease_abbr = "ALZ", ot_gene_set_label = "Genetic")
generate_sliding_windows <- function(char_vector, window_size, step_size, disease_abbr, ot_gene_set_label) {
  windows <- list()
  
  for (i in seq(1, length(char_vector) - window_size + 1, by = step_size)) {
    windows[[length(windows) + 1]] <- char_vector[i:(i + window_size - 1)]
  }
  
  Rank_names <- rep(paste0(disease_abbr, "_", ot_gene_set_label, "_", "Rank"), length(windows))
  Rank_names <- paste0(Rank_names, sep = 1:length(Rank_names))
  Rank_names <- paste0(Rank_names, sep = "_")
  names(windows) <- Rank_names
  
  return(windows)
}