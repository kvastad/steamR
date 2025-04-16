#' Rank Score Distribution Plots for Enrichment Structures
#'
#' Generates line plots of median enrichment scores across ranked windows 
#' for each cluster. This helps visualize how enrichment scores vary 
#' by window rank and how they compare to permutation null distributions.
#'
#' @param se A Seurat object containing cluster metadata and enrichment scores.
#' @param perm.mat A permutation matrix containing global median scores for each cluster.
#' @param perm.mat.window50 A permutation matrix for window-based (e.g., 50-gene) median scores.
#' @param window_rank_list A list or vector of window ranks used to generate enrichment scores.
#' @param ot_gene_set_label A string identifying the gene set type (e.g., "Genetic").
#' @param disease_abbr A short abbreviation for the disease or trait (e.g., "SCZ").
#' @param cluster_col The name of the metadata column specifying clusters 
#'        (default is `"seurat_clusters"`).
#'
#' @returns Line plots per cluster showing the structure of rank-wise median scores.
#' @export
#'
#' @examples
#' rankPlotForClusters(
#'   se = se,
#'   perm.mat = perm.mat.genetic.data,
#'   perm.mat.window50 = perm.mat.window50.data,
#'   window_rank_list = window50_rank_list_SCZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "SCZ",
#'   cluster_col = "supercluster_term"
#' )
rankPlotForClusters <- function(se, perm.mat, perm.mat.window50, window_rank_list, 
                                ot_gene_set_label, disease_abbr, cluster_col = "seurat_clusters") {
  
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  se_subset <- list()
  meta_data_subset <- list()
  meta_data_abr_label <- list()
  
  # Ensure cluster_col exists
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "does not exist in the Seurat object metadata."))
  }
  
  cluster_ids <- unique(se@meta.data[[cluster_col]])
  
  # Temporarily set identity class to the selected clustering column
  se <- SetIdent(se, value = cluster_col)
  
  for (i in cluster_ids) {
    subset_name <- paste0("se_", i)
    se_subset[[subset_name]] <- subset(se, idents = i)
    meta_data_subset[[subset_name]] <- se_subset[[subset_name]][[]]
    start_with_abr_label <- grep(term, colnames(meta_data_subset[[subset_name]]), value = TRUE)
    meta_data_abr_label[[subset_name]] <- subset(meta_data_subset[[subset_name]],
                                                 select = start_with_abr_label)
  }
  
  for (i in cluster_ids) {
    subset_name <- paste0("se_", i)
    abr_label_structure_i <- NULL
    
    for (k in seq_len(ncol(meta_data_abr_label[[subset_name]]))) {
      meta_data_vector <- data.frame(
        Rank = as.character(k),
        Median = median(meta_data_abr_label[[subset_name]][, k], na.rm = TRUE)
      )
      abr_label_structure_i <- rbind(abr_label_structure_i, meta_data_vector)
    }
    
    rownames(abr_label_structure_i) <- seq_len(nrow(abr_label_structure_i))
    abr_label_structure_i <- as.data.frame(abr_label_structure_i)
    
    # Plot
    plot(
      abr_label_structure_i$Rank,
      abr_label_structure_i$Median,
      main = paste(disease_abbr, ot_gene_set_label, "Structure", i),
      xlab = "Rank",
      ylab = "Score",
      type = "o",
      ylim = c(-0.3, 0.3)
    )
    axis(1, at = seq(1, 200, by = 1), cex.axis = 0.5)
    
    abline(h = max(perm.mat.window50[, 1]), col = "darkred", lwd = 2, lty = 2)
    abline(h = quantile(as.numeric(perm.mat.window50[, 1]), probs = 0.95), col = "red", lwd = 2, lty = 2)
    abline(h = median(perm.mat.window50[, 1]), col = "blue", lwd = 2, lty = 2)
    abline(h = quantile(as.numeric(perm.mat.window50[, 1]), probs = 0.05), col = "red", lwd = 2, lty = 2)
    abline(h = min(perm.mat.window50[, 1]), col = "darkred", lwd = 2, lty = 2)
  }
}