#' Plot Enrichment Score Distribution Against Null
#'
#' Generates a density plot of the NULL gene set median scores for each cluster,
#' and overlays the actual enrichment score to assess statistical significance.
#'
#' @param se A Seurat object containing enrichment scores in metadata.
#' @param perm.mat A permutation matrix of actual enrichment scores per cluster.
#' @param perm.mat.window50 A window-level permutation matrix (e.g. 50-gene ranks).
#' @param window_rank_list A list of ranked gene windows used for enrichment.
#' @param ot_gene_set_label A character label for the gene set type (e.g., `"Genetic"`).
#' @param disease_abbr A character abbreviation for the trait/disease (e.g., `"SCZ"`).
#' @param cluster_col Metadata column to use for clustering (default is `"seurat_clusters"`).
#'
#' @returns No return value. Generates enrichment score plots for each cluster.
#' @export
#'
#' @examples
#' plotScoreDist(
#'   se = se,
#'   perm.mat = perm.mat.genetic.data,
#'   perm.mat.window50 = perm.mat.window50.data,
#'   window_rank_list = window50_rank_list_SCZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "SCZ",
#'   cluster_col = "seurat_clusters"
#' )
plotScoreDist <- function(se, 
                          perm.mat, 
                          perm.mat.window50, 
                          window_rank_list, 
                          ot_gene_set_label = "Genetic", 
                          disease_abbr = "ALZ",
                          cluster_col = "seurat_clusters") {
  
  theme_set(ggplot2::theme_classic())
  
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  cluster_labels <- unique(se@meta.data[[cluster_col]])
  
  for (cluster_label in cluster_labels) {
    message("\nProcessing cluster: ", cluster_label)
    
    cells_in_cluster <- colnames(se)[se@meta.data[[cluster_col]] == cluster_label]
    if (length(cells_in_cluster) == 0) {
      warning("No cells found for cluster: ", cluster_label)
      next
    }
    
    se_subset <- subset(se, cells = cells_in_cluster)
    meta_data_subset <- se_subset[[]]
    
    start_with_abr_label <- grep(term, colnames(meta_data_subset), value = TRUE)
    if (length(start_with_abr_label) == 0) {
      warning("No matching columns found for term in cluster: ", cluster_label)
      next
    }
    
    meta_data_abr_label <- subset(meta_data_subset, select = start_with_abr_label)
    abr_label_structure_i <- data.frame(
      Rank = as.character(1:ncol(meta_data_abr_label)),
      Median = apply(meta_data_abr_label, 2, median)
    )
    
    cluster_name <- as.character(cluster_label)
    
    if (!(cluster_name %in% colnames(perm.mat))) {
      warning("Cluster name not found in permutation matrix: ", cluster_name)
      next
    }
    
    Rank_names <- paste0("Rank_nr", 1:length(window_rank_list))
    p.mat.cluster_i <- matrix(NA, nrow = length(window_rank_list), ncol = 2,
                              dimnames = list(Rank_names, c("p.val", "p.val.adj")))
    
    permutation.nr <- nrow(perm.mat)
    nr_of_tests <- length(cluster_labels)
    
    for (j in seq_along(window_rank_list)) {
      OT_genetic_window <- abr_label_structure_i[j, "Median"]
      perm_bigger <- perm.mat.window50[perm.mat.window50[[cluster_name]] >= OT_genetic_window, ]
      
      p.raw <- (nrow(perm_bigger) + 1) / (permutation.nr + 1)
      p.adj <- p.adjust(p.raw, method = "bonferroni", n = nr_of_tests)
      
      p.mat.cluster_i[j, 1] <- p.raw
      p.mat.cluster_i[j, 2] <- p.adj
    }
    
    median_score_NULL <- data.frame(Median_scores = perm.mat[[cluster_name]])
    ot_name <- paste0("OpenTargets_", disease_abbr, "_", ot_gene_set_label, "_1")
    
    if (!(ot_name %in% colnames(se_subset[[]]))) {
      warning("Missing OT score column: ", ot_name, " in cluster: ", cluster_label)
      next
    }
    
    actual_median <- median(as.numeric(unlist(se_subset[[ot_name]])))
    perm_bigger <- perm.mat[perm.mat[[cluster_name]] >= actual_median, ]
    p.global <- (nrow(perm_bigger) + 1) / (permutation.nr + 1)
    p.global.adj <- p.adjust(p.global, method = "bonferroni", n = nr_of_tests)
    
    annotation_text <- sprintf("p.val.adj = %.4f", p.global.adj)
    
    p <- ggplot2::ggplot(median_score_NULL, ggplot2::aes(x = Median_scores)) +
      ggplot2::geom_density(fill = "#69b3a2", alpha = 0.8) +
      ggplot2::geom_vline(xintercept = actual_median, color = "purple", size = 1) +
      ggplot2::annotate("text", x = 0, y = 0.01, label = annotation_text, hjust = 0) +
      ggplot2::ggtitle(paste(ot_gene_set_label, disease_abbr, "Structure", cluster_label)) +
      ggplot2::xlab("Median score NULL gene set") +
      ggplot2::ylab("Density")
    
    print(p)
  }
}
