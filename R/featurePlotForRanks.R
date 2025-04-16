#' Visualize Rank-Based Enrichment Scores with Spatial Feature Plots
#'
#' Generates spatial `ST.FeaturePlot()` visualizations for enrichment scores computed 
#' over ranked gene windows. Allows plotting either all clusters together or one 
#' at a time.
#'
#' @param se A Seurat object with spatial transcriptomics data and precomputed enrichment scores.
#' @param perm.mat A permutation matrix with baseline null distribution of median scores.
#' @param perm.mat.window50 A window-level permutation matrix for enrichment score thresholds.
#' @param window_rank_list A list of gene rank windows used for enrichment calculation.
#' @param ot_gene_set_label Character. The gene set label (e.g., `"Genetic"`).
#' @param disease_abbr Character. Abbreviation for the disease or trait (e.g., `"SCZ"`).
#' @param cluster_number Optional. A specific cluster to visualize (used when `all_clusters = FALSE`).
#' @param ranks_per_plot Number of rank plots to show per panel (default is `6`).
#' @param spot_alpha Alpha value for spot transparency in plots (default is `0.5`).
#' @param point_size Point size for `ST.FeaturePlot()` (default is `1`).
#' @param all_clusters Logical. Whether to visualize all clusters together (default is `FALSE`).
#'
#' @returns No return value. Generates spatial enrichment plots using `ST.FeaturePlot()` per cluster/rank.
#' @export
#'
#' @examples
#' featurePlotForRanks(
#'   se = se,
#'   perm.mat = perm.mat.genetic.data,
#'   perm.mat.window50 = perm.mat.window50.data,
#'   window_rank_list = window50_rank_list_SCZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "SCZ",
#'   cluster_number = 5,
#'   ranks_per_plot = 4,
#'   all_clusters = FALSE
#' )
featurePlotForRanks <- function(se, 
                                perm.mat, 
                                perm.mat.window50, 
                                window_rank_list, 
                                ot_gene_set_label, 
                                disease_abbr, 
                                cluster_number = NULL, 
                                ranks_per_plot = 6,
                                spot_alpha = 0.5, 
                                point_size = 1,
                                all_clusters = FALSE) {
  
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  
  se_subset <- list()
  meta_data_subset <- list()
  meta_data_abr_label <- list()
  
  if (all_clusters) {
    cluster_numbers <- unique(se$seurat_clusters)
    combined_meta_data <- do.call(rbind, lapply(cluster_numbers, function(i) {
      subset_name <- paste0("se_", i)
      se_subset[[subset_name]] <- SubsetSTData(se, expression = seurat_clusters == i)
      meta_data_subset[[subset_name]] <- se_subset[[subset_name]][[]]
      start_with_abr_label <- grep(term, colnames(meta_data_subset[[subset_name]]), value = TRUE)
      meta_data_abr_label[[subset_name]] <- subset(meta_data_subset[[subset_name]], 
                                                   select = start_with_abr_label)
      meta_data_abr_label[[subset_name]]
    }))
    
    abr_label_structure_i <- NULL
    for (k in 1:ncol(combined_meta_data)) {
      meta_data_vector <- data.frame(
        Rank = as.character(k), 
        Median = median(combined_meta_data[, k])
      )
      abr_label_structure_i <- rbind(abr_label_structure_i, meta_data_vector)
    }
    
    rownames(abr_label_structure_i) <- 1:ncol(combined_meta_data)
    abr_label_structure_i <- as.data.frame(abr_label_structure_i)
    
    Rank_names <- paste0("Rank_nr", 1:length(window_rank_list))
    p.mat.cluster_i <- matrix(NA, nrow = length(window_rank_list), ncol = 2, 
                              dimnames = list(Rank_names, c("p.val", "p.val.adj")))
    
    permutation.nr <- dim(perm.mat)[1]
    nr_of_tests <- 1
    
    for (j in 1:length(window_rank_list)) {
      OT_genetic_window <- as.numeric(abr_label_structure_i[j, 2])
      
      perm.mat.window50.cluster_i_bigger_window <- subset(
        perm.mat.window50, apply(perm.mat.window50, 1, function(row) any(row >= OT_genetic_window))
      )
      
      cluster_i_w_p.val <- (nrow(perm.mat.window50.cluster_i_bigger_window) + 1) / (permutation.nr + 1)
      p.mat.cluster_i[j, 1] <- as.numeric(cluster_i_w_p.val)
      
      cluster_i_w_p.val.adj <- p.adjust(cluster_i_w_p.val, method = "bonferroni", n = nr_of_tests)
      p.mat.cluster_i[j, 2] <- as.numeric(cluster_i_w_p.val.adj)
    }
    
    num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
    
    for (plot_idx in 1:num_plots) {
      start_rank <- (plot_idx - 1) * ranks_per_plot + 1
      end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
      
      features <- paste(disease_abbr, ot_gene_set_label, paste("Rank", start_rank:end_rank, sep = ""), "1", sep = "_")
      
      print(ST.FeaturePlot(se, 
                           features = features,
                           pt.alpha = spot_alpha,
                           pt.size = point_size,
                           grid.ncol = 4,
                           cols = c("black", "darkblue", "cyan", "yellow", "red", "darkred"),
                           value.scale = "all"))
    }
  } else {
    # Original logic for single cluster
    cluster_numbers <- cluster_number
    for (i in cluster_numbers) {
      subset_name <- paste0("se_", i)
      se_subset[[subset_name]] <- SubsetSTData(se, expression = seurat_clusters == i)
      meta_data_subset[[subset_name]] <- se_subset[[subset_name]][[]]
      
      start_with_abr_label <- grep(term, colnames(meta_data_subset[[subset_name]]), value = TRUE)
      meta_data_abr_label[[subset_name]] <- subset(meta_data_subset[[subset_name]], 
                                                   select = start_with_abr_label)
      
      abr_label_structure_i <- NULL
      for (k in 1:ncol(meta_data_abr_label[[subset_name]])) {
        meta_data_vector <- data.frame(
          Rank = as.character(k), 
          Median = median(meta_data_abr_label[[subset_name]][, k])
        )
        abr_label_structure_i <- rbind(abr_label_structure_i, meta_data_vector)
      }
      
      rownames(abr_label_structure_i) <- 1:ncol(meta_data_abr_label[[subset_name]])
      abr_label_structure_i <- as.data.frame(abr_label_structure_i)
      
      Rank_names <- paste0("Rank_nr", 1:length(window_rank_list))
      p.mat.cluster_i <- matrix(NA, nrow = length(window_rank_list), ncol = 2, 
                                dimnames = list(Rank_names, c("p.val", "p.val.adj")))
      
      permutation.nr <- dim(perm.mat)[1]
      nr_of_tests <- 1
      
      for (j in 1:length(window_rank_list)) {
        OT_genetic_window <- as.numeric(abr_label_structure_i[j, 2])
        cluster_name <- paste("cluster", i, sep = "_")
        
        perm.mat.window50.cluster_i_bigger_window <- subset(
          perm.mat.window50, perm.mat.window50[[cluster_name]] >= OT_genetic_window
        )
        
        cluster_i_w_p.val <- (nrow(perm.mat.window50.cluster_i_bigger_window) + 1) / (permutation.nr + 1)
        p.mat.cluster_i[j, 1] <- as.numeric(cluster_i_w_p.val)
        
        cluster_i_w_p.val.adj <- p.adjust(cluster_i_w_p.val, method = "bonferroni", n = nr_of_tests)
        p.mat.cluster_i[j, 2] <- as.numeric(cluster_i_w_p.val.adj)
      }
      
      num_plots <- ceiling(length(window_rank_list) / ranks_per_plot)
      
      for (plot_idx in 1:num_plots) {
        start_rank <- (plot_idx - 1) * ranks_per_plot + 1
        end_rank <- min(plot_idx * ranks_per_plot, length(window_rank_list))
        
        features <- paste(disease_abbr, ot_gene_set_label, paste("Rank", start_rank:end_rank, sep = ""), "1", sep = "_")
        
        print(ST.FeaturePlot(se_subset[[subset_name]], 
                             features = features,
                             pt.alpha = spot_alpha,
                             pt.size = point_size,
                             grid.ncol = 4,
                             cols = c("black", "darkblue", "cyan", "yellow", "red", "darkred"),
                             value.scale = "all"))
      }
    }
  }
}
