findClusterMarkers <- function(se, 
                               nfactors, 
                               dims = NULL, 
                               M = 10) {
  if (is.null(dims)) {
    dims <- 1:nfactors
  }
  
  #NMF
  se <- RunNMF(se, nfactors = nfactors)
  
  #find neighbors and clusters
  se <- FindNeighbors(object = se, verbose = FALSE, reduction = "NMF", dims = dims)
  se <- FindClusters(object = se, verbose = FALSE)
  
  seurat_clusters <- as.data.frame(se$seurat_clusters)
  seurat_clusters$Barcode <- colnames(se)
  seurat_clusters$seurat_clusters <- seurat_clusters$`se$seurat_clusters`
  seurat_clusters$`se$seurat_clusters` <- NULL
  
  DefaultAssay(se) <- "RNA"
  
  # normalize and scale data
  se <- se %>%
    NormalizeData() %>%
    ScaleData()
  
  se <- SetIdent(se, value = "seurat_clusters")
  
  N.markers <- FindAllMarkers(se, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, slot = "data")
  
  #top markers
  topGenes <- N.markers %>% 
    group_by(cluster) %>% 
    arrange(p_val_adj) %>%
    slice_head(n = M)
  
  #UMAP
  se <- RunUMAP(se, reduction = "NMF", dims = dims, n.neighbors = 10)
  
  #marker genes for each cluster
  all_gene_sets_structure_markers <- lapply(0:(length(unique(Idents(se))) - 1), function(i) {
    top_genes <- as.character(t(subset(topGenes, cluster == i, select = gene)))
    return(top_genes)
  })
  
  names(all_gene_sets_structure_markers) <- paste0("cluster", 0:(length(unique(Idents(se))) - 1), "_marker_genes_")
  
  return(list(updated_se = se, marker_genes = all_gene_sets_structure_markers))
}

generate_permutation_matrix_parallel <- function(se, gnum = 50, permutation_nr = 10000, workers = 10, cluster_col = "seurat_clusters") {
  library(future.apply)
  library(dplyr)
  library(Seurat)
  
  options(future.globals.maxSize = 4 * 1024^3)
  plan(multisession, workers = workers)
  
  # Check if the specified cluster column exists
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  execution_time <- system.time({
    
    ref_genes_se <- rownames(se)
    genes.in.NULL.set <- gnum
    
    # Use the specified cluster column
    cluster_names <- unique(se@meta.data[[cluster_col]])
    Perm_names <- paste0("Perm_nr", seq_len(permutation_nr))
    perm.mat <- matrix(data = NA, nrow = permutation_nr, ncol = length(cluster_names), dimnames = list(Perm_names, cluster_names))
    
    perm_results <- future_lapply(seq_len(permutation_nr), function(i) {
      NULL_gene_set_ <- sample(ref_genes_se, genes.in.NULL.set)
      genes <- unlist(list(NULL_gene_set_))
      
      se_temp <- AddModuleScore(se, list(genes), ctrl = 100, name = "NULL_gene_set")
      # Access the dynamic cluster column name
      se_metadata_NULL <- se_temp@meta.data[, c(cluster_col, "NULL_gene_set1")]
      
      # Rename the cluster column for consistent processing
      colnames(se_metadata_NULL)[1] <- "cluster"
      
      se_metadata_by_cluster_NULL <- se_metadata_NULL %>%
        group_by(cluster) %>%
        summarise(median_score = median(NULL_gene_set1, na.rm = TRUE))
      
      return(se_metadata_by_cluster_NULL$median_score)
    }, future.seed = TRUE)
    
    perm.mat <- do.call(rbind, perm_results)
    rownames(perm.mat) <- Perm_names
    colnames(perm.mat) <- cluster_names
    perm.mat <- as.data.frame(perm.mat)
  })
  
  print(execution_time)
  return(perm.mat)
}


SpatialTraitEnrichmentAnalysis <- function(
    se,
    perm.mat.label.data,
    perm.mat.window50.data,
    window_rank_list_abr_label,
    ot_gene_set_label = "Genetic",
    disease_abbr = "ALZ",
    cluster_col = "seurat_clusters"
) {
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  se_subset <- list()
  meta_data_subset <- list()
  meta_data_abr_label <- list()
  
  # Check that clustering column exists
  if (!(cluster_col %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_col, "is not found in the Seurat object metadata."))
  }
  
  cluster_numbers <- sort(unique(se@meta.data[[cluster_col]]))
  
  for (i in cluster_numbers) {
    subset_name <- paste0("se_", i)
    
    cells_i <- rownames(se@meta.data)[se@meta.data[[cluster_col]] == i]
    se_subset[[subset_name]] <- subset(se, cells = cells_i)
    
    meta_data_subset[[subset_name]] <- se_subset[[subset_name]][[]]
    start_with_abr_label <- grep(term, colnames(meta_data_subset[[subset_name]]), value = TRUE)
    meta_data_abr_label[[subset_name]] <- subset(meta_data_subset[[subset_name]], select = start_with_abr_label)
  }
  
  permutation_nr <- dim(perm.mat.label.data)[1]
  nr_of_tests <- length(cluster_numbers)
  cluster_names <- paste0("cluster_", cluster_numbers)
  p_val_mat <- matrix(NA, nrow = 2, ncol = length(cluster_names),
                      dimnames = list(c("-log10(p.val)", "-log10(p.val.adj)"), cluster_names))
  
  for (index in seq_along(cluster_numbers)) {
    i <- cluster_numbers[index]
    subset_name <- paste0("se_", i)
    cluster_name <- paste0("cluster_", i)
    
    abr_label_structure_i <- data.frame()
    for (k in 1:length(meta_data_abr_label[[subset_name]])) {
      median_value <- median(meta_data_abr_label[[subset_name]][, k], na.rm = TRUE)
      abr_label_structure_i <- rbind(abr_label_structure_i, data.frame(Rank = as.character(k), Median = median_value))
    }
    
    for (j in seq_along(window_rank_list_abr_label)) {
      OT_label_window <- abr_label_structure_i[j, "Median"]
      
      perm_mat_window_cluster_bigger <- subset(perm.mat.window50.data,
                                               perm.mat.window50.data[[cluster_name]] >= OT_label_window)
      
      cluster_i_w_p_val <- (nrow(perm_mat_window_cluster_bigger) + 1) / (permutation_nr + 1)
      adjusted_p_val <- p.adjust(cluster_i_w_p_val, method = "bonferroni", n = nr_of_tests)
      
      ot_name <- paste("OpenTargets", disease_abbr, ot_gene_set_label, "1", sep = "_")
      perm_mat_label_bigger <- subset(perm.mat.label.data,
                                      perm.mat.label.data[[cluster_name]] >= median(as.numeric(unlist(se_subset[[subset_name]][[ot_name]]))))
      
      cluster_i_p_val <- (nrow(perm_mat_label_bigger) + 1) / (permutation_nr + 1)
      cluster_i_p_val_adj <- p.adjust(cluster_i_p_val, method = "bonferroni", n = nr_of_tests)
      
      p_val_mat[1, index] <- -log10(cluster_i_p_val)
      p_val_mat[2, index] <- -log10(cluster_i_p_val_adj)
    }
  }
  
  return(p_val_mat)
}


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


plotScoreDist <- function(se, 
                          perm.mat, 
                          perm.mat.window50, 
                          window_rank_list, 
                          ot_gene_set_label = "Genetic", 
                          disease_abbr = "ALZ",
                          cluster_col = "seurat_clusters",
                          use_cluster_prefix = TRUE) {
  
  theme_set(theme_classic())
  
  term <- paste0("^", disease_abbr, "_", ot_gene_set_label)
  cluster_labels <- unique(se@meta.data[[cluster_col]])
  
  for (cluster_label in cluster_labels) {
    message("\nProcessing cluster: ", cluster_label)
    
    # Filter cells by metadata manually
    cells_in_cluster <- colnames(se)[se@meta.data[[cluster_col]] == cluster_label]
    if (length(cells_in_cluster) == 0) {
      warning("No cells found for cluster: ", cluster_label)
      next
    }
    
    se_subset <- subset(se, cells = cells_in_cluster)
    meta_data_subset <- se_subset[[]]
    
    # Find all relevant OT enrichment columns
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
    
    # Determine cluster name used in permutation matrices
    cluster_name <- if (use_cluster_prefix) {
      paste0("cluster_", cluster_label)
    } else {
      as.character(cluster_label)
    }
    
    if (!(cluster_name %in% colnames(perm.mat))) {
      warning("Cluster name not found in permutation matrix: ", cluster_name)
      next
    }
    
    # Prepare for p-value computation
    Rank_names <- paste0("Rank_nr", 1:length(window_rank_list))
    p.mat.cluster_i <- matrix(NA, nrow = length(window_rank_list), ncol = 2,
                              dimnames = list(Rank_names, c("p.val", "p.val.adj")))
    
    permutation.nr <- nrow(perm.mat)
    nr_of_tests <- length(cluster_labels)
    
    for (j in 1:length(window_rank_list)) {
      OT_genetic_window <- abr_label_structure_i[j, "Median"]
      perm_bigger <- perm.mat.window50[perm.mat.window50[[cluster_name]] >= OT_genetic_window, ]
      
      p.raw <- (nrow(perm_bigger) + 1) / (permutation.nr + 1)
      p.adj <- p.adjust(p.raw, method = "bonferroni", n = nr_of_tests)
      
      p.mat.cluster_i[j, 1] <- p.raw
      p.mat.cluster_i[j, 2] <- p.adj
    }
    
    # Plot global distribution
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
    
    p <- ggplot(median_score_NULL, aes(x = Median_scores)) +
      geom_density(fill = "#69b3a2", alpha = 0.8) +
      geom_vline(xintercept = actual_median, color = "purple", size = 1) +
      annotate("text", x = 0, y = 0.01, label = annotation_text, hjust = 0) +
      ggtitle(paste(ot_gene_set_label, disease_abbr, "Structure", cluster_label)) +
      xlab("Median score NULL gene set") +
      ylab("Density")
    
    print(p)
  }
}

generate_metacell_clusters <- function(se, 
                                       cluster_col = "subclass_label",
                                       min_cells_per_type = 50,
                                       num_cells_per_metacell = 10,
                                       pca_dims = 10,
                                       initial_resolution = 0.1,
                                       verbose = FALSE) {
  
  # Filter cell types with too few cells
  label_counts <- table(se@meta.data[[cluster_col]])
  valid_labels <- names(label_counts[label_counts >= min_cells_per_type])
  se <- subset(se, cells = colnames(se)[se@meta.data[[cluster_col]] %in% valid_labels])
  
  # Get remaining cell types
  cell_types <- unique(se@meta.data[[cluster_col]])
  
  # Split by cell type
  cell_type_subsets <- lapply(cell_types, function(ct) {
    subset(se, cells = colnames(se)[se@meta.data[[cluster_col]] == ct])
  })
  
  # Preprocess and cluster into meta groups
  cell_type_clusters <- lapply(cell_type_subsets, function(ct_subset) {
    ct_subset <- SCTransform(ct_subset, verbose = verbose)
    ct_subset <- FindVariableFeatures(ct_subset, verbose = verbose)
    ct_subset <- RunPCA(ct_subset, verbose = verbose)
    ct_subset <- FindNeighbors(ct_subset, dims = 1:pca_dims)
    
    desired_clusters <- floor(ncol(ct_subset) / num_cells_per_metacell)
    resolution <- initial_resolution
    ct_subset <- FindClusters(ct_subset, resolution = resolution, verbose = verbose)
    num_clusters <- length(unique(ct_subset$seurat_clusters))
    
    # Adjust resolution
    while (num_clusters < desired_clusters) {
      resolution <- resolution + 0.1
      ct_subset <- FindClusters(ct_subset, resolution = resolution, verbose = verbose)
      num_clusters <- length(unique(ct_subset$seurat_clusters))
    }
    while (num_clusters > desired_clusters) {
      resolution <- resolution - 0.1
      ct_subset <- FindClusters(ct_subset, resolution = resolution, verbose = verbose)
      num_clusters <- length(unique(ct_subset$seurat_clusters))
    }
    return(ct_subset)
  })
  
  return(cell_type_clusters)
}

aggregate_metacells <- function(cell_type_clusters, cluster_col = "subclass_label") {
  
  aggregate_metacells_sum_with_metadata <- function(ct_subset) {
    clusters <- unique(ct_subset$seurat_clusters)
    meta_cell_data <- list()
    meta_cell_metadata <- list()
    
    for (cluster_id in clusters) {
      cluster_cells <- ct_subset[, ct_subset$seurat_clusters == cluster_id]
      summed_expression <- rowSums(cluster_cells@assays$RNA@counts)
      meta_cell_data[[as.character(cluster_id)]] <- summed_expression
      
      subclass_label_values <- cluster_cells@meta.data[[cluster_col]]
      most_common_ident <- names(sort(table(subclass_label_values), decreasing = TRUE))[1]
      meta_cell_metadata[[as.character(cluster_id)]] <- most_common_ident
    }
    
    meta_cell_matrix <- do.call(cbind, meta_cell_data)
    colnames(meta_cell_matrix) <- names(meta_cell_data)
    meta_cell_seurat <- CreateSeuratObject(counts = meta_cell_matrix)
    meta_cell_seurat@meta.data[[cluster_col]] <- factor(unlist(meta_cell_metadata), levels = unique(unlist(meta_cell_metadata)))
    
    return(meta_cell_seurat)
  }
  
  # Apply aggregation to all
  meta_cell_seurat_objects <- lapply(cell_type_clusters, aggregate_metacells_sum_with_metadata)
  merged_meta_cells <- Reduce(function(x, y) merge(x, y), meta_cell_seurat_objects)
  
  return(merged_meta_cells)
}

run_fisher_test_on_DEGs <- function(se, DEG_genes, OT_genes, adjust_method = "bonferroni", test_alternative = "greater") {
  
  total_genes <- rownames(se)
  contingency_results <- data.frame(cluster = character(), A = integer(), B = integer(), C = integer(), D = integer())
  
  for (cluster in unique(DEG_genes$cluster)) {
    cluster_data <- DEG_genes[DEG_genes$cluster == cluster, ]
    cluster_data$Genetically_Associated <- cluster_data$gene %in% OT_genes
    
    A <- sum(cluster_data$Genetically_Associated)
    C <- length(cluster_data$gene) - A
    not_DE_genes <- setdiff(total_genes, cluster_data$gene)
    B <- sum(not_DE_genes %in% OT_genes)
    D <- length(not_DE_genes) - B
    
    contingency_results <- rbind(contingency_results, data.frame(cluster = cluster, A = A, B = B, C = C, D = D))
  }
  
  fisher_results <- data.frame(cluster = character(), p_value = numeric(), odds_ratio = numeric())
  
  for (i in 1:nrow(contingency_results)) {
    A <- contingency_results$A[i]
    B <- contingency_results$B[i]
    C <- contingency_results$C[i]
    D <- contingency_results$D[i]
    
    fisher_matrix <- matrix(c(A, B, C, D), nrow = 2, byrow = TRUE)
    fisher_test <- fisher.test(fisher_matrix, alternative = test_alternative)
    
    fisher_results <- rbind(fisher_results, data.frame(
      cluster = contingency_results$cluster[i],
      p_value = fisher_test$p.value,
      odds_ratio = fisher_test$estimate
    ))
  }
  
  fisher_results$adj_p_value <- p.adjust(fisher_results$p_value, method = adjust_method)
  rownames(fisher_results) <- NULL
  
  return(fisher_results)
}

load_permutation_matrices <- function(file_paths) {
  matrices <- lapply(file_paths, readRDS)
  names(matrices) <- basename(file_paths)
  return(matrices)
}

calculate_gene_percentage <- function(se, pattern) {
  genes <- grep(pattern, rownames(se), value = TRUE)
  (Matrix::colSums(GetAssayData(se, slot = "counts")[genes, ]) /
      Matrix::colSums(GetAssayData(se, slot = "counts"))) * 100
}

