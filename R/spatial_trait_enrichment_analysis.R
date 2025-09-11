#' Spatial Trait Enrichment Analysis
#'
#' Performs permutation-based enrichment analysis of a trait-specific gene score 
#' (e.g. from OpenTargets) across spatial or single-cell clusters, comparing 
#' observed medians to a null distribution.
#'
#' @param se A Seurat object containing gene scores and clustering metadata.
#' @param perm.mat.label.data A data frame of null median scores generated from random gene lists of the same size as the gene list being tested, where columns correspond to cluster names.
#' @param perm.mat.window.data A data frame of null median scores for ranked gene sets (e.g., sliding windows), with columns corresponding to cluster names.
#' @param window_rank_list_abr_label A list or vector of window rank identifiers 
#'        used to reference the rank-specific gene sets.
#' @param gene_list_score A character string matching the name of the gene score column in `se@meta.data` to test for enrichment (e.g., "OpenTargets_SCZ_Genetic_1").
#' @param cluster_anno Column name in `se@meta.data` specifying the clustering to use 
#'        (e.g., "seurat_clusters" or "supercluster_term").
#' @param imputation Strategy for handling p-value calculation when no permutations exceed the observed value.
#'        Options are:
#'        - "all": Always add 1 to numerator and denominator (default)
#'        - "none": No imputation, can result in p=0
#'        - "dynamic": Only impute when no permutations exceed observed value
#' @param log_file Path to the log file for imputed p-values
#'
#' @returns A list containing three data frames:
#'          - "p_val_mat": Nominal (unadjusted) p-values for each cluster
#'          - "impute": Information about imputed p-values
#'          - "imputation_details": A list of data frames, one per cluster, containing detailed imputation metrics:
#'            - n_more_extreme: Number of null median scores that are as extreme or more extreme than the observed median score
#'            - n_more_extreme_full: Total number of null median scores (used for p-value calculation)
#'            - imputed: Logical indicating whether the p-value was imputed
#'            - p-value: The calculated nominal p-value
#'
#' @section Debug Log:
#' The log file (if provided) contains detailed debug information for each cluster and window:
#' - OT_label_window: The observed median score for the window
#' - Null distribution summary: Summary statistics (Min, 1st Qu, Median, Mean, 3rd Qu, Max) of the null distribution
#' - n_more_extreme: Number of null median scores that are as extreme or more extreme than the observed median score
#' - n_more_extreme_full: Total number of null median scores (used for p-value calculation)
#' - imputed: Logical indicating whether the p-value was imputed
#' - p-value: The calculated nominal p-value
#'
#' @export
#'
#' @examples
#' pval_mat <- spatial_trait_enrichment_analysis(
#'   se = se,
#'   perm.mat.label.data = perm.mat.genetic.data,
#'   perm.mat.window.data = perm.mat.window.data,
#'   window_rank_list_abr_label = window_rank_list_SCZ_Genetic,
#'   gene_list_score = "OpenTargets_SCZ_Genetic_1",
#'   cluster_anno = "supercluster_term",
#'   imputation = "dynamic"
#' )
spatial_trait_enrichment_analysis <- function(
    se,
    perm.mat.label.data,
    perm.mat.window.data,
    window_rank_list_abr_label,
    gene_list_score,
    cluster_anno = "seurat_clusters",
    imputation = "all",
    log_file = NULL
) {
  # Validate imputation parameter
  if (!imputation %in% c("all", "none", "dynamic")) {
    stop("imputation must be one of: 'all', 'none', 'dynamic'")
  }
  
  # Check that clustering column exists
  if (!(cluster_anno %in% colnames(se@meta.data))) {
    stop(paste("The specified cluster column", cluster_anno, "is not found in the Seurat object metadata."))
  }
  
  # Check gene list column exists
  if (!(gene_list_score %in% colnames(se@meta.data))) {
    stop(paste("The specified gene list column", gene_list_score, "is not found in the metadata."))
  }
  
  cluster_numbers <- sort(unique(se@meta.data[[cluster_anno]]))
  cluster_names <- paste0("cluster_", cluster_numbers)
  
  # Initialize results matrix
  permutation_nr <- dim(perm.mat.label.data)[1]
  nr_of_tests <- length(cluster_numbers)
  p_val_mat <- matrix(NA, nrow = 1, ncol = length(cluster_names),
                      dimnames = list("p.val", cluster_names))
  
  # Initialize matrix to track imputed p-values
  imputed_mat <- matrix(FALSE, nrow = 1, ncol = length(cluster_names),
                       dimnames = list("imputed", cluster_names))
  
  # Open a connection to the log file if it's provided
  if (!is.null(log_file)) {
    log_con <- file(log_file, open = "wt")
    on.exit(close(log_con), add = TRUE)
  }
  
  # Process each cluster
  for (index in seq_along(cluster_numbers)) {
    i <- cluster_numbers[index]
    cluster_name <- as.character(i)
    
    if (!cluster_name %in% colnames(perm.mat.label.data)) {
      warning(paste("Cluster", cluster_name, "not found in permutation matrix. Skipping."))
      next
    }
    
    # Get cells for this cluster
    cells_i <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] == i]
    
    # Calculate median for the gene list in this cluster
    actual_median <- median(se@meta.data[cells_i, gene_list_score], na.rm = TRUE)
    
    # Get matched columns for this cluster
    matched_cols <- grep(gene_list_score, colnames(se@meta.data), value = TRUE)
    
    if (length(matched_cols) == 0) {
      warning(paste("No matching columns found for:", gene_list_score, "in cluster:", i))
      next
    }
    
    # Calculate medians for matched columns
    cluster_medians <- sapply(matched_cols, function(col) {
      median(se@meta.data[cells_i, col], na.rm = TRUE)
    })
    names(cluster_medians) <- matched_cols
    
    # Process windows
    for (j in seq_along(window_rank_list_abr_label)) {
      # Use the actual median value for this window
      OT_label_window <- cluster_medians[1]
      
      # Skip individual window logging - will log cluster summary at the end
      
      # Get null distribution for this cluster and window
      if (!cluster_name %in% colnames(perm.mat.window.data)) {
        warning(paste("Cluster", cluster_name, "not found in permutation matrix. Skipping."))
        next
      }
      
      null_dist <- perm.mat.window.data[, cluster_name]
      
      # Remove NA values from null distribution
      null_dist <- null_dist[!is.na(null_dist)]
      if (length(null_dist) == 0) {
        warning(paste("Null distribution for cluster", cluster_name, "contains only NA values. Skipping."))
        next
      }
      
      # Calculate number of more extreme values
      n_more_extreme <- sum(null_dist > OT_label_window, na.rm = TRUE)
      
      # Only log warnings if there's an issue
      if (n_more_extreme == 0 || is.na(OT_label_window)) {
        if (imputation == "none") {
          warning_msg <- sprintf("No median scores in the null distribution were larger than the queried median score for cluster %s window %d. Consider using imputation 'all' or 'dynamic'.", cluster_name, j)
          if (!is.null(log_file)) {
            writeLines(warning_msg, log_con)
          } else {
            warning(warning_msg)
          }
        }
        
        if (imputation == "all") {
          cluster_i_w_p_val <- 1 / (length(null_dist) + 1)
          imputed_mat[1, index] <- TRUE
        } else if (imputation == "dynamic") {
          cluster_i_w_p_val <- 1 / (length(null_dist) + 1)
          imputed_mat[1, index] <- TRUE
        } else { # imputation == "none"
          cluster_i_w_p_val <- 0
          imputed_mat[1, index] <- FALSE
        }
      } else {
        if (imputation == "all") {
          cluster_i_w_p_val <- (n_more_extreme + 1) / (length(null_dist) + 1)
          imputed_mat[1, index] <- TRUE
        } else {
          cluster_i_w_p_val <- n_more_extreme / length(null_dist)
          imputed_mat[1, index] <- FALSE
        }
      }
      
      # Get null distribution for the full gene list
      if (!cluster_name %in% colnames(perm.mat.label.data)) {
        warning(paste("Cluster", cluster_name, "not found in permutation matrix. Skipping."))
        next
      }
      null_dist_full <- perm.mat.label.data[, cluster_name]
      
      # Remove NA values from null distribution
      null_dist_full <- null_dist_full[!is.na(null_dist_full)]
      if (length(null_dist_full) == 0) {
        warning(paste("Null distribution for cluster", cluster_name, "contains only NA values. Skipping."))
        next
      }
      
      # Calculate number of more extreme values for full gene list
      n_more_extreme_full <- sum(null_dist_full > actual_median, na.rm = TRUE)
      
      # Calculate p-value with specified imputation strategy
      if (n_more_extreme_full == 0) {
        if (imputation == "none") {
          warning_msg <- sprintf("No median scores in the null distribution were larger than the queried median score for cluster %s full gene list. Consider using imputation 'all' or 'dynamic'.", cluster_name)
          if (!is.null(log_file)) {
            writeLines(warning_msg, log_con)
          } else {
            warning(warning_msg)
          }
        }
        
        if (imputation == "all") {
          cluster_i_p_val <- 1 / (length(null_dist_full) + 1)
          imputed_mat[1, index] <- TRUE
        } else if (imputation == "dynamic") {
          cluster_i_p_val <- 1 / (length(null_dist_full) + 1)
          imputed_mat[1, index] <- TRUE
        } else { # imputation == "none"
          cluster_i_p_val <- 0
          imputed_mat[1, index] <- FALSE
        }
      } else {
        if (imputation == "all") {
          cluster_i_p_val <- (n_more_extreme_full + 1) / (length(null_dist_full) + 1)
          imputed_mat[1, index] <- TRUE
        } else {
          cluster_i_p_val <- n_more_extreme_full / length(null_dist_full)
          imputed_mat[1, index] <- FALSE
        }
      }
      
      # Add cluster-level summary to log file
      if (!is.null(log_file)) {
        cluster_summary <- sprintf("\nCluster: %s\n", cluster_name)
        cluster_summary <- paste0(cluster_summary, sprintf("Gene set: %s\n", gene_list_score))
        cluster_summary <- paste0(cluster_summary, sprintf("Median score: %.3f\n", actual_median))
        cluster_summary <- paste0(cluster_summary, sprintf("P-value: %.6f\n", cluster_i_p_val))
        cluster_summary <- paste0(cluster_summary, sprintf("Imputation: %s\n", imputation))
        cluster_summary <- paste0(cluster_summary, sprintf("Was imputed: %s\n", imputed_mat[1, index]))
        cluster_summary <- paste0(cluster_summary, sprintf("Number of windows: %d\n", length(window_rank_list_abr_label)))
        
        writeLines(cluster_summary, log_con)
      }
      
      p_val_mat[1, index] <- cluster_i_p_val
    }
  }
  
  # Combine p-values and imputation information
  result <- list(
    p_val_mat = as.data.frame(p_val_mat),
    impute = as.data.frame(imputed_mat),
    imputation_details = data.frame(
      cluster = cluster_names,
      p_value = as.numeric(p_val_mat[1,]),
      was_imputed = as.logical(imputed_mat[1,]),
      imputation_type = ifelse(imputed_mat[1,], imputation, "none"),
      stringsAsFactors = FALSE
    )
  )
  
  return(result)
}