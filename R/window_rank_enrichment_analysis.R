#' Window Rank Enrichment Analysis
#'
#' Calculates enrichment statistics for each window rank across clusters,
#' comparing observed median scores to a null distribution from permutations.
#'
#' @param se A Seurat object containing gene scores and clustering metadata.
#' @param perm.mat.window.data A data frame of null median scores for ranked
#'   gene sets, with columns corresponding to cluster names.
#' @param window_rank_list A list or vector of window rank identifiers used
#'   to define the number and order of rank-specific gene sets.
#' @param ot_gene_set_label A string identifying the gene set type
#'   (e.g., "Genetic" or "Drugs").
#' @param disease_abbr A short abbreviation for the disease/trait
#'   (e.g., "SCZ" or "ALZ").
#' @param cluster_anno Column name in metadata specifying clusters
#'   (default: "seurat_clusters").
#' @param imputation Strategy for handling p-value calculation when no
#'   permutations are equal to or more extreme than the observed value.
#'   Options are:
#'   - `"all"`: Always add 1 to the numerator and denominator.
#'   - `"none"`: No imputation; p-values can be 0.
#'   - `"dynamic"`: Add 1 to the numerator and denominator only when no
#'     permutations are equal to or more extreme than the observed value.
#'
#' @returns A data frame containing for each cluster and window:
#' \describe{
#'   \item{cluster}{Cluster identifier.}
#'   \item{window}{Window rank number.}
#'   \item{median_score}{Observed median score for the window.}
#'   \item{p_value}{Nominal p-value for the window.}
#'   \item{q05}{5% quantile from null distribution.}
#'   \item{q95}{95% quantile from null distribution.}
#'   \item{is_significant}{Logical indicating if score is outside
#'     5-95% quantiles.}
#'   \item{was_imputed}{Logical indicating if p-value was imputed.}
#'   \item{imputation_type}{Type of imputation used.}
#' }
#'
#' @section Imputation Behavior:
#' For imputation = "none", if no median scores in the null distribution
#' are equal to or larger than the queried median score, the p-value is
#' set to 0. For imputation = "all", one is always added to the numerator
#' and denominator. For imputation = "dynamic", one is added only when
#' no null scores are equal to or larger than the observed score.
#'
#' @export
#'
#' @examples
#' window_results <- window_rank_enrichment_analysis(
#'   se = se,
#'   perm.mat.window.data = perm.mat.window.data,
#'   window_rank_list = window_rank_list_ALZ_Genetic,
#'   ot_gene_set_label = "Genetic",
#'   disease_abbr = "ALZ",
#'   cluster_anno = "cluster_anno",
#'   imputation = "dynamic"
#' )
window_rank_enrichment_analysis <- function(
    se,
    perm.mat.window.data,
    window_rank_list,
    ot_gene_set_label,
    disease_abbr,
    cluster_anno = "seurat_clusters",
    imputation = "all"
) {
  
  # Validate imputation parameter
  if (!imputation %in% c("all", "none", "dynamic")) {
    stop("imputation must be one of: 'all', 'none', 'dynamic'")
  }
  
  # Validate cluster annotation
  if (!(cluster_anno %in% colnames(se@meta.data))) {
    stop(
      paste(
        "The specified cluster column",
        cluster_anno,
        "is not found in the metadata."
      )
    )
  }
  
  # Check that window_rank_list is provided
  if (length(window_rank_list) == 0) {
    stop("window_rank_list must contain at least one window.")
  }
  
  # Find window-specific metadata columns
  pattern <- paste0(
    "^",
    disease_abbr,
    "_",
    ot_gene_set_label,
    "_Rank[0-9]+_1$"
  )
  
  all_window_cols <- grep(
    pattern,
    colnames(se@meta.data),
    value = TRUE
  )
  
  if (length(all_window_cols) == 0) {
    stop(
      paste(
        "No matching window score columns were found in the metadata for:",
        paste0(disease_abbr, "_", ot_gene_set_label)
      )
    )
  }
  
  # Extract rank numbers from column names
  window_ranks <- as.integer(
    sub(
      paste0(
        "^",
        disease_abbr,
        "_",
        ot_gene_set_label,
        "_Rank([0-9]+)_1$"
      ),
      "\\1",
      all_window_cols
    )
  )
  
  # Check that all ranks were successfully identified
  if (any(is.na(window_ranks))) {
    stop("Could not determine window rank numbers from the window score column names.")
  }
  
  # Sort window columns by rank
  all_window_cols <- all_window_cols[order(window_ranks)]
  window_ranks <- window_ranks[order(window_ranks)]
  
  # Check that the number of windows matches window_rank_list
  if (length(all_window_cols) != length(window_rank_list)) {
    stop(
      paste0(
        "The number of window score columns in the metadata (",
        length(all_window_cols),
        ") does not match the number of windows in window_rank_list (",
        length(window_rank_list),
        ")."
      )
    )
  }
  
  # Sort clusters using the package helper
  cluster_numbers <- sort_cluster_names(
    unique(se@meta.data[[cluster_anno]])
  )
  
  # Initialize results
  results <- list()
  result_index <- 1
  
  # Process each cluster
  for (cluster_id in cluster_numbers) {
    
    cluster_name <- as.character(cluster_id)
    cluster_key <- paste0("cluster_", cluster_name)
    
    # Check that cluster exists in permutation matrix
    if (!cluster_name %in% colnames(perm.mat.window.data)) {
      warning(
        paste(
          "Cluster",
          cluster_name,
          "not found in permutation matrix. Skipping."
        )
      )
      next
    }
    
    # Get cells belonging to cluster
    cells_in_cluster <- rownames(se@meta.data)[
      se@meta.data[[cluster_anno]] == cluster_id
    ]
    
    # Calculate observed median score for each window
    median_scores <- sapply(
      all_window_cols,
      function(col) {
        median(
          se@meta.data[cells_in_cluster, col],
          na.rm = TRUE
        )
      }
    )
    
    # Get null distribution for this cluster
    null_dist <- perm.mat.window.data[, cluster_name]
    null_dist <- null_dist[!is.na(null_dist)]
    
    if (length(null_dist) == 0) {
      warning(
        paste(
          "Null distribution for cluster",
          cluster_name,
          "contains only NA values. Skipping."
        )
      )
      next
    }
    
    # Calculate null distribution quantiles
    quantiles <- quantile(
      null_dist,
      probs = c(0.05, 0.95),
      na.rm = TRUE
    )
    
    q05 <- quantiles[1]
    q95 <- quantiles[2]
    
    # Process each window
    for (i in seq_along(median_scores)) {
      
      median_score <- median_scores[i]
      
      # Calculate number of null values equal to or more extreme
      # than the observed score
      n_more_extreme <- sum(
        null_dist >= median_score,
        na.rm = TRUE
      )
      
      # Calculate p-value
      if (is.na(median_score)) {
        
        p_value <- NA_real_
        was_imputed <- FALSE
        imputation_type <- "none"
        
      } else if (imputation == "all") {
        
        p_value <- (
          n_more_extreme + 1
        ) / (
          length(null_dist) + 1
        )
        
        was_imputed <- TRUE
        imputation_type <- "all"
        
      } else if (
        imputation == "dynamic" &&
        n_more_extreme == 0
      ) {
        
        p_value <- 1 / (length(null_dist) + 1)
        
        was_imputed <- TRUE
        imputation_type <- "dynamic"
        
      } else {
        
        p_value <- n_more_extreme / length(null_dist)
        
        was_imputed <- FALSE
        imputation_type <- "none"
      }
      
      # Determine whether observed score falls outside
      # the central 90% of the null distribution
      is_significant <- (
        !is.na(median_score) &&
          (
            median_score < q05 ||
              median_score > q95
          )
      )
      
      # Store result
      results[[result_index]] <- data.frame(
        cluster = cluster_key,
        window = window_ranks[i],
        median_score = median_score,
        p_value = p_value,
        q05 = q05,
        q95 = q95,
        is_significant = is_significant,
        was_imputed = was_imputed,
        imputation_type = imputation_type,
        stringsAsFactors = FALSE
      )
      
      result_index <- result_index + 1
    }
  }
  
  # Combine results
  if (length(results) == 0) {
    return(
      data.frame(
        cluster = character(),
        window = integer(),
        median_score = numeric(),
        p_value = numeric(),
        q05 = numeric(),
        q95 = numeric(),
        is_significant = logical(),
        was_imputed = logical(),
        imputation_type = character(),
        stringsAsFactors = FALSE
      )
    )
  }
  
  results <- do.call(rbind, results)
  
  # Set unique row names
  rownames(results) <- paste0(
    disease_abbr,
    "_",
    ot_gene_set_label,
    "_",
    results$cluster,
    "_Rank",
    results$window
  )
  
  # Convert cluster to factor to maintain cluster ordering
  results$cluster <- factor(
    results$cluster,
    levels = paste0("cluster_", cluster_numbers)
  )
  
  # Sort results by cluster and window
  results <- results[
    order(results$cluster, results$window),
    ,
    drop = FALSE
  ]
  
  return(results)
}
