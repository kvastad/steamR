#' Spatial Trait Enrichment Analysis
#'
#' Performs permutation-based enrichment analysis of a trait-specific gene
#' score across spatial or single-cell clusters by comparing observed median
#' scores to a null distribution generated from random gene sets.
#'
#' @param se A Seurat object containing gene scores and clustering metadata.
#' @param perm.mat.label.data A data frame of null median scores generated
#'   from random gene lists of the same size as the gene list being tested.
#'   Columns should correspond to cluster names.
#' @param gene_list_score A character string matching the name of the gene
#'   score column in `se@meta.data` to test for enrichment.
#' @param cluster_anno Column name in `se@meta.data` specifying the clustering
#'   to use (default: `"seurat_clusters"`).
#' @param imputation Strategy for handling p-value calculation when no
#'   permutations are equal to or more extreme than the observed value.
#'   Options are:
#'   - `"all"`: Always add 1 to the numerator and denominator.
#'   - `"none"`: No imputation; p-values can be 0.
#'   - `"dynamic"`: Add 1 to the numerator and denominator only when no
#'     permutations are equal to or more extreme than the observed value.
#' @param perm.mat.window.data Deprecated. Previously used for window-level
#'   enrichment analysis. Retained for backward compatibility but ignored.
#'   Use `window_rank_enrichment_analysis()` for window-level analysis.
#' @param window_rank_list_abr_label Deprecated. Previously used to specify
#'   window rank identifiers. Retained for backward compatibility but ignored.
#'   Use `window_rank_enrichment_analysis()` for window-level analysis.
#'
#' @returns A list containing three data frames:
#'   \describe{
#'     \item{p_val_mat}{Nominal (unadjusted) p-values for each cluster.}
#'     \item{impute}{Logical matrix indicating whether the p-value was
#'       imputed.}
#'     \item{imputation_details}{Data frame containing cluster, p-value,
#'       imputation status, and imputation type.}
#'   }
#'
#' @section Permutation p-value:
#' For the one-sided enrichment test, null scores equal to or greater than
#' the observed score are considered as extreme or more extreme. The
#' permutation p-value is calculated according to the selected imputation
#' strategy.
#'
#' @section Deprecated arguments:
#' `perm.mat.window.data` and `window_rank_list_abr_label` are retained for
#' backward compatibility with versions <= 0.1.0 but are no longer used.
#' Window-level enrichment analysis is now performed separately with
#' `window_rank_enrichment_analysis()`. These arguments may be removed in
#' a future release.
#'
#' @export
#'
#' @examples
#' result <- spatial_trait_enrichment_analysis(
#'   se = se,
#'   perm.mat.label.data = perm.mat.genetic.data,
#'   gene_list_score = "OpenTargets_SCZ_Genetic_1",
#'   cluster_anno = "supercluster_term",
#'   imputation = "dynamic"
#' )
spatial_trait_enrichment_analysis <- function(
    se,
    perm.mat.label.data,
    gene_list_score,
    cluster_anno = "seurat_clusters",
    imputation = "all",
    perm.mat.window.data = NULL,
    window_rank_list_abr_label = NULL
) {
  
  # Deprecated arguments retained for backward compatibility
  if (!is.null(perm.mat.window.data)) {
    warning(
      "`perm.mat.window.data` is deprecated and is no longer used. ",
      "Window enrichment analysis is now performed separately using ",
      "`window_rank_enrichment_analysis()`."
    )
  }
  
  if (!is.null(window_rank_list_abr_label)) {
    warning(
      "`window_rank_list_abr_label` is deprecated and is no longer used. ",
      "Window enrichment analysis is now performed separately using ",
      "`window_rank_enrichment_analysis()`."
    )
  }
  
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
        "is not found in the Seurat object metadata."
      )
    )
  }
  
  # Validate gene score
  if (!(gene_list_score %in% colnames(se@meta.data))) {
    stop(
      paste(
        "The specified gene list column",
        gene_list_score,
        "is not found in the metadata."
      )
    )
  }
  
  # Sort cluster names using the package helper
  cluster_numbers <- sort_cluster_names(
    unique(se@meta.data[[cluster_anno]])
  )
  
  cluster_names <- paste0("cluster_", cluster_numbers)
  
  # Initialize result matrices
  p_val_mat <- matrix(
    NA_real_,
    nrow = 1,
    ncol = length(cluster_names),
    dimnames = list("p.val", cluster_names)
  )
  
  imputed_mat <- matrix(
    FALSE,
    nrow = 1,
    ncol = length(cluster_names),
    dimnames = list("imputed", cluster_names)
  )
  
  # Process each cluster
  for (index in seq_along(cluster_numbers)) {
    
    cluster_id <- cluster_numbers[index]
    cluster_name <- as.character(cluster_id)
    
    # Check that cluster is represented in permutation matrix
    if (!cluster_name %in% colnames(perm.mat.label.data)) {
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
    
    # Calculate observed median
    observed_score <- median(
      se@meta.data[cells_in_cluster, gene_list_score],
      na.rm = TRUE
    )
    
    # Get null distribution
    null_dist <- perm.mat.label.data[, cluster_name]
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
    
    # Calculate number of null values equal to or more extreme
    # than the observed score
    n_more_extreme <- sum(
      null_dist >= observed_score,
      na.rm = TRUE
    )
    
    # Calculate p-value according to imputation strategy
    if (is.na(observed_score)) {
      
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
      
    } else if (imputation == "dynamic" && n_more_extreme == 0) {
      
      p_value <- 1 / (length(null_dist) + 1)
      
      was_imputed <- TRUE
      imputation_type <- "dynamic"
      
    } else {
      
      p_value <- n_more_extreme / length(null_dist)
      
      was_imputed <- FALSE
      imputation_type <- "none"
    }
    
    # Store results
    p_val_mat[1, index] <- p_value
    imputed_mat[1, index] <- was_imputed
  }
  
  # Combine results
  result <- list(
    p_val_mat = as.data.frame(p_val_mat),
    impute = as.data.frame(imputed_mat),
    imputation_details = data.frame(
      cluster = cluster_names,
      p_value = as.numeric(p_val_mat[1, ]),
      was_imputed = as.logical(imputed_mat[1, ]),
      imputation_type = ifelse(
        imputed_mat[1, ],
        imputation,
        "none"
      ),
      stringsAsFactors = FALSE
    )
  )
  
  return(result)
}
