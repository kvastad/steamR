#' Window Rank Enrichment Analysis
#'
#' Calculates enrichment statistics for each window rank across clusters,
#' comparing observed median scores to a null distribution from permutations.
#'
#' @param se A Seurat object containing gene scores and clustering metadata.
#' @param perm.mat.window50.data A data frame of null median scores for ranked gene sets
#'        (e.g., sliding windows), with columns corresponding to cluster names.
#' @param window_rank_list A list or vector of window rank identifiers 
#'        used to reference the rank-specific gene sets.
#' @param ot_gene_set_label A string identifying the gene set type (e.g., "Genetic" or "Drugs").
#' @param disease_abbr A short abbreviation for the disease/trait (e.g., "SCZ" or "ALZ").
#' @param cluster_anno Column name in metadata specifying clusters (default: "seurat_clusters").
#' @param imputation Strategy for handling p-value calculation when no permutations exceed the observed value.
#'        Options are:
#'        - "all": Always add 1 to numerator and denominator (default)
#'        - "none": No imputation, can result in p=0
#'        - "dynamic": Only impute when no permutations exceed observed value
#'
#' @returns A data frame containing for each cluster and window:
#' \describe{
#'   \item{cluster}{Cluster identifier}
#'   \item{window}{Window rank number}
#'   \item{observed_score}{Observed median score for the window}
#'   \item{p_value}{Nominal p-value for the window}
#'   \item{q05}{5% quantile from null distribution}
#'   \item{q95}{95% quantile from null distribution}
#'   \item{is_significant}{Logical indicating if score is outside 5-95% quantiles}
#' }
#'
#' @export
#'
#' @examples
#' window_results <- WindowRankEnrichmentAnalysis(
#'   se = se,
#'   perm.mat.window50.data = perm.mat.window50.data,
#'   window_rank_list = window50_rank_list_ALZ_Drugs,
#'   ot_gene_set_label = "Drugs",
#'   disease_abbr = "ALZ",
#'   imputation = "dynamic"
#' )
WindowRankEnrichmentAnalysis <- function(
    se,
    perm.mat.window50.data,
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
    
    # Input validation
    if (!(cluster_anno %in% colnames(se@meta.data))) {
        stop(paste("The specified cluster column", cluster_anno, "is not found in the metadata."))
    }
    
    # Initialize results data frame
    results <- data.frame()
    
    # Get unique clusters
    cluster_numbers <- unique(se@meta.data[[cluster_anno]])
    
    # Try to convert clusters to numeric, handling non-numeric gracefully
    numeric_clusters <- suppressWarnings(as.numeric(cluster_numbers))
    numeric_indices <- !is.na(numeric_clusters)
    
    # Split clusters into numeric and non-numeric
    numeric_clusters <- sort(numeric_clusters[numeric_indices])
    non_numeric_clusters <- sort(cluster_numbers[!numeric_indices])
    
    # Combine them in the correct order
    sorted_clusters <- c(numeric_clusters, non_numeric_clusters)
    
    # Pattern for matching window-specific columns
    pattern <- paste0("^", disease_abbr, "_", ot_gene_set_label, "_Rank")
    
    # Get all window-specific columns once
    all_window_cols <- grep(pattern, colnames(se@meta.data), value = TRUE)
    if (length(all_window_cols) == 0) {
        stop(paste("No matching columns found for pattern:", pattern))
    }
    
    # Process each cluster
    for (cluster_id in sorted_clusters) {
        cluster_name <- as.character(cluster_id)
        cluster_key <- paste0("cluster_", cluster_name)
        
        # Get cells for this cluster
        cells_in_cluster <- rownames(se@meta.data)[se@meta.data[[cluster_anno]] == cluster_id]
        
        # Calculate observed median scores for each window using the metadata directly
        observed_scores <- sapply(all_window_cols, function(col) {
            median(se@meta.data[cells_in_cluster, col], na.rm = TRUE)
        })
        
        # Get null distribution for this cluster
        null_dist <- perm.mat.window50.data[[cluster_name]]
        if (is.null(null_dist)) {
            warning(paste("No null distribution found for cluster:", cluster_name))
            next
        }
        
        # Calculate quantiles once for this cluster
        quantiles <- quantile(null_dist, probs = c(0.05, 0.95), na.rm = TRUE)
        q05 <- quantiles[1]
        q95 <- quantiles[2]
        
        # Calculate p-values and create results for this cluster
        for (i in seq_along(observed_scores)) {
            observed <- observed_scores[i]
            
            # Calculate number of more extreme values
            n_more_extreme <- sum(null_dist >= observed)
            
            # Calculate p-value with specified imputation strategy
            if (imputation == "all") {
                p_value <- (n_more_extreme + 1) / (length(null_dist) + 1)
            } else if (imputation == "none") {
                p_value <- n_more_extreme / length(null_dist)
            } else { # dynamic
                if (n_more_extreme == 0) {
                    p_value <- 1 / (length(null_dist) + 1)
                    message(sprintf("Imputed p-value for cluster %s window %d: %.2e", 
                                  cluster_name, i, p_value))
                } else {
                    p_value <- n_more_extreme / length(null_dist)
                }
            }
            
            # Check if score is significant (outside 5-95% quantiles)
            is_significant <- !is.na(observed) && 
                            (observed < q05 || observed > q95)
            
            # Add to results
            results <- rbind(results, data.frame(
                cluster = cluster_key,
                window = i,
                observed_score = observed,
                p_value = p_value,
                q05 = q05,
                q95 = q95,
                is_significant = is_significant,
                stringsAsFactors = FALSE
            ))
        }
    }
    
    # Set unique row names by combining cluster and window information
    rownames(results) <- paste0(disease_abbr, "_", ot_gene_set_label, "_", results$cluster, "_Rank", results$window)
    
    # Convert cluster to factor with levels in the correct order to maintain sorting
    results$cluster <- factor(results$cluster, 
                            levels = paste0("cluster_", sorted_clusters))
    
    # Sort the results by cluster (now a factor) and window
    results <- results[order(results$cluster, results$window), ]
    
    return(results)
} 