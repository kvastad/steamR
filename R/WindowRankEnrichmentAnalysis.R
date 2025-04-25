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
#' @param cluster_col Column name in metadata specifying clusters (default: "seurat_clusters").
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
#'   disease_abbr = "ALZ"
#' )
WindowRankEnrichmentAnalysis <- function(
    se,
    perm.mat.window50.data,
    window_rank_list,
    ot_gene_set_label,
    disease_abbr,
    cluster_col = "seurat_clusters"
) {
    # Input validation
    if (!(cluster_col %in% colnames(se@meta.data))) {
        stop(paste("The specified cluster column", cluster_col, "is not found in the metadata."))
    }
    
    # Initialize results data frame
    results <- data.frame()
    
    # Get unique clusters
    cluster_numbers <- sort(unique(se@meta.data[[cluster_col]]))
    
    # Pattern for matching window-specific columns
    pattern <- paste0("^", disease_abbr, "_", ot_gene_set_label, "_Rank")
    
    # Process each cluster
    for (cluster_id in cluster_numbers) {
        cluster_name <- as.character(cluster_id)
        cluster_key <- paste0("cluster_", cluster_name)
        
        # Subset cells for this cluster
        cells_in_cluster <- colnames(se)[se@meta.data[[cluster_col]] == cluster_id]
        se_subset <- subset(se, cells = cells_in_cluster)
        
        # Get window-specific columns
        window_cols <- grep(pattern, colnames(se_subset@meta.data), value = TRUE)
        
        if (length(window_cols) == 0) {
            warning(paste("No matching columns found for pattern:", pattern, "in cluster:", cluster_id))
            next
        }
        
        # Calculate observed median scores for each window
        observed_scores <- sapply(window_cols, function(col) {
            median(se_subset@meta.data[[col]], na.rm = TRUE)
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
            
            # Calculate p-value
            n_more_extreme <- sum(null_dist >= observed)
            p_value <- (n_more_extreme + 1) / (length(null_dist) + 1)
            
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