#' Run Fisher's Exact Test on Differentially Expressed Genes
#'
#' Performs Fisher's exact test to evaluate enrichment of a gene set
#' (e.g., genetically associated genes) among differentially expressed
#' genes (DEGs) for each cluster.
#'
#' @param se A Seurat object containing the full gene expression dataset.
#' @param DEG_genes A data frame of differentially expressed genes with
#'   columns `gene` and `cluster`.
#' @param OT_genes A character vector of genes of interest
#'   (e.g., OpenTargets genetic hits).
#' @param adjust_method Method for p-value adjustment (default:
#'   `"bonferroni"`). Passed to `p.adjust()`.
#' @param test_alternative Alternative hypothesis for Fisher's exact test
#'   (default: `"greater"`). Passed to `fisher.test()`.
#'
#' @returns A data frame containing the following columns:
#'   \describe{
#'     \item{cluster}{Cluster identifier.}
#'     \item{A}{Number of DEGs that are also in `OT_genes`.}
#'     \item{B}{Number of DEGs that are not in `OT_genes`.}
#'     \item{C}{Number of non-DEGs that are in `OT_genes`.}
#'     \item{D}{Number of non-DEGs that are not in `OT_genes`.}
#'     \item{prop_DEG}{Proportion of DEGs that are in `OT_genes`.}
#'     \item{prop_background}{Proportion of non-DEGs that are in
#'       `OT_genes`.}
#'     \item{odds_ratio}{Odds ratio estimated by Fisher's exact test.}
#'     \item{log2_OR}{Base-2 logarithm of the odds ratio.}
#'     \item{CI_lower}{Lower limit of the 95 percent confidence interval for the odds ratio.}
#'     \item{CI_upper}{Upper limit of the 95 percent confidence interval for the odds ratio.}
#'     \item{p_value}{Nominal p-value from Fisher's exact test.}
#'     \item{adj_p_value}{P-value adjusted using `adjust_method`.}
#'   }
#'
#' @export
#'
#' @examples
#' fisher_results <- run_fisher_test_degs(
#'   se,
#'   DEG_genes = N.markers,
#'   OT_genes = OpenTargets_SCZ_Genetic
#' )
run_fisher_test_degs <- function(
    se,
    DEG_genes,
    OT_genes,
    adjust_method = "bonferroni",
    test_alternative = "greater"
) {
  
  total_genes <- rownames(se)
  
  # Prepare list to store results
  fisher_results_list <- list()
  
  # Process each cluster
  for (cluster in unique(DEG_genes$cluster)) {
    
    cluster_data <- DEG_genes[DEG_genes$cluster == cluster, ]
    cluster_data$Genetically_Associated <- cluster_data$gene %in% OT_genes
    
    # Counts
    A <- sum(cluster_data$Genetically_Associated) # DEG & OT
    B <- length(cluster_data$gene) - A             # DEG & not OT
    
    not_DE_genes <- setdiff(total_genes, cluster_data$gene)
    
    C <- sum(not_DE_genes %in% OT_genes)           # non-DEG & OT
    D <- length(not_DE_genes) - C                  # non-DEG & not OT
    
    # 2x2 contingency table
    fisher_matrix <- matrix(
      c(A, B, C, D),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(
        c("DEG", "Not_DEG"),
        c("OT_gene", "Not_OT_gene")
      )
    )
    
    # Fisher's exact test
    fisher_test <- fisher.test(
      fisher_matrix,
      alternative = test_alternative
    )
    
    # Proportions
    prop_DEG <- A / (A + B)
    prop_background <- C / (C + D)
    
    # Odds ratio from Fisher's exact test
    odds_ratio <- unname(fisher_test$estimate)
    log2_OR <- log2(odds_ratio)
    
    # Store results
    fisher_results_list[[as.character(cluster)]] <- data.frame(
      cluster = cluster,
      A = A,
      B = B,
      C = C,
      D = D,
      prop_DEG = prop_DEG,
      prop_background = prop_background,
      odds_ratio = odds_ratio,
      log2_OR = log2_OR,
      CI_lower = fisher_test$conf.int[1],
      CI_upper = fisher_test$conf.int[2],
      p_value = fisher_test$p.value,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  fisher_results <- do.call(
    rbind,
    fisher_results_list
  )
  
  # Adjust p-values
  fisher_results$adj_p_value <- p.adjust(
    fisher_results$p_value,
    method = adjust_method
  )
  
  # Reset row names
  rownames(fisher_results) <- NULL
  
  return(fisher_results)
}
