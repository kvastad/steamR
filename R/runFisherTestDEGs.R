#' Run Fisher's Exact Test on Differentially Expressed Genes
#'
#' Performs Fisher's exact test to evaluate enrichment of a gene set (e.g., genetically associated genes)
#' among differentially expressed genes (DEGs) for each cluster.
#'
#' @param se A Seurat object containing the full gene expression dataset.
#' @param DEG_genes A data frame of differentially expressed genes with at least columns `gene` and `cluster`.
#' @param OT_genes A character vector of genes of interest (e.g., OpenTargets genetic hits).
#' @param adjust_method Method for p-value adjustment (default: `"bonferroni"`). Passed to `p.adjust()`.
#' @param test_alternative Alternative hypothesis for Fisher's test (default: `"greater"`).
#'
#' @returns A data frame with columns: `cluster`, `p_value`, `odds_ratio`, and `adj_p_value`.
#' @export
#'
#' @examples
#' fisher_results <- runFisherTestDEGs(se, DEG_genes = N.markers, OT_genes = OpenTargets_SCZ_Genetic)
runFisherTestDEGs <- function(se, DEG_genes, OT_genes, adjust_method = "bonferroni", test_alternative = "greater") {
  
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