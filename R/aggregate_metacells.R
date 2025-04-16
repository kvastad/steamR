#' Aggregate Gene Expression into Metacells
#'
#' Aggregates expression data across clusters (metacells) by summing counts
#' within each cluster of each cell type. Preserves the most common cell type label as metadata.
#'
#' @param cell_type_clusters A list of clustered Seurat objects, typically output from `generate_metacell_clusters()`.
#' @param cluster_col The metadata column used to track original cell type annotations (default: `"subclass_label"`).
#'
#' @returns A merged Seurat object representing all metacells with summed gene expression and assigned metadata.
#' @export
#'
#' @examples
#' clustered_se_list <- generate_metacell_clusters(se, cluster_col = "subclass_label")
#' metacell_seurat <- aggregate_metacells(clustered_se_list)
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