#' Generate Metacell Clusters by Cell Type
#'
#' Preprocesses and clusters cells within each cell type into metacells of approximately equal size
#' using PCA and clustering on the specified number of dimensions.
#'
#' @param se A Seurat object containing single-cell data.
#' @param cluster_anno The metadata column specifying cell types to split on (default: `"subclass_label"`).
#' @param min_cells_per_type Minimum number of cells required for a cell type to be included (default: 50).
#' @param num_cells_per_metacell Approximate number of cells per metacell cluster (default: 10).
#' @param pca_dims Number of PCA dimensions to use for clustering (default: 10).
#' @param initial_resolution Starting resolution for clustering (default: 0.1). Will be adjusted to reach desired cluster number.
#' @param verbose Logical, whether to print Seurat messages (default: `FALSE`).
#'
#' @returns A list of Seurat objects, each corresponding to a clustered cell type.
#' @export
#'
#' @examples
#' metacells <- generate_metacell_clusters(
#'   se = my_seurat_object,
#'   cluster_anno = "subclass_label",
#'   num_cells_per_metacell = 10
#' )
generate_metacell_clusters <- function(se, 
                                       cluster_anno = "subclass_label",
                                       min_cells_per_type = 50,
                                       num_cells_per_metacell = 10,
                                       pca_dims = 10,
                                       initial_resolution = 0.1,
                                       verbose = FALSE) {
  
  label_counts <- table(se@meta.data[[cluster_anno]])
  valid_labels <- names(label_counts[label_counts >= min_cells_per_type])
  removed_labels <- setdiff(names(label_counts), valid_labels)
  print("Removed cell types with too few cells:")
  print(removed_labels)
  se <- subset(se, cells = colnames(se)[se@meta.data[[cluster_anno]] %in% valid_labels])
  
  cell_types <- unique(se@meta.data[[cluster_anno]])
  
  cell_type_subsets <- lapply(cell_types, function(ct) {
    subset(se, cells = colnames(se)[se@meta.data[[cluster_anno]] == ct])
  })
  
  cell_type_clusters <- lapply(cell_type_subsets, function(ct_subset) {
    cell_type <- unique(ct_subset@meta.data[[cluster_anno]])[1]
    message("\nProcessing cell type: ", cell_type)
    
    ct_subset <- SCTransform(ct_subset, verbose = verbose)
    ct_subset <- FindVariableFeatures(ct_subset, verbose = verbose)
    ct_subset <- RunPCA(ct_subset, verbose = verbose)
    ct_subset <- FindNeighbors(ct_subset, dims = 1:pca_dims)
    
    desired_clusters <- floor(ncol(ct_subset) / num_cells_per_metacell)
    resolution <- initial_resolution
    ct_subset <- FindClusters(ct_subset, resolution = resolution, verbose = verbose)
    num_clusters <- length(unique(ct_subset$seurat_clusters))
    
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
