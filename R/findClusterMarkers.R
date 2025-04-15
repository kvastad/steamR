#' Find cluster markers
#'
#' @param se A Seurat object
#' @param nfactors Number of factors to use for NMF
#' @param dims Dimensions to use for clustering
#' @param M Number of top markers to return for each cluster
#'
#' @returns A list containing marker genes for each cluster
#' @export 
#'
#' @examples
#' clusterMarkers <- FindClusterMarkers(se, nfactors = 10, dims = 1:10, M = 10)
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