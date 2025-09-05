#' Add Cluster Annotations to Seurat Object
#'
#' Adds cluster annotations from a data frame to the Seurat object's metadata.
#'
#' @param se A Seurat object
#' @param AnnoDataframe A data frame containing cluster annotations
#' @param BarcodeName Name of the column in AnnoDataframe containing cell barcodes
#' @param AnnoName Name of the column in AnnoDataframe containing cluster annotations
#' @param slot_name Name of the metadata slot to store the annotations (default: "cluster_anno")
#'
#' @returns Seurat object with added cluster annotations
#' @export
#'
#' @examples
#' seurat_clusters <- read.csv("path/to/clusters.csv", sep = ";")
#' se <- add_cluster_annotations(se, 
#'                            AnnoDataframe = seurat_clusters, 
#'                            BarcodeName = "Barcode", 
#'                            AnnoName = "seurat_clusters")
add_cluster_annotations <- function(se, 
                                AnnoDataframe, 
                                BarcodeName, 
                                AnnoName,
                                slot_name = "cluster_anno") {
    
    # Check if required columns exist
    if (!BarcodeName %in% colnames(AnnoDataframe)) {
        stop(paste("Column", BarcodeName, "not found in AnnoDataframe"))
    }
    if (!AnnoName %in% colnames(AnnoDataframe)) {
        stop(paste("Column", AnnoName, "not found in AnnoDataframe"))
    }
    
    # Create a named vector of cluster annotations
    cluster_anno <- setNames(AnnoDataframe[[AnnoName]], AnnoDataframe[[BarcodeName]])
    
    # Add to Seurat object metadata
    se[[slot_name]] <- cluster_anno
    
    return(se)
}
