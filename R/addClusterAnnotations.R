#' Add Cluster Annotations to Seurat Object
#'
#' This function extracts barcode and annotation columns from a provided dataframe
#' and adds them directly to the metadata of a Seurat object under the column `seurat_clusters`.
#' It ensures that annotations are correctly aligned by barcode.
#'
#' @param se A Seurat object.
#' @param AnnoDataframe A dataframe containing cell or spatial barcodes and their corresponding annotations.
#' @param BarcodeName Name of the column in `AnnoDataframe` containing barcodes (must match `colnames(se)`).
#' @param AnnoName Name of the column in `AnnoDataframe` containing cluster annotations.
#'
#' @return A Seurat object with added metadata column `seurat_clusters`.
#' @export
#'
#' @examples
#' seurat_clusters <- read.csv("data/annotations/seurat_clusters.csv", sep = ";")
#' se <- addClusterAnnotations(se, AnnoDataframe = seurat_clusters, BarcodeName = "Barcode", AnnoName = "seurat_clusters")
addClusterAnnotations <- function(se,
                                  AnnoDataframe,
                                  BarcodeName,
                                  AnnoName) {
  seurat_clusters <- AnnoDataframe[, c(BarcodeName, AnnoName)]
  colnames(seurat_clusters) <- c("Barcode", "seurat_clusters")
  rownames(seurat_clusters) <- seurat_clusters$Barcode
  seurat_clusters$Barcode <- NULL
  se <- AddMetaData(se, metadata = seurat_clusters)
  return(se)
}
