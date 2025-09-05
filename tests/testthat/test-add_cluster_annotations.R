test_that("add_cluster_annotations correctly adds annotations to Seurat object", {
  # Load required packages
  library(Seurat)
  library(steam)
  
  #mock Seurat object
  mock_counts <- matrix(
    data = rpois(50, lambda = 5),
    nrow = 10,
    ncol = 5,
    dimnames = list(paste0("Gene", 1:10), paste0("Cell", 1:5))
  )
  se_mock <- CreateSeuratObject(counts = mock_counts)
  
  #mock annotation dataframe
  annotations <- data.frame(
    Barcode = colnames(se_mock),
    seurat_clusters = sample(0:2, 5, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  se_annotated <- add_cluster_annotations(
    se = se_mock,
    AnnoDataframe = annotations,
    BarcodeName = "Barcode",
    AnnoName = "seurat_clusters"
  )
  
  #chck if metadata column was added
  expect_true("seurat_clusters" %in% colnames(se_annotated@meta.data))
  
  #check that cluster values match
  expect_equal(
    se_annotated$seurat_clusters,
    setNames(annotations$seurat_clusters, annotations$Barcode)[colnames(se_annotated)]
  )
})
