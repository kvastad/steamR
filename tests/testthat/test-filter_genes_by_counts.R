test_that("filter_genes_by_counts() filters genes correctly", {
  # Generate a mock sparse matrix for testing
  set.seed(123)
  mock_data <- Matrix::Matrix(
    matrix(rpois(500, lambda = c(rep(2, 40), rep(50, 10))), nrow = 50, ncol = 10),
    sparse = TRUE
  )
  rownames(mock_data) <- paste0("Gene", 1:50)
  colnames(mock_data) <- paste0("Spot", 1:10)
  
  # Create a Seurat object
  se_mock <- Seurat::CreateSeuratObject(counts = mock_data)
  
  # Apply the filtering function
  se_filtered <- filter_genes_by_counts(se_mock, minCounts = 20, minSpots = 3)
  
  # Tests
  expect_s4_class(se_filtered, "Seurat") # Should return a Seurat object
  expect_lt(nrow(se_filtered), nrow(se_mock)) # Filtered object should have fewer genes
  expect_true(all(Matrix::rowSums(Seurat::GetAssayData(se_filtered, layer = "counts")) >= 20)) # Check count filter
  expect_true(all(Matrix::rowSums(Seurat::GetAssayData(se_filtered, layer = "counts") > 0) >= 3)) # Check spots filter
})
