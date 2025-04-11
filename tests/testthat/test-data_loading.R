test_that("load_spatial_data() validates files and loads data correctly", {
  # Resolve the absolute path to the infoTable
  infoTablePath <- system.file("extdata", "infoTable.csv", package = "STEAM")
  
  # Ensure the function runs without error
  expect_no_error({
    se <- load_spatial_data(infoTablePath)
    
    # Check if the returned object is of class "Seurat"
    expect_s4_class(se, "Seurat")
    
    # Check if the embedded Staffli object exists
    expect_true("Staffli" %in% names(se@tools))
    expect_s4_class(se@tools$Staffli, "Staffli")
    
    # Validate data dimensions
    expect_gt(ncol(se), 0)  # Ensure there are spots
    expect_gt(nrow(se), 0)  # Ensure there are genes
  })
})
