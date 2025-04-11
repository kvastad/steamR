test_that("permutation_null_matrix works correctly", {
  # Mock inputs
  info_table_path <- "/inst/extdata/infoTable.csv"
  seurat_clusters_path <- "/inst/extdata/seurat_clusters.csv"
  gnum_path <- "/inst/extdata/UPPMAX_runs_of_permutations/gnum.txt"
  output_dir <- tempdir()  # Use a temporary directory for testing output
  idx <- 3
  num_permutations <- 10  # Keep small for testing purposes
  workers <- 2
  
  # Run the function
  result <- permutation_null_matrix(
    info_table_path = info_table_path,
    seurat_clusters_path = seurat_clusters_path,
    gnum_path = gnum_path,
    output_dir = output_dir,
    idx = idx,
    num_permutations = num_permutations,
    workers = workers
  )
  
  # Check the structure of the result
  expect_true(is.list(result))
  expect_true("perm.mat" %in% names(result))
  expect_true("execution_time" %in% names(result))
  
  # Check that the permutation matrix has the correct dimensions
  expect_equal(nrow(result$perm.mat), num_permutations)
  expect_true(ncol(result$perm.mat) > 0)
  
  # Check if the output files were saved correctly
  rds_file <- file.path(output_dir, paste0("perm.mat_p", num_permutations, "_g", idx, ".RDS"))
  csv_file <- file.path(output_dir, paste0("perm.mat_p", num_permutations, "_g", idx, ".csv"))
  expect_true(file.exists(rds_file))
  expect_true(file.exists(csv_file))
  
  # Clean up
  unlink(rds_file)
  unlink(csv_file)
})
