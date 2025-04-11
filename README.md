<<<<<<< HEAD
# steam
STEAM is an algorithm for Spatial Trait Enrichment Analysis with perMutation testing, a robust computational approach to measure the enrichment of average gene expression across clusters in a dataset from a given gene list.
=======
# STEAM Package

The STEAM package is designed for Spatial Trait Enrichment Analysis with Permutation testing. Below is the directory structure for the package.


## Components

- **DESCRIPTION**: Package metadata, dependencies, and author information.
- **NAMESPACE**: Lists the functions and methods that are exported by the package.
- **R/**: Contains R scripts with core functionality.
  - **data_loading.R**: Functions for loading and preparing data.
  - **spatial_structures.R**: Functions for defining and subsetting spatial structures.
  - **opentarget_gene_sets.R**: Functions for managing and filtering gene sets.
  - **permutation_matrices.R**: Functions for calculating permutation matrices.
  - **module_scores.R**: Functions for calculating module scores.
  - **p_value_calculations.R**: Functions for p-value calculation and adjustment.
  - **visualization.R**: Functions for generating plots and visualizations.
  - **utilities.R**: Helper functions (e.g., `set_ident`, `add_annotations`).
- **man/**: Contains documentation for each function.
- **vignettes/**: Contains example workflows, including `STEAM_workflow.Rmd`.
- **tests/**: Includes test scripts to verify package functionality.
  - **testthat/**: Organized tests for each R script.
  - **testthat.R**: Entry point for running all tests in the package.

## Current Functions Overview

### Data Loading and Preparation

- **`load_spatial_data()`**: Loads and filters ST data from input tables, allowing specification of minimum counts and spots for genes.
- **`load_images()`**: Loads images associated with ST data for overlay and visualization.
- **`filter_genes_by_counts()`**: Filters genes by minimum counts and spots to ensure data quality.
- **`calculate_qc_metrics()`**: Calculates quality control metrics such as mitochondrial and ribosomal content percentages.

### Spatial Structure Definition

- **`define_spatial_structures()`**: Defines spatial structures (e.g., clusters) based on clustering or other segmentation criteria.
- **`subset_by_structure()`**: Subsets data by spatial structures for focused analysis of specific regions or clusters.

### Gene Set Management

- **`load_opentarget_gene_sets()`**: Loads and processes OpenTargets gene sets, including full and sliding-window sets.
- **`filter_gene_sets_by_st_data()`**: Filters gene sets to retain only genes present in the ST data, ensuring dataset relevance.
- **`select_top_genes()`**: Selects the top genes from OpenTargets gene sets for refined analysis.

#### Future work
- Functions for downloading and processing gene sets from OpenTargets and other sources, such as MSigDB or Gene Ontology.

### Permutation Matrix Calculation

- **`generate_permutation_matrix()`**: Calculates permutation matrices by simulating null distributions for each gene set across spatial structures. **Parallel processing** is supported for efficiency.
- **`calculate_null_scores()`**: Generates null distributions for each gene set through random sampling.
- **`save_permutation_matrix()`**: Saves calculated permutation matrices for future reuse, improving efficiency in downstream analyses.

### Permutation Testing

- **`load_permutation_matrices()`**: Loads saved permutation matrices to retrieve null distributions for Alzheimer’s gene set analysis.
- **`filter_observed_genes()`**: Filters Alzheimer’s-related gene sets for comparison with permutation matrices.
- **`calculate_observed_scores()`**: Calculates observed scores for Alzheimer’s gene sets, to be compared with null distributions for enrichment testing.
- **`compare_observed_to_null()`**: Compares observed scores to null distributions, calculating p-values to assess statistical significance.

### Add Module Scores

- **`calculate_module_scores()`**: Calculates module scores for each gene set across spatial structures using `AddModuleScore` to assess enrichment.

### P-Value Calculation and Adjustment

- **`calculate_p_values()`**: Computes p-values by comparing observed gene set scores to permutation null distributions.
- **`adjust_p_values()`**: Applies Bonferroni correction or other multiple testing adjustments to p-values.

### Visualization Functions

- **`plot_feature_distribution()`**: Plots distributions of selected features, such as mitochondrial and ribosomal content, across spatial clusters.
- **`plot_violin_plots()`**: Creates violin plots for gene set scores across spatial clusters, allowing comparison of expression distributions.
- **`plot_spatial_clusters()`**: Visualizes spatial clusters on the ST image, overlaying cluster information on the tissue image.
- **`feature_overlay_plot()`**: Overlays gene set scores on the spatial image for visual inspection of spatial enrichment.
- **`plot_observed_vs_null_distribution()`**: Visualizes observed gene set scores against null distributions to assess enrichment significance.
- **`plot_adjusted_p_values()`**: Displays adjusted p-values for observed gene set scores, highlighting statistically significant findings.
- **`plot_permutation_results()`**: Plots results from permutation testing for spatial structures, allowing insight into spatial enrichment.

### Utility Functions

- **`set_ident()`**: Sets identification levels in the data, such as cluster IDs, for use in downstream analysis.
- **`add_annotations()`**: Adds custom annotations (e.g., `seurat_clusters`) to the data for enhanced visualization and filtering.

>>>>>>> f9ed8e6 (Initial commit)
