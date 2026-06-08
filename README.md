
# STEAM R-package

<!-- badges: start -->
<!-- badges: end -->

<div style="text-align: right;"> 
  <img src="images/STEAM_logo.png" alt="STEAM logo" width="100"/> 
</div>

STEAM is a framework for **S**patial **T**rait **E**nrichment **A**nalysis with per**M**utation testing, a robust computational approach for measuring the enrichment of average gene expression across clusters in a dataset from a given gene list. It calculates a permutation p-value and performs multiple-testing corrections based on the number of clusters. For ranked gene lists, STEAM enables interrogation of the topmost relevant sets of genes based on their combined average enrichment. The STEAM framework also includes an approach using MetaSpots (or MetaCells or MetaBins) within assigned clusters for sparse and high-resolution data. Using permutations to guide the selection of topmost relevant ranked genes and then testing for trait gene set enrichment among differentially expressed genes.

Here, we provide the source code for the R-package. If you encounter any issues with the R-package, please report those here.

## Installation

You can install the latest version of STEAM with:

``` r
# install.packages("devtools")
devtools::install_github("kvastad/steamR")
```

Once installed, load the package in R:

``` r
library(STEAM)
```

To install the latest stable STEAM R-package version:

``` r
devtools::install_github("kvastad/steamR@v0.1.0")
```


## Vignettes

For STEAM R-package vignettes, visit: https://kvastad.github.io/steamR/

## Trait gene lists

For ranked genetically associated gene lists, visit the Open Targets platform: https://platform.opentargets.org/

