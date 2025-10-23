
# STEAM R-package

<!-- badges: start -->
<!-- badges: end -->

STEAM is an algorithm for Spatial Trait Enrichment Analysis with perMutation testing, a robust computational approach to measure the enrichment of average gene expression across clusters in a dataset from a given gene list; it calculates a permutation p-value and performs multiple testing corrections based on the number of clusters. For ranked gene lists, STEAM enables interrogation of the topmost relevant sets of genes based on their combined average enrichment.

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

## Vignettes

Vignettes can be found here [link](https://github.com/kvastad/STEAM).
