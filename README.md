
# alnDotPlot

<!-- badges: start -->
<!-- badges: end -->

The goal of alnDotPlot is to generate and plot dot matrices of pairwise sequence alignments.

## Installation

You can install the development version of alnDotPlot from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jgarciamesa/alnDotPlot")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(alnDotPlot)
# load sample data
sample_data = system.file("extdata", "1.fa". package = alnDotPlot)
plot_dot_matrix(sample_data)
```

