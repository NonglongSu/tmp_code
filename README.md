
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
# load and plot fasta sample
sample_data = system.file("extdata/samples", "1.fa", package = "alnDotPlot")
plot_dot_matrix(sample_data)

# load and plot json sample using ggplot
sample_data = system.file("extdata", "sample.json", package = "alnDotPlot")
plot_dot_matrix(input = sample_data, use_ggplot = TRUE)
```

