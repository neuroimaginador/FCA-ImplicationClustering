
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FCA-ImplicationClustering

<!-- badges: start -->

![R](https://img.shields.io/badge/R-%23f0f0ff.svg?&style=flat&logo=r&logoColor=2065ba)
<!-- badges: end -->

The goal of FCA-ImplicationClustering is to present examples of how to
perform clustering on implications. This is intended to be a repository
with code in the `R` language to let users of the `fcaR` package
replicate the results of a paper submitted to ICFCA 2021.

## Setup

Some packages, besides `fcaR`, are needed. Users can install them by
doing:

``` r
install.packages(c("tidyverse", "cluster", "ggrepel", "MASS"))
```

## Examples

We have included two examples in this repository:

-   `planets_example.R`, and the accompanying `.md` and `.html` files,
    present a simple example on how to perform implication clustering
    with a small dataset.
-   `replication_script.R` and the accompanying `.md` and `.html` files
    are the code and results of clustering the implications valid in the
    MONKS datasets. It replicates the results presented in the paper
    mentioned above.

If you want to check the results from this Github repository, you can
inspect the `.md` files directly. Otherwise, you can clone or download
the repository and check the `.html` files in your computer.
