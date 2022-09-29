
<!-- README.md is generated from README.Rmd. Please edit that file -->

# sgolay

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/sgolay)](https://CRAN.R-project.org/package=sgolay)
[![R-CMD-check](https://github.com/zeehio/sgolay/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zeehio/sgolay/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/zeehio/sgolay/branch/main/graph/badge.svg)](https://app.codecov.io/gh/zeehio/sgolay?branch=main)
<!-- badges: end -->

The goal of sgolay is to offer efficient and vectorized Savitzky-Golay
filters.

## Installation

You can install the CRAN version with

``` r
# install.packages("sgolay")
```

Or you can install the development version of sgolay from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("zeehio/sgolay")
```

## Benchmark

`sgolay` is faster than `signal`, especially on either larger filter
lengths or when applied on matrices, since it uses the Fast Fourier
Transform and avoids several memory copies and extra allocations.

``` r
library(sgolay)
x <- matrix(runif(1000*1000), nrow = 1000, ncol = 1000)

filt <- signal::sgolay(p = 2, n = 51)

timing <- bench::mark(
  signal = apply(x, 2L, function(s) signal::sgolayfilt(s, filt)),
  sgolay = sgolay::sgolayfilt(x, filt), 
  min_iterations = 50
)
#> Warning: Some expressions had a GC in every iteration; so filtering is disabled.
plot(timing, type = 'ridge')
#> Loading required namespace: tidyr
#> Picking joint bandwidth of 0.0135
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />
