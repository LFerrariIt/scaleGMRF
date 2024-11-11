
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scaleGMRF

<!-- badges: start -->
<!-- badges: end -->

scaleGMRF is a package that provides useful function to standardize
Gaussian Markov Random Fields effects used in Latent Gaussian models.
The standardization guarantees that the corresponding scale parameters
of the effects match their intuitive interpretation defined as the
variance contribution of the effect as intended by the user.
Additionally, the package also contains function to define the modified
version of P-Spline effects: specifically, unidimensional P-Spline
effects of first- and second-order, and two-dimensional P-Splines of
first-order.

## Installation

You can install the development version of scaleGMRF from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("LFerrariIt/scaleGMRF")
```

## Main function: standardize_GMRF()

???
