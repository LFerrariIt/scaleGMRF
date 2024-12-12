
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scaleGMRF

<!-- badges: start -->
<!-- badges: end -->

scaleGMRF is a package that provides useful functions to standardize
Gaussian Markov Random Fields (GMRF) effects used in Latent Gaussian
models.

The reason to standardize GMRF effects is to guarantee that their scale
parameters match their intuitive interpretation, defined as the variance
contribution of the effects as intended by the user. This is important
to correctly reflect prior information about the variance contributions
of the different effects of Latent Gaussian models through prior
specification on the scale parameters.

The package also contains functions to specify a modified version of
P-Spline effects (either in one or two dimensions).

A working paper about the new tools implemented in this package is under
development.

## Installation

You can install the scaleGMRF package from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("LFerrariIt/scaleGMRF")
```

## Overview

A Latent Gaussian model (Hrafnkelsson, 2023) is defined as a model in
which the linear predictor contains one or more Gaussian Markov Random
Field (GMRF) effects for different covariates $X_1,X_2,...$
(Hrafnkelsson, 2023, Rue & Held, 2005). A GMRF effect is defined in this
package as: $$f(X)=\mathbf{D}(X)\mathbf{u}$$
$$\mathbf{u}|\sigma^2 \sim N(\mathbf{0},\sigma^2 \mathbf{Q}^{-1})$$
$$ X\sim \pi(x)$$

$X$ is the covariate of the effect, $\mathbf{D}(X)$ is called the basis
and contains $K$ basis functions of $X$, $\mathbf{u}$ is called the set
of coefficients, $\mathbf{Q}$ is called the precision matrix, and
$\sigma^2$ is called the variance parameter. The effect evaluated at
different values of the covariate follows a Multivariate Normal
distribution.

The package contains a function `r_GMRF()` that generates realizations
of GMRF effects over a given sample $\mathbf{x}$ through the
specification of $\mathbf{Q}$,$\mathbf{D}(X)$.

The main function of the package is `standardize_GMRF()`, which takes as
argument a GMRF effect via the specification of
$\mathbf{Q}$,$\mathbf{D}(X)$,$\pi(X)$ and returns its standardized
version as a list of elements. The function is built such that the
elements of its output can be directly used to specify an effect of an
`INLA` model using the `inla.stack()` and the model `"generic0"`
(www.r-inla.org, Rue et al. 2009).

The concept of “standardization” and the usage of `standardize_GMRF()`
are thoroughly discussed in
`vignette("standardization",package="scaleGMRF")`.

Another useful function of the package is `f_Xunif()`, which is a
user-friendly wrapper of `standardize_GMRF()`. The function takes as
argument simply a string indicating the common name used to describe the
effect (e.g. `"linear","iid","besag"`, etc.) and it returns a
standardized version of many popular effects under the convenient
assumption that the covariate $X$ has a Uniform distribution. The
effects implemented in `f_Xunif()` are listed and exemplified in
`vignette("f_Xunif",package="scaleGMRF")`.

Setting the argument `plot_check=TRUE` in both `standardize_GMRF()` and
`f_Xunif()` causes the functions to print plots that can be used to
visually check whether the standardization procedure has been
successful.

An important class of effects implemented in the `f_Xunif()` function
are the P-Spline effects (Fahrmeir et al. 2004). Applying the
standardization procedure to these effects require a slight modification
of the precision matrices traditionally used for their specification.
The motivation and the design of this modified version of P-Splines is
presented in `vignette("psplines",package="scaleGMRF")`.

## Example of usage with `INLA`

## Getting help

If you need help or find any bug, feel free to send an email to
<luisa.ferrari5@unibo.it>.

## References

- Hrafnkelsson, B. (2023). *Statistical Modeling Using Bayesian Latent
  Gaussian Models.* Springer International Publishing.

- Rue, H., & Held, L. (2005). *Gaussian Markov random fields: theory and
  applications*. Chapman and Hall/CRC.

- Rue, H., Martino, S., & Chopin, N. (2009). Approximate Bayesian
  inference for latent Gaussian models by using integrated nested
  Laplace approximations. *Journal of the Royal Statistical Society
  Series B: Statistical Methodology*, 71(2), 319-392.

- Fahrmeir, L., Kneib, T., & Lang, S. (2004). Penalized structured
  additive regression for space-time data: a Bayesian perspective.
  *Statistica Sinica*, 731-761.
