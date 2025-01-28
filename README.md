
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `scaleGMRF`

<!-- badges: start -->
<!-- badges: end -->

`scaleGMRF` is an R package that provides useful functions to standardize Gaussian Markov Random Fields (GMRF) effects in Latent Gaussian Models (LGM) as proposed by  [Ferrari & Ventrucci (2025)](https://arxiv.org/abs/2501.16057).

Standardizing GMRF effects guarantees that their variance parameters
match their intuitive interpretation, defined as the variance
contribution of the effects as intended by the user: this greatly
simplifies the task of prior specification of the variance parameters of
LGMs.

A working paper about the new tools implemented in this package is under
development.

## Installation

The `scaleGMRF` package can be installed from
[GitHub](https://github.com/):

``` r
# install.packages("devtools")
devtools::install_github("LFerrariIt/scaleGMRF")
```

## Context

A Latent Gaussian model (Hrafnkelsson, 2023) is defined as a Bayesian
Hierarchical model in which the linear predictor contains one or more
Gaussian effects for different covariates (Hrafnkelsson, 2023). A
Gaussian effect for a given covariate can be specified by the choice of:

1.  a basis, i.e. a set of known functions of the covariate;

2.  a set of coefficients following a multivariate Gaussian
    distribution, with null mean and known precision matrix up to a
    scalar variance parameter;

3.  a distribution for the covariate.

The effect evaluated at different values of the covariate also follows a
multivariate Gaussian distribution. Another name for multivariate
Gaussian distributions is Gaussian Markov Random Fields (GMRF), which
have been thoroughly studied in Rue & Held, 2005.

In this context, we define “standardization” a procedure that ensures
that the variance parameters of GMRFs match quantities that are
intuitive for the users, namely the variance contributions of the
effects. This is useful as it greatly simplifies the process of prior
specification of the variance parameters, as prior beliefs about the
variance contribution of an effects can be directly reflected through a
prior on the corresponding variance parameters. The standardization
procedure that this package proposes differs for fixed and random
effects, as two different definitions of “variance contribution of an
effect” are adopted for the two categories.

## Overview

### Main function: `standardize_GMRF()`

The main function of the package is `standardize_GMRF()`, which takes as
arguments a GMRF effect (via its 3 components) and its categorization as
either fixed or random. `standardize_GMRF()` implements standardization
and returns a standardized version of the GMRF as a list of elements.
The function is built such that the elements of its output can be
directly used to specify an effect of an `INLA` model using the
`inla.stack()` function and the model `"generic0"` (www.r-inla.org, Rue
et al. 2009). If the argument `plot_check=TRUE`, the function also
prints a graphical output that can be used to visually assess whether
the standardization procedure has been successful.

### User-friendly wrapper: `f_Xunif()`

Another useful function of the package is `f_Xunif()`, which is a
user-friendly wrapper of `standardize_GMRF()` that returns standardized
versions of popular GMRF effects, under the convenient assumption that
the covariate follows a Uniform distribution. The function takes as
arguments simply a sample of values of the covariate and a string
indicating the type of effect, e.g. `"linear","iid","besag"`, etc. 

<!--### Modified P-Splines
&#10;An important class of effects implemented in the `f_Xunif()` function are the P-Spline effects, which are popularly used in LGMs (Fahrmeir et al. 2004). Applying the standardization procedure to these effects require a slight modification of the precision matrices traditionally used for their specification. The motivation and the design of this modified version of P-Splines is presented in `vignette("psplines", package = "scaleGMRF")`.---->

### Worked example

The folder “Leuk_application” contains a worked example on the usage of
the package on a real dataset: the functions of `scaleGMRF` are used to
standardize the effects of a LGM model designed for the INLA dataset
`Leuk`.
<!--The folder contains both a commented `.R` script file and a `.Rmd` walk-through version.---->

## Example of usage

The standardization of a first-order random walk effect and its
consequent inclusion in a simple Gaussian model is reported here to
illustrate the usage of `scaleGMRF` and its integrability within the
INLA framework.

``` r
rm(list = ls())
library(INLA)
library(scaleGMRF)
# Number of values of a discrete covariate X
K <- 50
# Values of the covariate as an ordered factor
x <- factor(1:K, ordered = TRUE)
# Standardized random walk effect of order 1
standard_rw1 <- f_Xunif(x, model = "rw1", fixed = TRUE, plot_check = T)
summary(standard_rw1)
# Realizations of a response Y
Y <- cumsum(rnorm(K)) + rnorm(K,sd = 2)
# INLA data stack
data_stack <- inla.stack(
  data = list(Y = Y), # response
  A = list(1, standard_rw1$basis), # basis matrices
  effects = list(
    data.frame(intercept = rep(1, length(Y))),
    list(u_rw1 = 1:K) # coefficients
  )
)
# INLA model
model <- inla(
  Y ~ -1 + intercept +
    f(u_rw1,
      model = "generic0", constr = FALSE,
      Cmatrix = standard_rw1$precision, # precision matrix
      extraconstr = # constraints
        list(
          A = t(standard_rw1$null_space),
          e = matrix(rep(0, ncol(standard_rw1$null_space)))
        )
    ),
  family = "gaussian", data = inla.stack.data(data_stack),
  control.predictor = list(A = inla.stack.A(data_stack))
)
```

<!-- ----# Plot of response, real random walk process, and estimated one
plot(as.numeric(x),Y)
lines(as.numeric(x),real_rw1_process)
lines(as.numeric(x),
      model$summary.fixed$mean+model$summary.random$u_rw1$mean,col=2) -->

## Getting help

If you need help or find any bug, feel free to send an email to
<luisa.ferrari5@unibo.it>.

## Acknowledgments

This package has been developed by [Luisa
Ferrari](https://www.unibo.it/sitoweb/luisa.ferrari5/en) and
[Prof. Massimo
Ventrucci](https://www.unibo.it/sitoweb/massimo.ventrucci/en) from the
University of Bologna (Italy).

The project has been supported by the research grant “PRIN 2022 - COD.
2022PA3BS2 - META 2 - METAbarcoding for METAcommunities: towards a
genetic approach to community ecology”
(<https://www.metasquared.unito.it/>), in collaboration with the
Department of Economics of the University of Modena & Reggio Emilia
(Italy) and the Department of Life Sciences and Systems Biology of the
University of Turin (Italy).

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
