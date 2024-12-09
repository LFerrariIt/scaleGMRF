#' Standardized effects for a Uniform covariate
#'
#' `f_Xunif()` returns a list to build various standardized effects (specified through the `model` argument) under the assumption that the corresponding covaiate follows a Uniform distribution (either, discrete or continuous). The first three elements of the list are designed to be used in INLA to specify an effect in its standardized version through the model `generic0`. The list also reports the scaling constant, a sample from the Uniform distribution on the covariate, and  the basis matrix evaluated at this sample: these elements can be useful for the user. The argument `plot_check=TRUE` returns a plot showing realizations of the standardized effect along with the histograms of the mean and variance of each realization.

#' @param x A numeric vector (for `model="linear","pspline1","pspline2","pspline_2D"`) or an ordered factor with more than 1 level (for `model="iid","rw1","rw2","besag"``).
#' @param model String to select the type of effect between `"iid","rw1","rw2","besag","linear","pspline1","pspline2","pspline_2D"`.
#' @param K A positive integer larger than 3, specifying the number of basis functions. Required for `model="pspline1","pspline2","pspline_2D"`.
#' @param adj_mat Matrix required for `model="besag"`, representing the adjacency matrix from the lattice. See Details.
#' @param m,M A pair of numbers indicating respectively the lower and upper boundaries of the support of X. If not provided, it is set to the extreme values from `x`. Required for `model="linear","pspline1","pspline2","pspline_2D"`.
#' @param fixed Logical, indicating whether the effect must be treated as fixed or random. By default, `TRUE`.
#' @param scale_Q Logical, indicating whether the precision matrix to be returned must already be scaled or not. By default, `TRUE`.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed or not. By default, `FALSE`.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot. By default, `n_sim=1000`.
#'
#'
#' @return A list of 6 elements, containing:
#' * `precision`: precision matrix
#' * `basis`:  basis matrix evaluated at `x`
#' * `scaling_constant`: a positive number, representing the appropriate scaling constant
#'  * `null_space`: a matrix, representing the null space of the precision matrix
#'  * `X_distribution`: a numeric vector (or a matrix of 2 columns) of values sampled from the Uniform distribution of X
#'  * `basis_distribution`: basis matrix evaluated at `X_distribution`#'
#' @examples
#' # Examples for discrete covariates --------------------------
#' K <- 20 # number of levels
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K)) # sample from X as an ordered factor object
#' f_Xunif(x, model = "iid", plot_check = TRUE) # fixed iid effect
#' f_Xunif(x, model = "iid", fixed = FALSE, plot_check = TRUE) # random iid effect
#' f_Xunif(x, model = "rw1", plot_check = TRUE) # random walk of order 1
#' f_Xunif(x, model = "rw2", plot_check = TRUE) # random walk of order 2
#'
#' # Examples for continuous covariates --------------------------
#' x <- runif(1000, min = 2, max = 5)
#' f_Xunif(x, model = "linear", plot_check = TRUE) # linear effect
#' K <- 20 # number of b-spline basis functions
#' f_Xunif(x, model = "pspline1", K = K, plot_check = TRUE) # p-spline of order 1
#' f_Xunif(x, model = "pspline2", K = K, plot_check = TRUE) # p-spline of order 2

f_Xunif <- function(
    x, model, K,
    adj_mat, m = NULL, M = NULL,
    fixed = TRUE, scale_Q = T,
    plot_check = F, n_sim = 1000) {
  # Discrete ------
  if (model == "iid") {
    result <- iid_standard(x = x, fixed = fixed)
  }
  if (model == "rw1") {
    result <- rw_standard(x = x, order = 1)
  }
  if (model == "rw2") {
    result <- rw_standard(x = x, order = 2)
  }
  if (model == "besag") {
    result <- besag_standard(x = x, adj_mat = adj_mat)
  }
  # Continuous -----
  if (model == "linear") {
    result <- linear_standard(x = x, m = m, M = M)
  }
  if (model == "pspline1") {
    result <- pspline_standard(x = x, K = K, order = 1, m = m, M = M)
  }
  if (model == "pspline2") {
    result <- pspline_standard(x = x, K = K, order = 2, m = m, M = M)
  }
  if (model == "pspline_2D") {
    result <- Pspline_2D_standard(x = x, K = K, m = m, M = M)
  }

  # Precision matrix unscaling -----
  if (!scale_Q) {
    result$precision <- result$precision / result$scaling_constant
  }

  # Check results via simulation-----------------------------

  if (plot_check) {
    realizations <- result$basis_distribution %*%
      t(mvtnorm::rmvnorm(n_sim, sigma = gen_inv(result$precision)))

    return(
      ggpubr::annotate_figure(
        check_GMRF(realizations,
          fixed = fixed,
          X_dist = result$X_distribution
        ),
        top = paste("Model:", model)
      )
    )
  }

  return(result)
}
