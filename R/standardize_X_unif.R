#' Standardized effects for a Uniform covariate
#'
#' `standardize_X_unif()` provides a list of elements to build various standardized effects (chosen through the `model` argument) under the assumption that the corresponding covaiate follows a Uniform distribution.

#' @param x A numeric vector (for `model="linear","pspline1","pspline2","pspline_2D"`) or an ordered factor with more than 1 level (for `model="iid","rw1","rw2","besag"``).
#' @param model String to select the type of effect between `"iid","rw1","rw2","besag","linear","pspline1","pspline2","pspline_2D"`.
#' @param K A positive integer larger than 3, specifying the number of basis functions. Required for `model="pspline1","pspline2","pspline_2D"`.
#' @param adj_mat Matrix required for `model="besag"`, representing the adjacency matrix from the lattice. See Details.
#' @param m,M A pair of numbers indicating respectively the lower and upper boundaries of the support of X. If not provided, they are set to the extreme values from `x`. Required for `model="linear","pspline1","pspline2","pspline_2D"`.
#' @param fixed Logical, indicating whether the effect must be treated as fixed or random. By default, `FALSE`.
#' @param scale_Q Logical, indicating whether the precision matrix to be returned must already be scaled or not. By default, `TRUE`.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed or not. By default, `FALSE`.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot. By default, `n_sim=100`.
#'
#' @details This function is designed so that the first three elements of the output list can directly be used to specify a model effect in INLA through the model `generic0` (see `inla.doc(generic0)`). More details about this function can be found in `vignette("standardize_X_unif",package="scaleGMRF")`.
#'
#' @inherit iid_standard return
#'
#' @examples
#' # Examples for discrete covariates --------------------------
#' K <- 20 # number of levels
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K))
#' standardize_X_unif(x, model = "iid", plot_check = TRUE) # fixed iid effect
#' standardize_X_unif(x, model = "rw1", plot_check = TRUE) # random walk of order 1
#'
#' # Examples for continuous covariates --------------------------
#' x <- runif(1000, min = 2, max = 5)
#' standardize_X_unif(x, model = "linear", plot_check = TRUE) # linear effect
#' K <- 20 # number of b-spline basis functions
#' standardize_X_unif(x, model = "pspline2", K = K, plot_check = TRUE) # p-spline of order 2
standardize_X_unif <- function(
    x, model, K,
    adj_mat, m = NULL, M = NULL,
    fixed = FALSE, scale_Q = T,
    plot_check = F, n_sim = 100) {
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
    result <- pspline_2D_standard(x = x, K = K, m = m, M = M)
  }

  # Precision matrix unscaling -----
  if (!scale_Q) {
    result$Q <- result$Q / result$C
  }

  # Check results via simulation-----------------------------

  if (plot_check) {
    realizations <- r_GMRF(
      Q = result$Q,
      D = result$basis_distribution,
      n_sim = n_sim
    )

    print(
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
