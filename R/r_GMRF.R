#' Generate realizations of a GMRF
#'
#' `r_GMRF()` returns a matrix where each column is a realization of the specified Gaussian effect.
#'
#' @inheritParams mean0_GMRF
#' @param sigma2 Positive number, indicating the value of the variance parameter. By default, 1.
#' @param n_sim Positive integer, indicating the number of realizations to generate. By default, `n_sim=100`.
#'
#' @return A numeric matrix, containing a realization of the process in each column.
#'
#' @details This function is used in the `standardize_GMRF()` and `standardize_X_unif()` functions of the `scaleGMRF` package to generate realizations of Gaussian effects.
#'
#' @examples
#' # Example for an i.i.d. effect --------------------
#' r_GMRF(Q = diag(10))
#'
#' # Example for a linear effect ------------------------
#' r_GMRF(Q = matrix(1), D = matrix(runif(100)))
#'
#' # Example for a random walk effect of order 1 --------
#' r_GMRF(Q = as.matrix(spam::precmat.RW1(10)), rank_def = 1)
#'
r_GMRF <- function(Q, D = NULL, sigma2 = 1,
                   rank_def = NULL, n_sim = 100) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  realizations <- D %*%
    t(mvtnorm::rmvnorm(n_sim, sigma = sigma2 * gen_inv(Q, rank_def = rank_def)))

  return(realizations)
}
