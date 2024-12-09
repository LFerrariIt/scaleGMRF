#' Expectation-based scaling for an LGM effect
#'
#' `scale_GMRF()` returns the appropriate scaling constant for a Gaussian effect, defined through a basis matrix and a set of Normally distributed coefficients with null mean and given precision matrix.
#'
#' @inheritParams mean0_GMRF
#'
#' @return A positive number to be used for scaling.
#'
#' @examples
#' scale_GMRF(Q = diag(10))
#' prec <- as.matrix(spam::precmat.RW1(10))
#' scale_GMRF(Q = prec, rank_def = 1)

scale_GMRF <- function(Q, D = NULL, rank_def = NULL) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  var_diag <- rowSums(D * (D %*% Sigma))

  if (any(var_diag < 0)) {
    stop("Invalid case: the conditional variance entries cannot be negative.")
  }

  if (all(var_diag == 0)) {
    stop("Invalid case: the conditional variance entries cannot be all null.")
  }

  C <- mean(var_diag)

  return(C)
}
