#' Scaling constant for a GMRF
#'
#' `scale_GMRF()` returns the appropriate scaling constant for a Gaussian effect.
#'
#' @inheritParams mean0_GMRF
#'
#' @return `scaling_constant`A positive number.
#'
#' @details This function is used in the `standardize_GMRF()` function of the `scaleGMRF` package to appropriately scale Gaussian effects multiplying their precision matrix by the scaling constant returned by the `scale_GMRF()` function.
#'
#' @examples
#' # Example for an i.i.d. effect --------------------
#' scale_GMRF(Q = diag(10))
#'
#' # Example for a linear effect ------------------------
#' scale_GMRF(Q = matrix(1), D = matrix(runif(100)))
#'
#' # Example for a random walk effect of order 1 --------
#' scale_GMRF(Q = as.matrix(spam::precmat.RW1(10)), rank_def = 1)
#'
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
  names(C) <- "scaling_constant"

  return(C)
}
