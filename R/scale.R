#' Expectation-based scaling for an LGM effect
#'
#' `scale_GMRF()` returns the appropriate scaling constant for a given effect of an LGM model, defined through a basis matrix and a precision matrix.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
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
