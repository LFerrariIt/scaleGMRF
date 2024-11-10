#' Standardize a GMRF
#'
#' `ex_scale()` returns the appropriate scaling constant for a given effect of an LGM model, defined through a basis matrix and a precision matrix.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
#' @param fixed logical indicating whether the effect is to be considered fixed.
#'
#' @return A square matrix representing the precision matrix after appropriate scaling.
#'
#' @examples
#' prec <- as.matrix(spam::precmat.RW1(10))
#' standardize_GMRF(prec, rank_def=1)
#' standardize_GMRF(prec, rank_def=1, fixed=T)
#' standardize_GMRF(diag(10))
#' standardize_GMRF(diag(10), fixed=T)
standardize_GMRF <- function(Q, D = NULL, rank_def = NULL, fixed = FALSE) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  if (fixed) {
    Q_constr <- prec_under_0mean(Q, D, rank_def)$Q_constr

    C <- ex_scale(Q_constr, D, rank_def = NULL)

    Q_standardized <- Q_constr * C

  } else {
    C <- ex_scale(Q, D, rank_def = rank_def)

    Q_standardized <- Q * C

  }

    return(Q_standardized)
}
