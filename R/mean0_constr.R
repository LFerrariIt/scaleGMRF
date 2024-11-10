#' Precision matrix of a GMRF subject to 0-mean constraint.
#'
#' `prec_under_0mean` returns the precision matrix subject to a 0-mean constraint, based on the provided precision and basis matrices.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
#'
#' @return A list of 2 elements containing:
#' \itemize{
#' \item{`Q_constr`: a matrix, representing the precision matrix under the 0-mean constraint.}
#' \item{`constr`: a numeric vector, representing the 0-mean constraint.}
#'}
#' @examples
#' prec_under_0mean(diag(10))
#' prec_under_0mean(as.matrix(spam::precmat.RW1(10)), rank_def = 1)
prec_under_0mean <- function(Q, D = NULL, rank_def = NULL) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  A <- matrix(colMeans(D), nrow = 1)

  if (all(Q %*% t(A) == 0)) {
    return(list(
      "Q_constr" = Q,
      "constr" = A
    ))
  } else {
    Sigma_constr <- Sigma - (Sigma %*% t(A) %*% A %*% Sigma) / as.numeric(A %*% Sigma %*% t(A))

    Q_constr <- gen_inv(Sigma_constr, rank_def = NULL)

    isSPSD(Q_constr)

    return(list(
      "Q_constr" = Q_constr,
      "constr" = A
    ))
  }
}
