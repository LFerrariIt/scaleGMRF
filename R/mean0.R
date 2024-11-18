#' Precision matrix of a GMRF subject to 0-mean constraint.
#'
#' `mean0_GMRF` returns the precision matrix subject to a 0-mean constraint, based on the provided precision and basis matrices.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
#'
#' @return A list of 3 elements containing:
#' \itemize{
#' \item{`Q`: a matrix, representing the precision matrix under the 0-mean constraint.}
#' \item{`D`: a matrix, representing the basis under the 0-mean constraint.}
#' #' \item{`A`: a matrix, representing the 0-mean constraint as the coefficients of a linear combination.}
#' }
#' @examples
#' mean0_GMRF(diag(10))
#' mean0_GMRF(as.matrix(spam::precmat.RW1(10)), rank_def = 1)
mean0_GMRF <- function(Q, D = NULL, rank_def = NULL) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  A <- matrix(colMeans(D), ncol = 1)

  if (all(Q %*% A == 0)) {
    return(list(
      "Q" = Q,
      "D" = D,
      "A" = A
    ))
  } else {
    if (nrow(A) == 1) {
      D_constr <- D - as.numeric(A)

      return(list(
        "Q" = Q,
        "D" = D_constr,
        "A" = A
      ))
    } else {
      Sigma_constr <- Sigma - (Sigma %*% A %*% t(A) %*% Sigma) / as.numeric(t(A) %*% Sigma %*% A)

      Q_constr <- gen_inv(Sigma_constr, rank_def = NULL)

      isSPSD(Q_constr)

      return(list(
        "Q" = Q_constr,
        "D" = D,
        "A" = A
      ))
    }
  }
}
