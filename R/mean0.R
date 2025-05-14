#' 0-mean constraint imposition on GMRF.
#'
#' `mean0_GMRF` returns the precision and basis matrices of a Gaussian effect subject to a 0-mean constraint, along with the constraint itself.
#'
#' @param Q A matrix, representing the precision matrix on the coefficients of the Gaussian effect.
#' @param D A matrix, representing the basis of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def A non-negative integer, indicating the rank deficiency of the matrix `Q`. If not provided, it is estimated.
#'
#' @return A list of 3 elements containing:
#' * `Q`: a matrix, representing the precision matrix subject to the 0-mean constraint.
#' * `D`: a matrix, representing the basis subject to the 0-mean constraint.
#' * `A`: a matrix with 1 row, representing the 0-mean constraint as the coefficients of a linear combination.
#'
#' @details This function is used in the `standardize_GMRF()` function of the `scaleGMRF` package to impose the 0-mean constraint, which is required for the standardization of fixed effects. The function first checks whether the effect already respects the 0-mean constraint. If not, the function proceeds by appropriately modifying either the basis (for effects with a single basis function) or the precision matrix (for effects with more than 1 basis function) following Equation 2.29 from Rue & Held, 2005. The correct imposition of the 0-mean constraint can be checked on the output of the function computing `mean(D)` if there is only one basis function or `Q%*%A` for effects with more than one basis function: the result should be a 0 vector of length equal to the number of basis functions.
#'
#' # Example for an i.i.d. effect -------------------------------
#' result <- mean0_GMRF(diag(10))
#' result$Q %*% t(result$A)
#'
#' # Example for a linear effect -------------------------------
#' result <- mean0_GMRF(Q=matrix(1),D=matrix(runif(100)))
#' mean(result$D)
#'
#' # Example for a random walk effect of order 1------------------------
#' result <- mean0_GMRF(as.matrix(spam::precmat.RW1(10)), rank_def = 1)
#' result$Q %*% t(result$A)
#'

mean0_GMRF <- function(Q, D = NULL, rank_def = NULL) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  A <- matrix(colMeans(D), nrow = 1)

  if (all(round(Q %*% t(A),10) == 0)) {
    return(list(
      "Q" = Q,
      "D" = D,
      "A" = A
    ))
  } else {
    if (ncol(A) == 1) {
      D_constr <- D - as.numeric(A)

      return(list(
        "Q" = Q,
        "D" = D_constr,
        "A" = A
      ))
    } else {
      Sigma_constr <- Sigma - (Sigma %*% t(A) %*% A %*% Sigma) / as.numeric(A %*% Sigma %*% t(A))

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
