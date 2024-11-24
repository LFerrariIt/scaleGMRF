#' Standardize a GMRF
#'
#' `scale_GMRF()` returns the appropriate scaling constant for a given effect of an LGM model, defined through a basis matrix and a precision matrix.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
#' @param fixed logical indicating whether the effect is to be considered fixed.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed (TRUE) or not (FALSE). By default, FALSE.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot.
#'
#' @return A list of 3 or 4 elements containing:
#' \itemize{
#' \item{`Q`: a matrix, representing the precision matrix under the 0-mean constraint.}
#' \item{`D`: a matrix, representing the basis under the 0-mean constraint.}
#' \item{`C`: a number, representing the scaling constant.}
#' \item{`A`: a matrix, representing the 0-mean constraint as the coefficients of a linear combination.}
#' }
#'
#' @examples
#' prec <- as.matrix(spam::precmat.RW1(10))
#' standardize_GMRF(prec, rank_def = 1)
#' standardize_GMRF(prec, rank_def = 1, fixed = FALSE)
#' standardize_GMRF(diag(10))
#' standardize_GMRF(diag(10), fixed = FALSE)
standardize_GMRF <- function(Q, D = NULL, rank_def = NULL, fixed = TRUE, plot_check = FALSE, n_sim = 1000) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  if (!is.logical(fixed)) {
    stop("fixed must be logical.")
  }

  if (fixed) {
    constr_results <- mean0_GMRF(Q, D, rank_def)
    Q <- constr_results$Q
    D <- constr_results$D
    A <- constr_results$A

    C <- scale_GMRF(Q, D)

    result <- list(
      "precision" = Q * C,
      "basis" = D,
      "scaling_constant" = C,
      "mean_constraint" = A
    )

    realizations <- result$basis %*%
      t(mvtnorm::rmvnorm(n_sim,
        sigma = gen_inv(
          result$precision,
          rank_def
        )
      ))
  } else {
    C <- scale_GMRF(Q, D, rank_def = rank_def)

    result <- list(
      "precision" = Q * C,
      "basis" = D,
      "scaling_constant" = C
    )

    realizations <- result$basis %*%
      t(mvtnorm::rmvnorm(n_sim,
        sigma = gen_inv(result$precision,
          rank_def = rank_def
        )
      ))
  }

  if (plot_check) {
    return(check_GMRF(realizations))
  }

  return(result)
}
