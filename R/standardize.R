#' Standardize a GMRF
#'
#' `scale_GMRF()` returns the appropriate scaling constant for a given Gaussian effect, defined through a basis matrix and a set of Normally distributed coefficients with null mean and given precision matrix.
#'
#' @param Q A square matrix, representing the precision matrix on the coefficients.
#' @param D A matrix, representing the basis matrix of the effect. If not provided, it is set to an identity matrix of dimension equal to `Q`.
#' @param rank_def The rank deficiency of the matrix `Q`. If not provided, it is estimated.
#' @param fixed Logical, indicating whether the effect must be treated as fixed or random. By default, `TRUE`.
#' @param scale_Q Logical, indicating whether the precision matrix to be returned must already be scaled or not. By default, `TRUE`.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed or not. By default, `FALSE`.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot. By default, `n_sim=1000`.
#'
#' @return A list of 4 elements containing:
#' * `Q`: a matrix, representing the precision subject to the 0-mean constraint.
#' * `D`: a matrix, representing the basis subject to the 0-mean constraint.
#' * `C`: a positive number, representing the scaling constant.
#' * `A`: a matrix with 1 column, representing the 0-mean constraint as the coefficients of a linear combination. `NULL` for random effects for which the 0-mean constraint is not required, i.e. `fixed=FALSE`.
#'
#'
#' @examples
#' # Random walk example --------------------------------------
#' prec <- as.matrix(spam::precmat.RW1(10))
#' standardize_GMRF(prec, rank_def = 1) # fixed effect
#' standardize_GMRF(prec, rank_def = 1, fixed = FALSE) # random effect
#'
#' # Cluster effect example -----------------------------------
#' standardize_GMRF(diag(10)) # fixed effect
#' standardize_GMRF(diag(10), fixed = FALSE) # random effect

standardize_GMRF <- function(
    Q, D = NULL, rank_def = NULL,
    fixed = TRUE, scale_Q = TRUE,
    plot_check = FALSE, n_sim = 1000) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(1,nrow(Q))
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

  } else {
    C <- scale_GMRF(Q, D, rank_def = rank_def)

    result <- list(
      "precision" = Q * C,
      "basis" = D,
      "scaling_constant" = C,
      "mean_constraint" = NULL
    )

  }

  if (!scale_Q) {
    result$precision <- result$precision/result$scaling_constant
  }


  if (plot_check) {
    realizations <- result$basis %*%
      t(mvtnorm::rmvnorm(n_sim,
                         sigma = gen_inv(result$precision,
                                         rank_def = rank_def)
      ))

    return(check_GMRF(realizations,fixed = fixed,X_dist = NULL))
  }

  return(result)
}
