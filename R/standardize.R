#' Standardization of a GMRF
#'
#' `standardize_GMRF()` standardizes a Gaussian effect and returns the elements needed to specify the standardized effect in INLA.
#'
#' @inheritParams mean0_GMRF
#' @param fixed Logical, indicating whether the effect must be treated as fixed or random. An effect is always considered fixed if it only has one basis function. By default, `FALSE`.
#' @param scale_Q Logical, indicating whether the precision matrix must be scaled or not. By default, `TRUE`.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed or not. By default, `FALSE`.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot. By default, `n_sim=100`.
#'
#' @return A list of 4 elements containing:
#' * `Q`: a matrix, representing the precision matrix.
#' * `D`: a matrix, representing the basis.
#' * `C`: a positive number, representing the scaling constant.
#' * `A`: a matrix with 1 column, representing the 0-mean constraint as the coefficients of a linear combination. `NULL` for `fixed=FALSE`.
#'
#'
#' @details This function is one of the two main functions of the `scaleGMRF` package. More details about this function can be found in `vignette("standardization",package="scaleGMRF")`.
#'
#' @examples
#'
#' # Example for an i.i.d effect ------------------------
#' prec <- diag(10)
#' # fixed effect
#' standardize_GMRF(prec)
#' # random effect
#' standardize_GMRF(prec, fixed = FALSE)
#'
#' # Example for a linear effect ------------------------
#' basis <- matrix(runif(100))
#' # fixed effect
#' standardize_GMRF(Q = matrix(1), D = basis)
#'
#' # Example for a random walk effect of order 1 ------------------------
#' prec <- as.matrix(spam::precmat.RW1(10))
#' # fixed effect
#' standardize_GMRF(prec, rank_def = 1)
#' # random effect
#' standardize_GMRF(prec, rank_def = 1, fixed = FALSE)
#'
standardize_GMRF <- function(
    Q, D = NULL, rank_def = NULL,
    fixed = FALSE, scale_Q = TRUE,
    plot_check = FALSE, n_sim = 100) {
  isSPSD(Q)

  if (is.null(D)) {
    D <- diag(1, nrow(Q))
  }

  isValidBasis(D, Q)

  if (!is.logical(fixed)) {
    stop("fixed must be logical.")
  }

  if (!is.logical(scale_Q)) {
    stop("fixed must be logical.")
  }

  if (!is.logical(plot_check)) {
    stop("`plot_check` must be logical.")
  }

  if (!is.numeric(n_sim)) {
    stop("`n_sim` must be a positive number.")
  }

  if (fixed | ncol(D)==1) {
    constr_results <- mean0_GMRF(Q, D, rank_def)

    C <- unname(scale_GMRF(constr_results$Q, constr_results$D))

    result <- list(
      "precision" = constr_results$Q * C,
      "basis" = constr_results$D,
      "scaling_constant" = C,
      "mean_constraint" = constr_results$A
    )
  } else {
    C <- unname(scale_GMRF(Q, D, rank_def = rank_def))

    result <- list(
      "precision" = Q * C,
      "basis" = D,
      "scaling_constant" = C,
      "mean_constraint" = NULL
    )
  }

  if (!scale_Q) {
    result$precision <- result$precision / result$scaling_constant
  }


  if (plot_check) {
    realizations <- r_GMRF(
      Q = result$precision,
      D = result$basis, n_sim = n_sim
    )

    print(check_GMRF(realizations, fixed = fixed, X_dist = NULL))
  }

  return(result)
}
