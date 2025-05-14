#' Standardized cluster effect
#'
#' `iid_standard()` provides a list of elements to build a standardized cluster effect.
#'
#' @param x An ordered factor, with more than 1 level.
#' @param fixed Logical, indicating whether the effect must be treated as fixed (TRUE) or random (FALSE). If it is treated as fixed, a 0 mean constraint is imposed By default, TRUE.
#'
#' @return A list of 6 elements, containing:
#' * `Q`: precision matrix
#' * `D`: basis matrix evaluated at `x`
#' * `C`: a positive number, representing the appropriate scaling constant
#' * `null_space`: a matrix, representing the null space of the precision matrix
#' * `X_distribution`: a numeric vector (or matrix for `pspline2D`) of values sampled from the Uniform distribution of X
#' * `basis_distribution`: basis matrix evaluated at `X_distribution`
#'
#'
#' @examples
#' K <- 20
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K))
#' #Example as fixed effect ------------
#' iid_standard(x)
#'
#' #Example as random effect ------------
#' iid_standard(x, fixed = FALSE)
iid_standard <- function(x, fixed = TRUE) {
  isOrderedFactor(x)

  K <- length(levels(x))
  X_dist <- 1:K
  D_dist <- diag(K)
  D <- as.matrix(
    fastDummies::dummy_cols(x,
      remove_selected_columns = T, ignore_na = T
    )
  )

  stand_iid <- standardize_GMRF(
    Q = diag(K), D = NULL,
    rank_def = 0, fixed = fixed
  )

  Q <- stand_iid$Q
  if (!is.null(stand_iid$A)) {
    S <- t(stand_iid$A)
  } else {
    S <- stand_iid$A
  }
  C <- stand_iid$C

  return(list(
    "Q" = Q,
    "D" = D,
    "C" = C,
    "null_space" = S,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
