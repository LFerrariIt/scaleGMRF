#' Standardized random walk effects
#'
#' @description
#' `rw1_standard()` provides a list of elements to build a standardized random walk effect of order 1.
#'
#' `rw2_standard()` provides a list of elements to build a standardized random walk effect of order 2.
#'
#' @param x An ordered factor, with more than 1 level.
#'
#' @return A list of 6 elements, containing:
#' \itemize{
#' \item{`precision`: precision matrix.}
#' \item{`basis`:  basis matrix evaluated at `x`}
#' \item{`scaling_constant`: a positive number, representing the appropriate scaling constant.}
#'  \item{`null_space`: a matrix, representing the null space of the precision matrix.}
#'  \item{`X_distribution`: a numeric vector sampled from the Uniform distribution on X.}
#'  \item{`basis_distribution`: basis matrix evaluated at `X_distribution`.}
#' }
#'
#' @examples
#' K <- 20
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K))
#'
#' rw1_standard(x)
#' rw1_standard(x)
#' rw2_standard(x)
#' rw2_standard(x)
#'
rw1_standard <- function(x) {
  isOrderedFactor(x)

  K <- length(levels(x))
  X_dist <- 1:K
  D_dist <- diag(K)
  D <- as.matrix(
    fastDummies::dummy_cols(x,
      remove_selected_columns = T, ignore_na = T
    )
  )
  Q <- as.matrix(spam::precmat.RW1(K))
  S <- matrix(1, ncol = 1, nrow = K)
  C <- scale_GMRF(Q = Q, D = NULL, rank_def = 1)

  return(list(
    "precision" = Q * C,
    "basis" = D,
    "scaling_constant" = C,
    "null_space" = S,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
