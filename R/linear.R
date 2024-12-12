#' Standardized linear effects
#'
#' `linear_standard()` provides a list of elements to build a standardized linear effect.
#'
#' @param x A numeric vector, indicating the values on which the basis must be evaluated.
#' @param m,M A pair of numbers indicating respectively the lower and upper boundary of the support of X. If not provided, they are set to the minimum and maximum values of `x`.
#'
#' @inherit iid_standard return
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#'
#' linear_standard(x)
#'
linear_standard <- function(x, m = NULL, M = NULL) {
  # Error messages on the arguments

  m_M <- xBoundaries(x = x, m = m, M = M)

  m <- m_M[1]

  M <- m_M[2]

  x_mean <- (m + M) / 2
  x_var <- (M - m)^2 / 12

  Q <- diag(1)
  D <- matrix((x - x_mean), ncol = 1, nrow = length(x))
  C <- x_var

  X_dist <- seq(m, M, length.out = 1000)
  D_dist <- matrix((X_dist - x_mean), ncol = 1, nrow = length(X_dist))

  return(list(
    "precision" = Q * C,
    "basis" = D,
    "scaling_constant" = C,
    "null_space" = NULL,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
