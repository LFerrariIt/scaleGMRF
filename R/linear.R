#' Standardized linear effect
#'
#' @description
#' `linear_standard()` provides a list of elements to build a standardized linear effect.
#'
#' @param x A numeric vector.
#' @param m A number indicating the lower boundary of the support of X. If not provided, it is set to the minimum value from `x`.
#' @param M A number indicating the upper boundary of the support of X. If not provided, it is set to the maximum value from `x`.
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
#' x <- seq(0, 1, length.out = 100)
#'
#' linear_standard(x)
#'
linear_standard <- function(x, m = NULL, M = NULL) {
  # Error messages on the arguments

  m_M <- checkBoundaries(x = x, m = m, M = M)

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
