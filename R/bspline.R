#' Create a basis of equidistant cubic B-Spline functions
#'
#' `bspline()` creates a basis matrix with `K` columns, each of them evaluating the numeric vector `x` at a cubic B-Spline function defined on equidistant nodes on the `m,M` interval.
#'
#' @inheritParams linear_standard
#' @param K A positive integer larger than 3, specifying the number of basis functions.
#'
#' @return A matrix with `K` columns and a number of rows equal to the length of `x`.
#'
#' @details This function is used in the `pspline_standard()` function of the `scaleGMRF` package to build a basis of equidistant cubic B-Splines on a given interval. It is also used in the function `bspline_2D()` to create bivariate B-Spline bases on a rectangular support.
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' bspline(x, K = 10)
bspline <- function(x, K, m = NULL, M = NULL) {
  m_M <- xBoundaries(x = x, m = m, M = M)
  m <- m_M[1]
  M <- m_M[2]
  x_norm <- (x - m) / (M - m)

  if (!is.numeric(K) | !is.vector(K) | length(K) != 1) {
    stop("`K` must be a number. Instead, it is a ", class(K), ".")
  }

  if (K < 4 | K %% 1 != 0) {
    stop("`K` must be a whole number larger than 3. Instead, K=", K, ".")
  }

  dx <- 1 / (K - 3)
  knots <- seq(-3 * dx, 1 + 3 * dx, by = dx)
  B <- splines::spline.des(knots, x_norm, ord = 4, 0 * x, outer.ok = T)$design

  return(B)
}
