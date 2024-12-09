#' Create a basis of cubic uniform B-Spline functions
#'
#' `bspline()` creates a basis matrix with `K` columns, each of them evaluating the numeric vector `x` at a cubic B-Spline function defined on equidistant nodes on the `m,M` interval.
#'
#' @param x A numeric vector, with values contained between 0 and 1.
#' @param K A positive integer larger than 3, specifying the number of basis functions.
#'
#' @param m,M A pair of numbers indicating respectively the lower and upper boundaries of the support. By default, 0 and 1.
#'
#' @return A matrix with `K` columns and a number of rows equal to the length of `x`.
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' bspline(x, K = 10)
bspline <- function(x, K, m = 0, M = 1) {
  if (K < 4) {
    stop("`K` must be larger than 3.")
  }

  if (!is.numeric(x)) {
    stop("`x` must be a numeric vector.")
  }

  if (any(x < m) | any(x > M) | any(is.na(x))) {
    stop("The vector `x` must be normalized to lie between m and M and not contain NA values.")
  }

  x_norm <- (x - m) / (M - m)

  dx <- 1 / (K - 3)
  knots <- seq(-3 * dx, 1 + 3 * dx, by = dx)
  B <- splines::spline.des(knots, x_norm, ord = 4, 0 * x, outer.ok = T)$design

  return(B)
}
