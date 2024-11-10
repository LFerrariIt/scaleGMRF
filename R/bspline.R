#' Create a basis of cubic uniform B-Spline functions
#'
#' `bspline()` creates a basis matrix with `K` columns, each of them evaluating the numeric vector `x` at a cubic B-Spline function defined on equidistant nodes on the [0,1] interval.
#'
#' @param x A numeric vector, with values contained between 0 and 1.
#' @param K A positive integer larger than 3, specifying the number of basis functions.
#'
#' @return A matrix with `K` columns and a number of rows equal to the length of `x`.
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' bspline(x, K = 10)

bspline <- function(x, K) {
  if (K < 4) {
    stop(paste0("`K` must be larger than 3."))
  }

  if (!is.numeric(x)) {
    stop(paste0("`x` must be a numeric vector."))
  }

  if (any(x < 0) | any(x > 1) | any(is.na(x))) {
    stop(paste0("The vector `x` must be normalized to lie between 0 and 1 and not contain NA values."))
  }

  dx <- 1 / (K - 3)
  knots <- seq(-3 * dx, 1 + 3 * dx, by = dx)
  B <- splines::spline.des(knots, x, ord = 4, 0 * x, outer.ok = T)$design

  return(B)
}
