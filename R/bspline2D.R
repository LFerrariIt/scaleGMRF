#' Create a basis of cubic uniform 2D B-Spline functions
#'
#' `bspline2D()` creates a basis matrix with `K1`x`K2` columns, each of them evaluating the numeric vector `x` at a 2-dimensional cubic B-Spline function, defined through the Kronecker product of 2 unidimensional B-Spline functions.
#'
#' @param x1 A numeric vector, with values contained between 0 and 1.
#' @param x2 A numeric vector, with values contained between 0 and 1, of the same length of `x1`.
#' @param K1 A positive integer larger than 3, specifying the number of basis functions on the first dimension.
#'
#' @param K2 A positive integer larger than 3, specifying the number of basis functions on the second dimension.
#'
#' @return A matrix with `K1`x`K2` columns and a number of rows equal to the length of `x1` and `x2`.
#'
#' @examples
#' x_grid <- expand.grid(
#'   seq(0, 1, length.out = 100),
#'   seq(0, 1, length.out = 100)
#' )
#'
#' bspline2D(x_grid[, 1], x_grid[, 2], K1 = 5, K2 = 5)
#'
bspline2D <- function(x1, x2, K1, K2) {
  if (K1 < 4) {
    stop(paste0("`K1` must be larger than 3."))
  }

  if (K2 < 4) {
    stop(paste0("`K2` must be larger than 3."))
  }

  if (!is.numeric(x1)) {
    stop(paste0("`x1` must be a numeric vector."))
  }

  if (!is.numeric(x2)) {
    stop(paste0("`x2` must be a numeric vector."))
  }

  if (length(x1) != length(x2)) {
    stop(paste0("The vectors `x1` and `x2` must have the same length. Instead,", length(x1), "!=", length(x2), "."))
  }

  B_1 <- bspline(x1, K1)
  B_2 <- bspline(x2, K2)

  B <- matrix(NA, nrow = length(x1), ncol = K1 * K2)

  for (i in 1:length(x1)) {
    B[i, ] <- kronecker(B_2[i, ], B_1[i, ])
  }
  return(B)
}
