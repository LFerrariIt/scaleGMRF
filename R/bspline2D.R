#' Create a basis of cubic uniform 2D B-Spline functions
#'
#' `bspline2D()` creates a basis matrix with `K1`x`K2` columns, each of them evaluating the numeric vector `x` at a 2-dimensional cubic B-Spline function, defined through the Kronecker product of 2 unidimensional B-Spline functions.
#'
#' @param x A matrix of 2 columns, containing values between 0 and 1.
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
#' bspline2D(x_grid, K1 = 5, K2 = 5)
#'
bspline2D <- function(x, K1, K2) {

  if (K1 < 4) {
    stop(paste0("`K1` must be larger than 3."))
  }

  if (K2 < 4) {
    stop(paste0("`K2` must be larger than 3."))
  }

  if (!is.matrix(x)) {
    stop(paste0("`x` must be a matrix."))
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop(paste0("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), "."))
  }

  if (any(x < 0) | any(x > 1) | any(is.na(x))) {
    stop(paste0("The entries of `x` must be normalized to lie between 0 and 1 and must not contain NA values."))
  }

  x1 <- x[,1]
  x2 <- x[,2]

  B_1 <- bspline(x1, K1)
  B_2 <- bspline(x2, K2)

  B <- matrix(NA, nrow = length(x1), ncol = K1 * K2)

  for (i in 1:length(x1)) {
    B[i, ] <- kronecker(B_2[i, ], B_1[i, ])
  }
  return(B)
}
