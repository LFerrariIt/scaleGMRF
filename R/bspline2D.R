#' Create a basis of cubic uniform 2D B-Spline functions
#'
#' `bspline_2D()` creates a basis matrix with `K1`x`K2` columns, each of them evaluating the numeric vector `x` at a 2-dimensional cubic B-Spline function, defined through the Kronecker product of 2 unidimensional B-Spline functions on the rectangular support delimited by `m,M`.
#'
#' @param x A matrix of 2 columns. The first column must contain values between `m[1]` and `M[1]`. The second column must contain values between `m[2]` and `M[2]`.
#' @param K Two positive integers larger than 3, specifying the number of basis functions on each dimension.

#' @param m,M Two pairs of numbers indicating respectively the lower and upper bounds of the rectangular support on which to define the basis functions. By default, `c(0,0)` and `c(1,1)`.
#'
#' @return A matrix with `K1`x`K2` columns and a number of rows equal to the length of `x1` and `x2`.
#'
#' @examples
#' x_grid <- as.matrix(expand.grid(
#'   seq(0, 1, length.out = 100),
#'   seq(0, 1, length.out = 100)
#' ))
#'
#' bspline_2D(x_grid, K=c(5,5))
#'
bspline_2D <- function(x, K, m = c(0, 0), M = c(1, 1)) {
  K1 <- K[1]
  K2 <- K[2]

  if (K[1] < 4) {
    stop("`K[1]` must be larger than 3.")
  }

  if (K[2] < 4) {
    stop("`K[2]` must be larger than 3.")
  }

  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), ".")
  }

  x1 <- x[, 1]
  x2 <- x[, 2]

  B_1 <- bspline(x1, K[1], m = m[1], M = M[1])
  B_2 <- bspline(x2, K[2], m = m[2], M = M[2])

  B <- matrix(NA, nrow = length(x1), ncol = K[1] * K[1])

  for (i in 1:length(x1)) {
    B[i, ] <- kronecker(B_2[i, ], B_1[i, ])
  }
  return(B)
}
