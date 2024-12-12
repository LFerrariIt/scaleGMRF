#' Create a basis of equidistant cubic 2D B-Spline functions
#'
#' `bspline_2D()` creates a basis matrix with `K[1]`x`K[2]` columns, each of them evaluating the values in `x` at a bivariate cubic B-Spline function, defined as the Kronecker product between 2 univariate B-Spline functions.
#'
#' @param x A numeric matrix with 2 columns, indicating the values on which the basis must be evaluated.
#' @param K Two positive integers larger than 3, indicating the number of basis functions for each dimension.
#' @param m,M Two pairs of numbers, indicating respectively the lower and upper boundaries of the rectangular support on which to define the basis functions. If not provided, they are estimated using the minimum and maximum values from the columns of `x`.
#'
#' @return A matrix with `K[1]`x`K[2]` columns and a number of rows equal to the one of `x`.
#'
#' @details This function is used in the `pspline_2D_standard()` function of the `scaleGMRF` package to build a basis of equidistant cubic 2D B-Splines on a rectangular support.
#'
#' @examples
#' x_grid <- as.matrix(expand.grid(
#'   seq(0, 1, length.out = 100),
#'   seq(0, 1, length.out = 100)
#' ))
#'
#' bspline_2D(x_grid, K = c(5, 5))
#'
bspline_2D <- function(x, K, m = NULL, M = NULL) {
  if (!is.numeric(K) | !is.vector(K) | length(K) != 2) {
    stop("`K` must be a vector containing 2 numbers. Instead, K=", K, ".")
  }

  if (K[1] < 4) {
    stop("`K[1]` must be larger than 3. Instead, K[1]=", K[1], ".")
  }

  if (K[2] < 4) {
    stop("`K[2]` must be larger than 3.Instead, K[2]=", K[2], ".")
  }

  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), ".")
  }

  if (is.null(m)) {
    m <- c(min(x[, 1]), min(x[, 2]))
  }
  if (is.numeric(m) & length(m) == 2) {
    m1 <- m[1]
    m2 <- m[2]
  } else {
    stop(
      "`m` must contain 2 numeric elements or be NULL. Instead, m=",
      paste0(m, collapse = ","), "."
    )
  }

  if (is.null(M)) {
    M <- c(max(x[, 1]), max(x[, 2]))
  }
  if (is.numeric(M) & length(M) == 2) {
    M1 <- M[1]
    M2 <- M[2]
  } else {
    stop(
      "`M` must contain 2 numeric elements. Instead, M=",
      paste0(M, collapse = ","), "."
    )
  }

  if (any(m >= M)) {
    stop(
      "`M` must be larger than `m`. Instead, m=",
      paste0(m, collapse = ","), " and M=",
      paste0(M, collapse = ","), "."
    )
  }

  B_1 <- bspline(x[, 1], K[1], m = m[1], M = M[1])
  B_2 <- bspline(x[, 2], K[2], m = m[2], M = M[2])

  B <- matrix(NA, nrow = nrow(x), ncol = K[1] * K[1])

  for (i in 1:nrow(x)) {
    B[i, ] <- kronecker(B_2[i, ], B_1[i, ])
  }
  return(B)
}
