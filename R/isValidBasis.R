#' Check if a matrix is an appropriate basis
#'
#' `isValidBasis()` throws an error message if the first argument is not a non-empty matrix of the correct dimension. The error message explains which property is not respected.
#'
#' @param D A matrix.
#' @param Q A square matrix.
#'
#' @return NULL
#'
#' @details This function is used in other functions of the `scaleGMRF` package to check whether a matrix is a valid basis matrix for a given precision matrix.
#'
#' @examples
#' # Example with basis of correct dimension ------------------------
#' x <- seq(0, 1, length.out = 100)
#' isValidBasis(D = bspline(x, K = 10), Q = diag(10)) # NULL
#'
#' # Example with basis of incorrect dimension -----------------------
#' try(isValidBasis(D = bspline(x, K = 20), Q = diag(10))) # error
isValidBasis <- function(D, Q) {
  isSPSD(Q)

  if (!is.matrix(D)) {
    stop("`D` must be a matrix object. Instead it is a ", class(D), " object.")
  }

  if (nrow(D) == 0 | ncol(D) == 0) {
    stop("`D` must be a non-empty matrix.  Instead, it has dimension ", nrow(D), "x", ncol(D), ".")
  }

  if (anyNA(D)) {
    stop("`D` must not contain NA values.")
  }

  if (!ncol(D) == ncol(Q)) {
    stop("The number of columns of `Q` and of `D` must be equal. Instead ncol(Q)=", ncol(Q), "!=ncol(D)=", ncol(D), ".")
  }
}
