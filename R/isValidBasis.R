#' Test if a matrix is an appropriate basis
#'
#' `isValidBasis()` throws an error message if the first argument is not a non-empty matrix of the correct dimension. The error message explains which property is not respected.
#'
#' @param D a matrix.
#' @param Q a square matrix.
#'
#' @return NULL
#' @export
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#' isValidBasis(D = bspline(x, K = 10), Q = diag(10)) # NULL
#'
#' try(isValidBasis(D = bspline(x, K = 20), Q = diag(10))) # error
isValidBasis <- function(D, Q) {
  tryCatch(
    {
      isSPSD(Q)
    },
    error = function(e) {
      print("`Q` must be a valid precision matrix:")
      print(e)
    }
  )

  if (!is.matrix(D)) {
    stop("`D` must be a matrix. Instead it is a", class(D), ".")
  }

  if (nrow(D) == 0 | ncol(D) == 0) {
    stop("`D` must be a non-empty matrix.")
  }

  if (anyNA(D)) {
    stop("`D` must not contain NA values.")
  }

  if (!ncol(D) == ncol(Q)) {
    stop("The number of columns of `Q` and of `D` must be equal. Instead ncol(Q)=,", ncol(Q), "!=ncol(D)=", ncol(D), ".")
  }
}
