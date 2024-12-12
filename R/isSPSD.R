#' Check if a Matrix is Symmetric Positive (Semi-)Definite
#'
#' `isSPSD()` throws an error message if the argument is not a non-empty Symmetric Positive (Semi-)Definite (SPSD) matrix, without missing values, and not nilpotent. The error message explains which property is not respected.
#'
#' @param Q A matrix object.
#'
#' @return NULL
#' @details This function is used in other functions of the `scaleGMRF` package to check whether a matrix is a valid precision matrix.
#'
#' @examples
#' # Example with SPSD matrix ---------------------------------
#' isSPSD(diag(10)) # NULL
#'
#' # Example with empty matrix ---------------------------------
#' try(isSPSD(diag(0))) # error
#'
#' # Example with missing values ---------------------------------
#' try(isSPSD(diag(NA, ncol = 3, nrow = 3))) # error
#'
#' # Example with non-square matrix ---------------------------------
#' try(isSPSD(matrix(1, ncol = 3, nrow = 2))) # error
#'
isSPSD <- function(Q) {
  if (!is.matrix(Q)) {
    stop("`Q` must be a matrix object, not a ", class(Q), " object.")
  }

  if (!ncol(Q) == nrow(Q)) {
    stop("`Q` must be a square matrix. Instead, it has dimension ", nrow(Q), "x", ncol(Q), ".")
  }

  if (ncol(Q) == 0) {
    stop("`Q` must be a non-empty matrix. Instead, it has dimension ", nrow(Q), "x", ncol(Q), ".")
  }

  if (anyNA(Q)) {
    stop("`Q` must not contain NA values.")
  }

  if (!isSymmetric(unname(Q))) {
    stop("`Q` must be symmetric.")
  }

  V <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values

  if (any(V < -0.001) | all(V == 0)) {
    stop("`Q` must be a positive definite or positive semi-definite matrix, but not nilpotent.")
  }
}
