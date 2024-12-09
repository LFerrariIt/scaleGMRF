#' Test if a Matrix is Symmetric Positive (Semi-)Definite
#'
#' `isSPSD()` throws an error message if the argument is not a non-empty symmetric positive (semi-)definite matrix without NA entries and not nilpotent. The error message explains which property is not respected.
#'
#' @param Q A matrix object.
#'
#' @return NULL
#'
#'
#' @examples
#' isSPSD(diag(10)) # NULL
#' try(isSPSD(diag(0))) # error
#' try(isSPSD(diag(NA, ncol = 3, nrow = 3))) # error
#' try(isSPSD(matrix(1, ncol = 3, nrow = 2))) # error

isSPSD <- function(Q) {
  if (!is.matrix(Q)) {
    stop("`Q` must be a matrix, not a ", class(Q), ".")
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

  if (!isSymmetric(Q)) {
    stop("`Q` must be symmetric.")
  }

  V <- eigen(Q, symmetric = T, only.values = T)$values

  if (any(V < -0.001) | all(V == 0)) {
    stop("`Q` must be a positive definite or positive semi-definite matrix, but not nilpotent.")
  }
}
