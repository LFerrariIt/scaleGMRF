#' Test if a Matrix is Symmetric Positive (Semi-)Definite
#'
#' `isSPSD()` throws an error message if the argument is not a non-empty symmetric positive (semi-)definite matrix without NA entries and not nilpotent. The error message explains which property is not respected.
#'
#' @param M A matrix object.
#'
#' @return NULL
#'
#'
#' @examples
#' isSPSD(diag(10)) # NULL
#' try(isSPSD(diag(0))) # error
#' try(isSPSD(diag(NA, ncol = 3, nrow = 3))) # error
#' try(isSPSD(matrix(1, ncol = 3, nrow = 2))) # error
isSPSD <- function(M) {
  if (!is.matrix(M)) {
    stop("`M` must be a matrix, not a ", class(M), ".")
  }

  if (!ncol(M) == nrow(M)) {
    stop("`M` must be a square matrix. Instead, it has dimension ", nrow(M), "x", ncol(M), ".")
  }

  if (ncol(M) == 0) {
    stop("`M` must be a non-empty matrix.")
  }

  if (anyNA(M)) {
    stop("`M` must not contain NA values.")
  }

  if (!isSymmetric(M)) {
    stop("`M` must be symmetric.")
  }

  V <- eigen(M, symmetric = T, only.values = T)$values

  if (any(V < -0.001) | all(V == 0)) {
    stop("`M` must be a positive definite or semi-definite matrix, but not nilpotent.")
  }
}
