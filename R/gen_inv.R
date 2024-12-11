#' Generalized Inverse of a Matrix
#'
#' `gen_inv()` calculates the Moore-Penrose generalized inverse of a Symmetric Positive (Semi-)Definite (SPSD) matrix, with the possibility of explicitly stating the rank-deficiency when the matrix is positive semi-definite.
#'
#' @param Q An SPSD matrix.
#' @param rank_def A non-negative integer, indicating the rank deficiency of the matrix. If `NULL`, the rank deficiency is estimated.
#'
#' @details This function is used in other functions of the `scaleGMRF` package to compute the inverse of rank-deficient precision matrices, for which the rank-deficiency is known and does not need to be estimated, as it happens in the function [MASS::ginv()].
#'
#' @examples
#' # Example with a diagonal matrix -----------------------------------
#' gen_inv(diag(10), rank_def = 0)
#'
#' # Example with a random walk of order 1 precision matrix ------------
#' gen_inv(as.matrix(spam::precmat.RW1(10)), rank_def = 1)
#'
#' # Example with a random walk of order 2 precision matrix ------------
#' gen_inv(as.matrix(spam::precmat.RW2(10)), rank_def = 2)
#'
#' @seealso [solve()],[MASS::ginv()]

gen_inv <- function(Q, rank_def = NULL) {
  isSPSD(Q)

  U <- eigen(Q, symmetric = T)$vectors
  V <- eigen(Q, symmetric = T)$values

  if (is.null(rank_def)) {
    rank_def <- nrow(Q) - Matrix::rankMatrix(Q)
  }

  if (!is.numeric(rank_def)) {
    stop("`rank_def` must be a non-negative integer, not a ", class(rank_def), ".")
  }

  if (rank_def %% 1 != 0 | rank_def < 0) {
    stop("`rank_def` must be a non-negative integer. Instead, rank_def=", rank_def, " .")
  }

  gen_inverse <- U %*% diag(x = c(
    1 / (V[1:(nrow(Q) - rank_def)]),
    rep(0, rank_def)
  ), nrow = ncol(Q)) %*% t(U)

  return(gen_inverse)
}
