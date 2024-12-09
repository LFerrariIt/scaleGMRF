#' Generalized Inverse of a Matrix
#'
#' `gen_inv()` calculates the Moore-Penrose generalized inverse of a symmetric positive (semi-)definite matrix, with the possibility of explicitly stating the rank-deficiency when the matrix is positive semi-definite.
#'
#' @param Q A matrix, which must be symmetric positive (semi-)definite.
#' @param rank_def A non-negative integer, indicating the rank deficiency of the matrix. If `NULL`, the rank deficiency is estimated.
#'
#'
#' @examples
#' gen_inv(diag(10), rank_def = 0)
#' gen_inv(as.matrix(spam::precmat.RW1(10)), rank_def = 1)
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
    stop("rank_def must be a non-negative integer.")
  }

  if (rank_def %% 1 != 0 | rank_def < 0) {
    stop("rank_def must be a non-negative integer.")
  }

  gen_inverse <- U %*% diag(c(1 / (V[1:(nrow(Q) - rank_def)]), rep(0, rank_def))) %*% t(U)

  return(gen_inverse)
}
