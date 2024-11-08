gen_inv <- function(M, rank_def = 1) {

  isSPSP(M)

  U <- eigen(M, symmetric = T)$vectors
  V <- eigen(M, symmetric = T)$values

  if (is.null(rank_def)) {
    rank_def <- nrow(M) - Matrix::rankMatrix(M)
  }

  gen_inverse <- U %*% diag(c(1 / (V[1:(nrow(M) - rank_def)]), rep(0, rank_def))) %*% t(U)

  return(gen_inverse)
}
