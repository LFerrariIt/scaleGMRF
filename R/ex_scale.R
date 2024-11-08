ex_scale <- function(Q, D = NULL, rank_def = NULL) {
  isSPSP(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D, Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  var_diag <- rowSums(D * (D %*% Sigma))

  if (any(var_diag < 0)) {
    stop(paste0("Invalid case: the conditional variance entries cannot be negative."))
  }

  if (all(var_diag == 0)) {
    stop(paste0("Invalid case: the conditional variance entries cannot be all null."))
  }

  C <- mean(var_diag)

  names(C) <- "scaling constant"

  return(C)
}
