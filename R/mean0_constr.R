prec_under_mean0 <- function(Q,D=NULL, rank_def = NULL) {

  isSPSP(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D)

  if (!ncol(D) == ncol(Q)) {
    stop(paste0("The number of columns of `Q` and `D` must be equal. Instead,", ncol(Q), "!=", ncol(D), "."))
  }

  Sigma <- gen_inverse_func(Q, rank_def = rank_def)

  A <- matrix(colMeans(D),nrow=1)

  Sigma_constr <- Sigma-Sigma%*%t(A)%*%A%*%Sigma/(A%*%Sigma%*%t(A))

  Q_constr <- gen_inverse_func(Q, rank_def=NULL)

  isSPSP(Q_constr)

  return(Q_constr)
}
