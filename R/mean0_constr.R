prec_under_mean0 <- function(Q,D=NULL, rank_def = NULL) {

  isSPSP(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D,Q)

  Sigma <- gen_inv(Q, rank_def = rank_def)

  A <- matrix(colMeans(D),nrow=1)

  Sigma_constr <- Sigma-Sigma%*%t(A)%*%A%*%Sigma/(A%*%Sigma%*%t(A))

  Q_constr <- gen_inv(Q, rank_def=NULL)

  isSPSP(Q_constr)

  return(Q_constr)
}
