ex_scale <- function(Q,D=NULL, rank_def = NULL) {

  isSPSP(Q)

  if (is.null(D)) {
    D <- diag(nrow(Q))
  }

  isValidBasis(D,Q)

  if (fixed) {
    Q <- prec_under_mean0(Q,D,rank_def)
  }

  C <- function(Q,D,rank_def=NULL)

  Q_standardized <- Q*C

  return(Q_standardized)
}
