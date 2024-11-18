#' @rdname rw1_standard
#'

rw2_standard <- function(x) {
  isOrderedFactor(x)

  K <- length(levels(x))
  X_dist <- 1:K
  D_dist <- diag(K)
  D <- as.matrix(
    fastDummies::dummy_cols(x,
      remove_selected_columns = T, ignore_na = T
    )
  )
  Q <- as.matrix(spam::precmat.RW2(K))
  S <- cbind(rep(1, K), 1:K)
  C <- scale_GMRF(Q = Q, D = NULL, rank_def = 2)

  return(list(
    "precision" = Q * C,
    "basis" = D,
    "scaling_constant" = C,
    "null_space" = S,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
