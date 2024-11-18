#' Standardized Besag or ICAR effects
#'
#' #'@description
#' `besag_standard()` provides a list of elements to build a standardized Besag effect.
#'
#' @param x An ordered factor, with more than 1 level.
#' @param adj_mat Matrix needed for `model="besag"`, representing the adjacency matrix from the lattice.
#'
#' @return A list of 6 elements, containing:
#' \itemize{
#' \item{`precision`: precision matrix.}
#' \item{`basis`:  basis matrix evaluated at `x`}
#' \item{`scaling_constant`: a positive number, representing the appropriate scaling constant.}
#'  \item{`null_space`: a matrix, representing the null space of the precision matrix.}
#'  \item{`X_distribution`: a numeric vector sampled from the Uniform distribution on X.}
#'  \item{`basis_distribution`: basis matrix evaluated at `X_distribution`.}
#' }
#'
besag_standard <- function(x, adj_mat) {
  isOrderedFactor(x)

  K <- length(levels(x))
  X_dist <- 1:K
  D_dist <- diag(K)
  D <- as.matrix(
    fastDummies::dummy_cols(x,
      remove_selected_columns = T, ignore_na = T
    )
  )

  if (!is.matrix(adj_mat)) {
    stop("`adj_mat` must be a matrix. Instead, it is a", class(adj_mat), ".")
  }

  if (!isSymmetric(adj_mat)) {
    stop("`adj_mat` must be a symmetric matrix.")
  }

  if (ncol(adj_mat) != K) {
    stop("The dimension of `adj_mat` must correspond to the number of levels in `x`. Instead, ncol(adj_mat)=", ncol(adj_mat), "!=length(levels(x))=", length(levels(x)), ".")
  }

  Q <- as.matrix(spam::precmat.IGMRFirreglat(adj_mat))
  S <- matrix(1, ncol = 1, nrow = K)
  C <- scale_GMRF(Q = Q, D = NULL, rank_def = 1)

  return(list(
    "precision" = Q * C,
    "basis" = D,
    "scaling_constant" = C,
    "null_space" = S,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
