#' Standardized random walk effects
#'
#' `rw_standard()` provides a list of elements to build a standardized random walk effect of order 1 or 2.
#'
#' @param x An ordered factor, with more than 1 level.
#' @param order Either 1 or 2, representing the order of the random walk.
#'
#' @inherit iid_standard return
#'
#' @examples
#' K <- 20
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K))
#'
#' # Example for a random walk of order 1 --------------------
#' rw_standard(x, order = 1)
#'
#' # Example for a random walk of order 2 --------------------
#' rw_standard(x, order = 2)
rw_standard <- function(x, order) {
  isOrderedFactor(x)

  K <- length(levels(x))
  X_dist <- 1:K
  D_dist <- diag(K)
  D <- as.matrix(
    fastDummies::dummy_cols(x,
      remove_selected_columns = T, ignore_na = T
    )
  )
  if (order == 1) {
    Q <- as.matrix(spam::precmat.RW1(K))
    S <- matrix(1, ncol = 1, nrow = K)
  } else if (order == 2) {
    if (K < 4) {
      stop("The number of levels 'K' should be larger than 3.")
    }
    Q <- as.matrix(spam::precmat.RW2(K))
    S <- cbind(rep(1, K), 1:K)
  } else {
    stop("'order' can only be equal to 1 or 2. Instead, order=", order, ".")
  }

  C <- scale_GMRF(Q = Q, D = NULL, rank_def = order)

  return(list(
    "Q" = Q * C,
    "D" = D,
    "C" = C,
    "null_space" = S,
    "X_distribution" = X_dist,
    "basis_distribution" = D_dist
  ))
}
