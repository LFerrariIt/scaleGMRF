#' Standardized 2D P-Spline effects
#'
#' `pspline_2D_standard()` provides a list of elements to build a standardized 2-dimensional P-Spline effect with an IGMRF prior of order 1 on the coefficients.
#'
#' @inheritParams bspline_2D
#' @param sparse_sol Logical, indicating whether the output is to be provided via sparse matrices. Otherwise, the solution is provided in a non-sparse format. By default, TRUE.
#'
#'
#' @inherit iid_standard return
#'
#' @details This function is used in the `f_Xunif()` function of the `scaleGMRF` package to provide the elements to build a standardized 2D P-Spline effect. This function provides a modified version of P-Spline effects. See more details in `vignette("psplines", package="scaleGMRF")`.
#'
#' @examples
#' x_grid <- as.matrix(expand.grid(
#'   seq(0, 1, length.out = 100),
#'   seq(0, 1, length.out = 100)
#' ))
#' pspline_2D_standard(x_grid, K = c(5, 5))
pspline_2D_standard <- function(x, K, m = NULL, M = NULL, sparse_sol = T) {
  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), ".")
  }

  if (is.null(m)) {
    m <- c(min(x[, 1]), min(x[, 2]))
  }
  if (is.numeric(m) & length(m) == 2) {
    m1 <- m[1]
    m2 <- m[2]
  } else {
    stop(
      "`m` must contain 2 numeric elements or be NULL. Instead, m=",
      paste0(m, collapse = ","), "."
    )
  }

  if (is.null(M)) {
    M <- c(max(x[, 1]), max(x[, 2]))
  }
  if (is.numeric(M) & length(M) == 2) {
    M1 <- M[1]
    M2 <- M[2]
  } else {
    stop(
      "`M` must contain 2 numeric elements. Instead, M=",
      paste0(M, collapse = ","), "."
    )
  }

  if (any(m >= M)) {
    stop(
      "`M` must be larger than `m`. Instead, m=",
      paste0(m, collapse = ","), " and M=",
      paste0(M, collapse = ","), "."
    )
  }

  B <- bspline_2D(x, K = K, m = m, M = M)

  x_unif <- as.matrix(
    expand.grid(
      seq(m[1], M[1], length.out = 200),
      seq(m[2], M[2], length.out = 200)
    )
  )
  B_unif <- bspline_2D(as.matrix(
    expand.grid(
      seq(0, 1, length.out = 200),
      seq(0, 1, length.out = 200)
    )
  ), K = K)


  # Original Q, G, W matrices and generalized inverse of Q
  Q <- as.matrix(spam::precmat.IGMRFreglat(K[1], K[2], order = 1))
  G <- diag(diag(Q))
  W <- G - Q
  Sigma <- gen_inv(Q, rank_def = 1)
  # Computation of the S tilde matrix
  Delta1 <- 1 / (K[1] - 3)
  if (K[1] == 4) {
    s1 <- c(1 / 24, 11 / 24, 11 / 24, 1 / 24) * Delta1
  } else if (K[1] == 5) {
    s1 <- c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) * Delta1
  } else if (K[1] > 5) {
    s1 <- c(1 / 24, 1 / 2, 23 / 24, rep(1, K[1] - 6), 23 / 24, 1 / 2, 1 / 24) * Delta1
  }
  Delta2 <- 1 / (K[2] - 3)
  if (K[2] == 4) {
    s2 <- c(1 / 24, 11 / 24, 11 / 24, 1 / 24) * Delta2
  } else if (K[2] == 5) {
    s2 <- c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) * Delta2
  } else if (K[2] > 5) {
    s2 <- c(1 / 24, 1 / 2, 23 / 24, rep(1, K[2] - 6), 23 / 24, 1 / 2, 1 / 24) * Delta2
  }
  S_tilde <- kronecker(s2, s1)

  # New process with symmetric lambdas
  tilde_f_fun <- function(lambdas) {
    lambdas <- abs(lambdas)
    # Creation of the diagonal matrix Lambda
    lambdas1 <- lambdas[1:K[1]]
    lambdas2 <- lambdas[(K[1] + 1):(K[1] + K[2])]
    Lambda <- diag(kronecker(
      c(lambdas2),
      c(lambdas1)
    ))
    # Null space for the R matrix
    new_S <- Lambda %*% S_tilde
    # New W,G,R,Q matrices
    W_tilde <- W / (new_S %*% t(new_S))
    G_tilde <- diag(diag(G) / c(new_S^2))
    R_tilde <- G_tilde - W_tilde

    return(list(
      "R_tilde" = R_tilde,
      "Lambda" = Lambda,
      "new_S" = new_S
    ))
  }
  # Optimization of the KLD function
  kld_fun <- function(lambdas) {
    tilde_f <- tilde_f_fun(lambdas)
    R_tilde <- tilde_f$R_tilde
    Lambda <- tilde_f$Lambda
    new_S <- tilde_f$new_S
    Sigma_tilde <- Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda
    # Computation of the KLD (only non-constant part)
    kld <- sum(colSums(Q * Sigma_tilde)) -
      sum(log(eigen(Sigma_tilde, only.values = T)$values[1:((K[1] * K[2]) - 1)]))
    return(kld)
  }
  # Optimization of the KLD function with symmetric entries for lambda
  optimal_lambdas <- stats::nlm(kld_fun,
    c(rep(1, K[1]), rep(1, K[2])),
    print.level = 2
  )$estimate
  # Save the lambda values that minimize the KLD
  tilde_f <- tilde_f_fun(optimal_lambdas)
  R_tilde <- tilde_f$R_tilde
  Lambda <- tilde_f$Lambda
  new_S <- tilde_f$new_S

  # Basis matrix evaluated at x
  scaling_constant <- scale_GMRF(R_tilde, B_unif %*% Lambda, rank_def = 1)

  if (sparse_sol) {
    # Sparse version
    basis <- B %*% Lambda
    precision <- R_tilde * scaling_constant
    null_space <- new_S
    basis_dist <- B_unif %*% Lambda
  } else {
    # Dense version
    Q_tilde <- gen_inv(Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda, rank_def = order)

    basis <- B
    precision <- Q_tilde * scaling_constant
    null_space <- S_tilde
    basis_dist <- B_unif
  }

  return(list(
    "precision" = precision,
    "basis" = basis,
    "scaling_constant" = scaling_constant,
    "null_space" = null_space,
    "X_distribution" = x_unif,
    "basis_distribution" = basis_dist
  ))
}
