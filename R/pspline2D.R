#' Standardized 2D P-Spline effect
#'
#' `pspline_2D_standard()` provides a list of elements to build a standardized 2-dimensional P-Spline effect. Specifically, the scaled precision matrix, the basis matrix, the null space, and the scaling constant are provided.
#'
#' @param x A matrix of 2 columns of numeric entries.
#' @param K A vector containing 2 positive integers larger than 3, specifying the number of basis functions respectively for the first and second dimensions.
#' @param m A vector of 2 numbers, indicating the lower extremes in the two dimensions. If not provided, it is set to the minimum values from the `x` columns.
#' @param M A vector of 2 numbers, indicating the upper extremes in the two dimensions. If not provided, it is set to the maximum values from the `x` columns.
#' @param sparse_sol Logical, indicating whether the solution is to be provided in its sparse version. Otherwise, the solution is provided in a non-sparse format. By default, TRUE.
#'
#'
#' @return A list of 6 elements, containing:
#' * `precision`: precision matrix
#' * `basis`:  basis matrix evaluated at `x`
#' * `scaling_constant`: a positive number, representing the appropriate scaling constant
#'  * `null_space`: a matrix, representing the null space of the precision matrix
#'  * `X_distribution`: a matrix with 2 columns of values sampled from the Uniform distribution of X
#'  * `basis_distribution`: basis matrix evaluated at `X_distribution`
#' @examples
#' x_grid <- as.matrix(expand.grid(
#'   seq(0, 1, length.out = 100),
#'   seq(0, 1, length.out = 100)
#' ))
#' pspline_2D_standard(x_grid, K = c(5, 5))
pspline_2D_standard <- function(x, K, m = NULL, M = NULL, sparse_sol = T) {
  # Error messages on the arguments

  if (!is.numeric(K) | !is.vector(K) | length(K) != 2) {
    stop("`K` must be a vector containing 2 numbers. Instead, K=", K, ".")
  }

  if (any(K < 4)) {
    stop("The entries of `K` must be larger than 3. Instead, K = ", K, ".")
  }

  if (!is.matrix(x)) {
    stop("`x` must be a matrix.")
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), ".")
  }

  if (any(is.na(x))) {
    stop("`x` must not contain NA values.")
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

  x_norm <- cbind(
    (x[, 1] - m[1]) / (M[1] - m[1]),
    (x[, 2] - m[2]) / (M[2] - m[2])
  )

  if (any(x_norm < 0) | any(x_norm > 1) | any(is.na(x_norm))) {
    stop("The entries of `x` must lie between the given boundaries and must not contain NA values.")
  }

  B <- bspline_2D(x, K=K, m = m, M = M)


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
  # Optimization of the KLD function with symmetric entries for lambda
  optim_function <- function(lambdas) {
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
    Sigma_tilde <- Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda
    # Computation of the KLD (only non-constant part)
    kld <- sum(colSums(Q * Sigma_tilde)) -
      sum(log(eigen(Sigma_tilde, only.values = T)$values[1:((K[1] * K[2]) - 1)]))
    return(kld)
  }
  # Optimization of the KLD function with symmetric entries for lambda
  results <- stats::nlm(optim_function,
    c(rep(1, K[1]), rep(1, K[2])),
    print.level = 2
  )
  # Save the lambda values that minimize the KLD
  lambdas <- abs(results$estimate)
  # Lambda matrix
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

  # Basis matrix evaluated at x
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
  ), K=K)

  # Basis matrix evaluated at x
  scaling_constant <- scale_GMRF(R_tilde, B_unif %*% Lambda, rank_def = 1)

  if (sparse_sol) {
    basis <- B %*% Lambda
    precision <- R_tilde * scaling_constant
    null_space <- new_S
    basis_dist <- B_unif %*% Lambda
  } else {
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
