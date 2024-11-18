#' @rdname Pspline1_standard
#'
Pspline2_standard <- function(x, K, m = NULL, M = NULL, sparse_sol = T) {
  # Error messages on the arguments

  if (!is.numeric(K) | !is.vector(K) | length(K) != 1) {
    stop("`K` must be a number.")
  }

  if (K < 4 | K %% 1 != 0) {
    stop("`K` must be an whole number larger than 4, instead K=", K, ".")
  }

  m_M <- checkBoundaries(x = x, m = m, M = M)
  m <- m_M[1]
  M <- m_M[2]

  B <- bspline((x - m) / (M - m), K = K)
  x_unif <- seq(m, M, length.out = 1000)
  B_unif <- bspline((x_unif - m) / (M - m), K = K)

  # P-Spline with Random Walk of order 2 on the coefficients
  # Original Q, G, W matrices and generalized inverse of Q
  Q <- as.matrix(spam::precmat.RW2(n = K))
  G <- diag(diag(Q))
  W <- G - Q
  Sigma <- gen_inv(Q, rank_def = 2)
  # Computation of the S tilde matrix for different values of K
  Delta <- 1 / (K - 3)
  if (K == 4) {
    S_tilde <- cbind(
      c(1 / 24, 11 / 24, 11 / 24, 1 / 24) * Delta,
      c(1 / 120, 11 / 60, 11 / 40, 1 / 30)
    )
  } else if (K == 5) {
    S_tilde <- cbind(
      c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) * Delta,
      c(1 / 480, 7 / 120, 11 / 48, 23 / 120, 3 / 160)
    )
  } else if (K == 6) {
    S_tilde <- cbind(
      c(1 / 24, 1 / 2, 23 / 24, rep(1, K - 6), 23 / 24, 1 / 2, 1 / 24) * Delta,
      c(
        Delta / 120,
        14 * Delta / 60,
        121 * Delta / 120,
        23 / 24 - 121 * Delta / 120,
        1 / 2 - 14 * Delta / 60,
        1 / 24 - Delta / 120
      ) * Delta
    )
  } else if (K >= 7) {
    S_tilde <- cbind(
      c(1 / 24, 1 / 2, 23 / 24, rep(1, K - 6), 23 / 24, 1 / 2, 1 / 24) * Delta,
      c(
        Delta / 120,
        14 * Delta / 60,
        121 * Delta / 120,
        c(2:(K - 5)) * Delta,
        23 / 24 - 121 * Delta / 120,
        1 / 2 - 14 * Delta / 60,
        1 / 24 - Delta / 120
      ) * Delta
    )
  }
  # Optimization of the KLD function with symmetric entries for lambdas
  optim_function <- function(lambdas) {
    lambdas <- abs(lambdas)
    # Creation of the diagonal matrix Lambda
    if (K %% 2 == 0) {
      Lambda <- diag(c(lambdas, rev(lambdas)))
    } else {
      Lambda <- diag(c(lambdas, rev(lambdas)[-1]))
    }
    # Null space for the R matrix
    new_S <- Lambda %*% S_tilde
    # New W, G, R matrix and generalized inverse of new Q
    W_tilde <- matrix(0, nrow = K, ncol = K)
    G_tilde <- matrix(0, nrow = K, ncol = K)
    for (k in 1:K) {
      for (l in 1:K) {
        W_tilde[k, l] <- (l - k) * W[k, l] / (new_S[k, 1] * new_S[l, 2] - new_S[k, 2] * new_S[l, 1])
      }
      W_tilde[k, k] <- 0
      G_tilde[k, k] <- W_tilde[k, ] %*% new_S[, 1] / new_S[k, 1]
    }
    R_tilde <- G_tilde - W_tilde
    Sigma_tilde <- Lambda %*% gen_inv(R_tilde, rank_def = 2) %*% Lambda
    # Computation of the KLD (only the non-constant part wrt to lambda)
    kld <- sum(colSums(Q * Sigma_tilde)) -
      sum(log(eigen(Sigma_tilde)$values[1:(K - 2)]))
    return(kld)
  }
  # Optimization of the KLD function with symmetric entries for lambda
  results <- stats::nlm(optim_function, rep(1, ceiling(K / 2)), print.level = 2)
  # Save the lambda values that minimize the KLD
  lambdas <- abs(results$estimate)
  # Lambda matrix
  if (K %% 2 == 0) {
    Lambda <- diag(c(lambdas, rev(lambdas)))
  } else {
    Lambda <- diag(c(lambdas, rev(lambdas)[-1]))
  }
  # Null space for R
  new_S <- Lambda %*% S_tilde
  # New W,G,R,Q matrices
  W_tilde <- matrix(0, nrow = K, ncol = K)
  G_tilde <- matrix(0, nrow = K, ncol = K)
  for (k in 1:K) {
    for (l in 1:K) {
      W_tilde[k, l] <- (l - k) * W[k, l] / (new_S[k, 1] * new_S[l, 2] - new_S[k, 2] * new_S[l, 1])
    }
    W_tilde[k, k] <- 0
    G_tilde[k, k] <- W_tilde[k, ] %*% new_S[, 1] / new_S[k, 1]
  }
  R_tilde <- G_tilde - W_tilde

  # Basis matrix evaluated at x
  scaling_constant <- scale_GMRF(R_tilde, B_unif %*% Lambda, rank_def = 2)

  if (sparse_sol) {
    basis <- B %*% Lambda
    precision <- R_tilde * scaling_constant
    null_space <- new_S
    basis_dist <- B_unif %*% Lambda
  } else {
    Q_tilde <- gen_inv(Lambda %*% gen_inv(R_tilde, rank_def = 2) %*% Lambda, rank_def = 2)

    basis <- B
    precision <- Q_tilde * scaling_constant
    null_space <- S_tilde
    basis_dist <- B_unif
  }

  return(list(
    "basis" = basis,
    "precision" = precision,
    "scaling_constant" = scaling_constant,
    "null_space" = null_space,
    "X_distribution" = x_unif,
    "basis_distribution" = basis_dist
  ))
}
