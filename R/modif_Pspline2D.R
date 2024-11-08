modif_Pspline2D <- function(x, K1, K2, sparse_sol = T) {
  # Error messages on the arguments

  if (!is.numeric(K1) | !is.numeric(K2)) {
    stop(paste0("`K1` and `K2` must be numeric."))
  }

  if (K1 < 4 | K2 < 4) {
    stop(paste0("`K1` and `K2` must be larger than 3, instead K1 = ", K1, " and K2 = ", K2, "."))
  }

  if (!is.matrix(x)) {
    stop(paste0("`x` must be a matrix."))
  }

  if (ncol(x) != 2 | nrow(x) == 0) {
    stop(paste0("`x` must have 2 columns and a positive number of rows. Instead, it has dimension ", nrow(x), "x", ncol(x), "."))
  }

  if (any(x < 0) | any(x > 1) | any(is.na(x))) {
    stop(paste0("The entries of `x` must be normalized to lie between 0 and 1 and must not contain NA values."))
  }

  # Original Q, G, W matrices and generalized inverse of Q
  Q <- as.matrix(spam::precmat.IGMRFreglat(K1, K2, order = 1))
  G <- diag(diag(Q))
  W <- G - Q
  Sigma <- gen_inv(Q, rank_def = 1)
  # Computation of the S tilde matrix
  Delta1 <- 1 / (K1 - 3)
  if (K1 == 4) {
    s1 <- c(1 / 24, 11 / 24, 11 / 24, 1 / 24) * Delta1
  } else if (K1 == 5) {
    s1 <- c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) * Delta1
  } else if (K1 > 5) {
    s1 <- c(1 / 24, 1 / 2, 23 / 24, rep(1, K1 - 6), 23 / 24, 1 / 2, 1 / 24) * Delta1
  }
  Delta2 <- 1 / (K2 - 3)
  if (K2 == 4) {
    s2 <- c(1 / 24, 11 / 24, 11 / 24, 1 / 24) * Delta2
  } else if (K2 == 5) {
    s2 <- c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) * Delta2
  } else if (K2 > 5) {
    s2 <- c(1 / 24, 1 / 2, 23 / 24, rep(1, K2 - 6), 23 / 24, 1 / 2, 1 / 24) * Delta2
  }
  S_tilde <- kronecker(s2, s1)
  # Optimization of the KLD function with symmetric entries for lambda
  optim_function <- function(lambdas) {
    lambdas <- abs(lambdas)
    # Creation of the diagonal matrix Lambda
    lambdas1 <- lambdas[1:ceiling(K2 / 2)]
    lambdas2 <- lambdas[(ceiling(K2 / 2) + 1):(ceiling(K2 / 2) + ceiling(K1 / 2))]
    Lambda <- diag(kronecker(
      c(lambdas2, rev(lambdas2)),
      c(lambdas1, rev(lambdas1))
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
      sum(log(eigen(Sigma_tilde, only.values = T)$values[1:((K1 * K2) - 1)]))
    return(kld)
  }
  # Optimization of the KLD function with symmetric entries for lambda
  results <- stats::nlm(optim_function,
    c(rep(1, ceiling(K2 / 2)), rep(1, ceiling(K1 / 2))),
    print.level = 2
  )
  # Save the lambda values that minimize the KLD
  lambdas <- abs(results$estimate)
  # Lambda matrix
  lambdas1 <- lambdas[1:ceiling(K2 / 2)]
  lambdas2 <- lambdas[(ceiling(K2 / 2) + 1):(ceiling(K2 / 2) + ceiling(K1 / 2))]
  Lambda <- diag(kronecker(
    c(lambdas2, rev(lambdas2)),
    c(lambdas1, rev(lambdas1))
  ))
  # Null space for the R matrix
  new_S <- Lambda %*% S_tilde
  # New W,G,R,Q matrices
  W_tilde <- W / (new_S %*% t(new_S))
  G_tilde <- diag(diag(G) / c(new_S^2))
  R_tilde <- G_tilde - W_tilde

  # Basis matrix evaluated at x
  B_unif <- bspline2D(expand.grid(
    seq(0, 1, length.out = 200),
    seq(0, 1, length.out = 200)
  ), K1, K2)
  C <- ex_scale(R_tilde, B_unif %*% Lambda, rank_def = 1)
  B <- bspline2D(x[, 1], x[, 2], K1, K2)

  if (sparse_sol) {
    return(list(
      "basis matrix" = B %*% Lambda,
      "precision matrix" = R_tilde,
      "constraint matrix" = new_S,
      "scaling constant" = C
    ))
  } else {
    Q_tilde <- gen_inv(Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda, rank_def = 1)

    return(list(
      "basis matrix" = B,
      "precision matrix" = Q_tilde,
      "constraint matrix" = S_tilde,
      "scaling constant" = C
    ))
  }
}
