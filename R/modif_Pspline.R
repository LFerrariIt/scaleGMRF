modif_Pspline <- function(x, K, order = 2, sparse_sol = T) {
  # Error messages on the arguments

  if (!is.numeric(K)) {
    stop(paste0("`K` must be numeric."))
  }

  if (K < 4) {
    stop(paste0("`K` must be at least than 4, instead K=", K, "."))
  }

  if (!is.numeric(x)) {
    stop(paste0("`x` must be a numeric vector."))
  }

  if (any(x < 0) | any(x > 1)| any(is.na(x))) {
    stop(paste0("The vector `x` must be normalized to lie between 0 and 1 and not contain NA values."))
  }

  # P-Spline with Random Walk of order 1 on the coefficients
  if (order == 1) {
    # Definition of original matrices
    Q <- as.matrix(spam::precmat.RW1(n = K))
    G <- diag(diag(Q))
    W <- G - Q

    # Computation of the S tilde matrix for different values of K
    if (K == 4) {
      S_tilde <- c(1 / 24, 11 / 24, 11 / 24, 1 / 24) / (K - 3)
    } else if (K == 5) {
      S_tilde <- c(1 / 24, 1 / 2, 22 / 24, 1 / 2, 1 / 24) / (K - 3)
    } else if (K > 5) {
      S_tilde <- c(1 / 24, 1 / 2, 23 / 24, rep(1, K - 6), 23 / 24, 1 / 2, 1 / 24) / (K - 3)
    }
    # Function computing the KLD for a given choice of lambdas
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
      W_tilde <- W / (new_S %*% t(new_S))
      G_tilde <- diag(as.vector(W_tilde %*% new_S / new_S))
      R_tilde <- G_tilde - W_tilde
      Sigma_tilde <- Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda
      # Computation of the KLD (only the non-constant part wrt to lambda)
      kld <- sum(colSums(Q * Sigma_tilde)) -
        sum(log(eigen(Sigma_tilde)$values[1:(K - 1)]))
      return(kld)
    }
    # Optimization of the KLD function with symmetric entries for lambda
    results <- nlm(optim_function, rep(1, ceiling(K / 2)), print.level = 2)
    # Save the lambda values that minimize the KLD
    lambdas <- abs(results$estimate)
    # Lambda matrix
    if (K %% 2 == 0) {
      Lambda <- diag(c(lambdas, rev(lambdas)))
    } else {
      Lambda <- diag(c(lambdas, rev(lambdas)[-1]))
    }
    # Null space of R matrix
    new_S <- Lambda %*% S_tilde
    # New W,G,R,Q matrices
    W_tilde <- W / (new_S %*% t(new_S))
    G_tilde <- diag(as.vector(W_tilde %*% new_S / new_S))
    R_tilde <- G_tilde - W_tilde

    # Basis matrix evaluated at x
    B_unif <- bspline(seq(0, 1, length.out = 1000), N_basis = K)
    C <- ex_scale(R_tilde, B_unif %*% Lambda, rank_def = order)
    B <- bspline(x, N_basis = K)

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

  # P-Spline with Random Walk of order 2 on the coefficients
  if (order == 2) {
    # Original Q, G, W matrices and generalized inverse of Q
    Q <- as.matrix(spam::precmat.RW2(n = K))
    G <- diag(diag(Q))
    W <- G - Q
    Sigma <- gen_inv(Q, rank_def = order)
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
      Sigma_tilde <- Lambda %*% gen_inv(R_tilde, rank_def = order) %*% Lambda
      # Computation of the KLD (only the non-constant part wrt to lambda)
      kld <- sum(colSums(Q * Sigma_tilde)) -
        sum(log(eigen(Sigma_tilde)$values[1:(K - 2)]))
      return(kld)
    }
    # Optimization of the KLD function with symmetric entries for lambda
    results <- nlm(optim_function, rep(1, ceiling(K / 2)), print.level = 2)
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
    B_unif <- bspline(seq(0, 1, length.out = 1000), N_basis = K)
    C <- ex_scale(R_tilde, B_unif %*% Lambda, rank_def = order)
    B <- bspline(x, N_basis = K)

    if (sparse_sol) {
      return(list(
        "basis matrix" = B %*% Lambda,
        "precision matrix" = R_tilde,
        "constraint matrix" = new_S,
        "scaling constant" = C
      ))
    } else {
      Q_tilde <- gen_inv(Lambda %*% gen_inv(R_tilde, rank_def = order) %*% Lambda, rank_def = order)

      return(list(
        "basis matrix" = B,
        "precision matrix" = Q_tilde,
        "constraint matrix" = S_tilde,
        "scaling constant" = C
      ))
    }
  }
}
