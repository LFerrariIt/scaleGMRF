#' Standardized P-Spline effects
#'
#' @description
#' `pspline_standard()` provides a list of elements to build a standardized P-Spline effect with an IGMRF prior of order 1 or 2 on the coefficients.
#'
#' @inheritParams  bspline
#' @param order Either 1 or 2, representing the order of the random walk.
#' @param sparse_sol Logical, indicating whether the output is to be provided via sparse matrices. Otherwise, the solution is provided in a non-sparse format. By default, TRUE.
#'
#' @inherit iid_standard return
#'
#' @details This function is used in the `f_Xunif()` function of the `scaleGMRF` package to provide the elements to build a standardized P-Spline effect. This function provides a modified version of P-Spline effects. See more details in `vignette("psplines", package="scaleGMRF")`.
#'
#' @examples
#' x <- seq(0, 1, length.out = 100)
#'
#' pspline_standard(x, K = 20, order = 1)
#'
pspline_standard <- function(x, K, order, m = NULL, M = NULL, sparse_sol = T) {
  m_M <- xBoundaries(x = x, m = m, M = M)
  m <- m_M[1]
  M <- m_M[2]

  B <- bspline(x, K = K, m = m, M = M)
  x_unif <- seq(m, M, length.out = 1000)
  B_unif <- bspline(x_unif, K = K, m = m, M = M)

  # Definition of original matrices

  if (order == 1) {
    Q <- as.matrix(spam::precmat.RW1(K))
  } else if (order == 2) {
    Q <- as.matrix(spam::precmat.RW2(K))
  } else {
    stop("'order' can only be equal to 1 or 2. Instead, order=", order, ".")
  }

  G <- diag(diag(Q))
  W <- G - Q

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

  if (order == 1) {
    S_tilde <- S_tilde[, 1]
  }

  # Design of the new process for a given choice of symmetric lambdas
  tilde_f_fun <- function(lambdas) {
    lambdas <- abs(lambdas)
    # Lambda matrix
    if (K %% 2 == 0) {
      Lambda <- diag(c(lambdas, rev(lambdas)))
    } else {
      Lambda <- diag(c(lambdas, rev(lambdas)[-1]))
    }
    # Null space of R matrix
    new_S <- Lambda %*% S_tilde
    # New W,G,R,Q matrices
    if (order == 1) {
      W_tilde <- W / (new_S %*% t(new_S))
      G_tilde <- diag(as.vector(W_tilde %*% new_S / new_S))
    } else if (order == 2) {
      W_tilde <- matrix(0, nrow = K, ncol = K)
      G_tilde <- matrix(0, nrow = K, ncol = K)
      for (k in 1:K) {
        for (l in 1:K) {
          W_tilde[k, l] <- (l - k) * W[k, l] / (new_S[k, 1] * new_S[l, 2] - new_S[k, 2] * new_S[l, 1])
        }
        W_tilde[k, k] <- 0
        G_tilde[k, k] <- W_tilde[k, ] %*% new_S[, 1] / new_S[k, 1]
      }
    }

    R_tilde <- G_tilde - W_tilde

    return(list(
      "R_tilde" = R_tilde,
      "Lambda" = Lambda,
      "new_S" = new_S
    ))
  }

  # Function computing the KLD for a given choice of lambdas
  kld_fun <- function(lambdas) {
    # New process
    tilde_f <- tilde_f_fun(lambdas)
    R_tilde <- tilde_f$R_tilde
    Lambda <- tilde_f$Lambda
    new_S <- tilde_f$new_S
    # Covariance matrix
    Sigma_tilde <- Lambda %*%
      gen_inv(R_tilde, rank_def = order) %*%
      Lambda
    # only the non-constant part of KLD wrt to lambdas
    kld <- sum(colSums(Q * Sigma_tilde)) -
      sum(log(eigen(Sigma_tilde)$values[1:(K - order)]))
    return(kld)
  }
  # Optimization of the KLD function with symmetric entries for lambda
  optimal_lambdas <- stats::nlm(kld_fun, rep(1, ceiling(K / 2)), print.level = 2)$estimate
  # Process for the optimal lambdas
  tilde_f <- tilde_f_fun(optimal_lambdas)
  R_tilde <- tilde_f$R_tilde
  Lambda <- tilde_f$Lambda
  new_S <- tilde_f$new_S
  # Scaling constant
  scaling_constant <- scale_GMRF(
    R_tilde,
    B_unif %*% Lambda,
    rank_def = order
  )

  if (sparse_sol) {
    # Sparse version
    basis <- B %*% Lambda
    precision <- R_tilde * scaling_constant
    null_space <- new_S
    basis_dist <- B_unif %*% Lambda
  } else {
    # Dense version
    Q_tilde <- gen_inv(Lambda %*%
      gen_inv(R_tilde, rank_def = order) %*%
      Lambda, rank_def = order)

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
