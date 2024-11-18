#' Standardized P-Spline effects
#'
#' @description
#' `Pspline1_standard()` provides a list of elements to build a standardized P-Spline effect with a random walk process of order 1 on the coefficients. Note that the modified version is used, as detailed in ???.
#'
#' `Pspline2_standard()` provides a list of elements to build a standardized P-Spline effect with a random walk process of order 2 on the coefficients. Note that the modified version is used, as detailed in ???.
#'
#' @param x A numeric vector.
#' @param m A number indicating the lower boundary of the support of X. If not provided, it is set to the minimum value from `x`.
#' @param M A number indicating the upper boundary of the support of X. If not provided, it is set to the maximum value from `x`.
#' @param K A positive integer larger than 3, specifying the number of basis functions.
#' @param sparse_sol Logical, indicating whether the solution is to be provided in the sparse version. Otherwise, the solution is provided in a non-sparse format. By default, TRUE.
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
#' @examples
#' x <- seq(0, 1, length.out = 100)
#'
#' Pspline1_standard(x, K = 10)
#' Pspline1_standard(x, K = 20)
#' Pspline2_standard(x, K = 10)
#' Pspline2_standard(x, K = 20)
#'
Pspline1_standard <- function(x, K, m = NULL, M = NULL, sparse_sol = T) {
  # Error messages on the arguments

  if (!is.numeric(K) | !is.vector(K) | length(K) != 1) {
    stop("`K` must be a number.")
  }

  if (K < 4 | K %% 1 != 0) {
    stop("`K` must be an whole number larger than 3. Instead, K=", K, ".")
  }

  m_M <- checkBoundaries(x = x, m = m, M = M)
  m <- m_M[1]
  M <- m_M[2]

  B <- bspline((x - m) / (M - m), K = K)
  x_unif <- seq(m, M, length.out = 1000)
  B_unif <- bspline((x_unif - m) / (M - m), K = K)

  # P-Spline with Random Walk of order 1 on the coefficients

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
  results <- stats::nlm(optim_function, rep(1, ceiling(K / 2)), print.level = 2)
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
  scaling_constant <- scale_GMRF(R_tilde, B_unif %*% Lambda, rank_def = 1)

  if (sparse_sol) {
    basis <- B %*% Lambda
    precision <- R_tilde * scaling_constant
    null_space <- new_S
    basis_dist <- B_unif %*% Lambda
  } else {
    Q_tilde <- gen_inv(Lambda %*% gen_inv(R_tilde, rank_def = 1) %*% Lambda, rank_def = 1)

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
