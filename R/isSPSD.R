isSPSP <- function(M) {
  M_name <- deparse(substitute(M))

  if (!is.matrix(M)) {
    stop(paste0("`", M_name, "` must be a matrix, not a ", class(M), "."))
  }

  if (!ncol(M) == nrow(M)) {
    stop(paste0("`", M_name, "` must be a square matrix, not a matrix of dimension ", nrow(M), "x", ncol(M), "."))
  }

  if (ncol(M) == 0) {
    stop(paste("`", M_name, "` must be a non-empty matrix."))
  }

  if (anyNA(M)) {
    stop(paste("`", M_name, "` must not contain NA values."))
  }

  if (!isSymmetric(M)) {
    stop(paste0("`", M_name, "` must be symmetric."))
  }

  V <- eigen(M, symmetric = T, only.values = T)$values

  if (any(V < 0) | all(V == 0)) {
    stop(paste("`M` must be a positive semi-definite matrix."))
  }
}

isValidBasis <- function(D, Q) {
  D_name <- deparse(substitute(D))
  Q_name <- deparse(substitute(Q))

  if (!is.matrix(D)) {
    stop(paste("`", D_name, "` must be a matrix, not a", class(D), "."))
  }

  if (nrow(D) == 0 | ncol(D) == 0) {
    stop(paste("`", D_name, "` must be a non-empty matrix."))
  }

  if (anyNA(D)) {
    stop(paste("`", D_name, "` must not contain NA values."))
  }

  if (!ncol(D) == ncol(Q)) {
    stop(paste0("The number of columns of `", Q_name, "` and `", D_name, "` must be equal. Instead,", ncol(Q), "!=", ncol(D), "."))
  }
}
