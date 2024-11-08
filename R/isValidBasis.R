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
