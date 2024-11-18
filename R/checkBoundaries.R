#' Returns the boundaries for the covariates
#'
#' `checkBoundaries()` tests if `x` is a numeric vector with no missing values, throwing an error if it is not, and returns the boundaries of the interval if all entries of `x` lie within this interval.
#'
#' @param x A numeric vector.
#' @param m A number indicating the lower boundary of the support of X. If not provided, it is set to the minimum value from `x`.
#' @param M A number indicating the upper boundary of the support of X. If not provided, it is set to the maximum value from `x`.
#'
#' @return A vector of 2 numeric entries, with the first being smaller than the second.
#'
#' @examples
#' x <- runif(100)
#' checkBoundaries(x)
#' try(checkBoundaries(x, m = 2, M = 3))
#' try(checkBoundaries(x, m = 1, M = 0))
checkBoundaries <- function(x, m = NULL, M = NULL) {
  if (!is.numeric(x) | !is.vector(x)) {
    stop("`x` must be a numeric vector.")
  }

  if (any(is.na(x))) {
    stop("`x` must not contain NA values.")
  }

  if (is.null(m)) {
    m <- min(x)
  }

  if (is.null(M)) {
    M <- max(x)
  }

  if (!is.numeric(m) | !is.vector(m) | length(m) != 1 |
    !is.numeric(M) | !is.vector(M) | length(M) != 1) {
    stop(
      "`m` and `M` must be numbers or NULL. Instead m=",
      paste0(m, collapse = ","), " and M=",
      paste0(M, collapse = ","), "."
    )
  }

  if (m >= M) {
    stop("`M` must be larger than `m`. Instead, m=", m, " and M=", M, ".")
  }

  if (any(x < m) | any(x > M)) {
    stop("`x` must lie between m=", m, " and M=", M, ".")
  }

  return(c(m, M))
}
