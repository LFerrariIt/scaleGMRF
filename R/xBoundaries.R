#' Valid boundaries for a numeric vector
#'
#' `xBoundaries()` returns valid boundaries for the specified numeric vector.
#'
#' @inheritParams linear_standard
#'
#' @return A numeric vector with 2 elements, the first being smaller than the second.
#'
#' @details This function is used in other functions of the `scaleGMRF` package to return appropriate boundaries for a numeric vector. If the boundaries are specified, the function simply checks that each entry of the numeric vector lies between the two extremes. Otherwise, the function first estimates the boundaries using the minimum and maximum values of the vector.
#'
#' @examples
#' # Valid example --------------------------------------------------
#' x <- runif(100)
#' xBoundaries(x)
#'
#' # Examples throwing errors ---------------------------------------
#' try(xBoundaries(x, m = 2, M = 3))
#' try(xBoundaries(x, m = 1, M = 0))
xBoundaries <- function(x, m = NULL, M = NULL) {
  if (!is.numeric(x) | !is.vector(x)) {
    stop("`x` must be a numeric vector. Instead, x is a ", class(x), ".")
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
