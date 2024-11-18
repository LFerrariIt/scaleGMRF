#' Check if the argument is an ordered factor.
#'
#' `isOrderedFactor()` throws an error message if the argument is not an ordered factor with more than 2 levels or contains NA values.
#'
#' @param x An ordered factor, with more than 1 level.
#'
#' @return NULL
#'
#' @examples
#' K <- 20
#' x <- factor(1:K, ordered = TRUE, levels = c(1:K))
#' isOrderedFactor(x) # NULL
#' try(isOrderedFactor(c(x, NA))) # error
#' try(isOrderedFactor(factor(1:K, levels = c(1:K)))) # error
#' try(isOrderedFactor(factor(rep(1, K), levels = 1))) # error
isOrderedFactor <- function(x) {
  if (!is.ordered(x)) {
    stop("`x` must be an ordered factor. Instead, it is a", class(x), ".")
  }

  if (any(is.na(x))) {
    stop("`x` must not contain NA values.")
  }

  if (length(levels(x)) < 2) {
    stop("The number of levels of `x` must be larger than 1. Instead, it is equal to", length(levels(x)), ".")
  }
}
