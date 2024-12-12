#' Check if the argument is an ordered factor.
#'
#' `isOrderedFactor()` throws an error message if the argument is not an ordered factor with more than 2 levels or contains NA values.
#'
#' @param x An ordered factor, with more than 1 level.
#'
#' @return NULL
#'
#' @details This function is used in other functions of the `scaleGMRF` package to check whether their first argument is an ordered factor.
#'
#' @examples
#' K <- 20
#' x <- factor(1:K, ordered = TRUE, levels = c(1:K))
#'
#' # Valid example -------------------------------------------------
#' isOrderedFactor(x) # NULL
#'
#' # Examples throwing errors ---------------------------------------
#' try(isOrderedFactor(c(x, NA)))
#' try(isOrderedFactor(factor(1:K, levels = c(1:K))))
#' try(isOrderedFactor(factor(rep(1, K), levels = 1)))
isOrderedFactor <- function(x) {
  if (!is.ordered(x)) {
    stop("`x` must be an ordered factor. Instead, it is a ", class(x), ".")
  }

  if (any(is.na(x))) {
    stop("`x` must not contain NA values.")
  }

  if (length(levels(x)) < 2) {
    stop("The number of levels of `x` must be larger than 1. Instead, it is equal to", length(levels(x)), ".")
  }
}
