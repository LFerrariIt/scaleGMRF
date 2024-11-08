bspline <- function(x, N_basis = 20) {

  if (!is.numeric(x)) {
    stop(paste0("`x` must be a numeric vector."))
  }

  if (any(x < 0) | any(x > 1) | any(is.na(x))) {
    stop(paste0("The vector `x` must be normalized to lie between 0 and 1 and not contain NA values."))
  }

  dx <- 1 / (N_basis - 3)
  knots <- seq(-3 * dx, 1 + 3 * dx, by = dx)
  B <- splines::spline.des(knots, x, ord = 4, 0 * x, outer.ok = F)$design
  return(B)
}

