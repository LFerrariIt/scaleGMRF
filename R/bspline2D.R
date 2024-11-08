
bspline2D <- function(x1, x2, N_basis1, N_basis2) {

  if (length(x1)!=length(x2)) {
    stop(paste0("The vectors `x1` and `x2` must have the same length. Instead,",length(x1),"!=",length(x2),"."))
  }

  B_1 <- bspline(x1, N_basis1)
  B_2 <- bspline(x2, N_basis2)

  B <- matrix(NA, nrow = length(x1), ncol = N_basis1 * N_basis2)

  for (i in 1:length(x1)) {
    B[i, ] <- kronecker(B_2[i, ], B_1[i, ])
  }
  return(B)
}
