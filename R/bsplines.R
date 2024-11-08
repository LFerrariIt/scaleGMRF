bspline <- function(x, N_basis=20) {

  dx <- 1 / (N_basis-3)
  knots <- seq(- 3* dx, 1 + 3* dx, by = dx)
  B <- splines::spline.des(knots, x, ord=4, 0 * x, outer.ok = F)$design
  return(B)

}

bspline_2D <- function(y1,y2, N_basis1,N_basis2) {

  B_1 <- bspline(y1,N_basis1)
  B_2 <- bspline(y2,N_basis2)

  B <- matrix(NA,nrow=length(y1),ncol=N_basis1*N_basis2)

  for (i in 1:length(y1)) {
    B[i,] <- kronecker(B_2[i,],B_1[i,])
  }
  return(B)
}
