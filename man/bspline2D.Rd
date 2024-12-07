% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline2D.R
\name{bspline2D}
\alias{bspline2D}
\title{Create a basis of cubic uniform 2D B-Spline functions}
\usage{
bspline2D(x, K1, K2)
}
\arguments{
\item{x}{A matrix of 2 columns, containing values between 0 and 1.}

\item{K1}{A positive integer larger than 3, specifying the number of basis functions on the first dimension.}

\item{K2}{A positive integer larger than 3, specifying the number of basis functions on the second dimension.}
}
\value{
A matrix with `K1`x`K2` columns and a number of rows equal to the length of `x1` and `x2`.
}
\description{
`bspline2D()` creates a basis matrix with `K1`x`K2` columns, each of them evaluating the numeric vector `x` at a 2-dimensional cubic B-Spline function, defined through the Kronecker product of 2 unidimensional B-Spline functions.
}
\examples{
x_grid <- as.matrix(expand.grid(
  seq(0, 1, length.out = 100),
  seq(0, 1, length.out = 100)
))

bspline2D(x_grid, K1 = 5, K2 = 5)

}
