% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspline.R
\name{bspline}
\alias{bspline}
\title{Create a basis of cubic uniform B-Spline functions}
\usage{
bspline(x, K, m = 0, M = 1)
}
\arguments{
\item{x}{A numeric vector, with values contained between 0 and 1.}

\item{K}{A positive integer larger than 3, specifying the number of basis functions.}

\item{m, M}{A pair of numbers indicating respectively the lower and upper boundaries of the support. By default, 0 and 1.}
}
\value{
A matrix with \code{K} columns and a number of rows equal to the length of \code{x}.
}
\description{
\code{bspline()} creates a basis matrix with \code{K} columns, each of them evaluating the numeric vector \code{x} at a cubic B-Spline function defined on equidistant nodes on the \verb{m,M} interval.
}
\examples{
x <- seq(0, 1, length.out = 100)
bspline(x, K = 10)
}
