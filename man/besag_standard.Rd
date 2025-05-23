% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/besag.R
\name{besag_standard}
\alias{besag_standard}
\title{Standardized Besag (ICAR) effects}
\usage{
besag_standard(x, adj_mat)
}
\arguments{
\item{x}{An ordered factor, with more than 1 level.}

\item{adj_mat}{Matrix representing the adjacency matrix of the cells of a lattice.}
}
\value{
A list of 6 elements, containing:
\itemize{
\item \code{Q}: precision matrix
\item \code{D}: basis matrix evaluated at \code{x}
\item \code{C}: a positive number, representing the appropriate scaling constant
\item \code{null_space}: a matrix, representing the null space of the precision matrix
\item \code{X_distribution}: a numeric vector (or matrix for \code{pspline2D}) of values sampled from the Uniform distribution of X
\item \code{basis_distribution}: basis matrix evaluated at \code{X_distribution}
}
}
\description{
#'@description
\code{besag_standard()} provides a list of elements to build a standardized Besag effect.
}
\examples{
# Example with the North West England districts ------------
data("nwEngland_adj_mat")
K <- nrow(nwEngland_adj_mat)
x <- factor(1:K, ordered = TRUE, levels = c(1:K))
besag_standard(x, adj_mat = nwEngland_adj_mat)
}
