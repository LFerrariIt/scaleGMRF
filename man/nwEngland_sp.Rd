% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nwEngland_sp.R
\docType{data}
\name{nwEngland_sp}
\alias{nwEngland_sp}
\title{North West England districts}
\format{
An object of class \code{SpatialPolygons} of length 24.
}
\source{
https://CRAN.R-project.org/package=spBayesSurv
}
\usage{
nwEngland_sp
}
\description{
\code{nwEngland_sp} is a \code{SpatialPolygons} object representing the 24 districts in North West England found in the \code{INLA::Leuk} dataset of INLA (also, the \code{spBayesSurv::LeukSurv}) from the data studied in Henderson et al. (2002). "Modeling spatial variation in leukaemia survival data". JASA, 97(460), 965-972. The object has been built using internal files from the \code{spBayesSurv} package.
}
\details{
The \code{INLA::Leuk} dataset is used to illustrate the usage of the \code{scaleGMRF} package in the vignette \code{vignette("Leuk_case_study",package="scaleGMRF")}. \code{nwEngland_sp} is used to map the results for the districts.
}
\keyword{datasets}
