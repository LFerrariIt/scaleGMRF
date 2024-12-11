#' Adjacency matrix for the North West England districts
#'
#' @description  `nwEngland_adj_mat` is the adjacency matrix for the 24 districts in North West England found in the `INLA::Leuk` dataset of INLA (also, the `spBayesSurv::LeukSurv`) from the data studied in Henderson et al. (2002). Modeling spatial variation in leukemia survival data. JASA, 97(460), 965-972. The object has been built using internal files from the `spBayesSurv` package.
#'
#' @format The matrix is a square, symmetric matrix made up by 24 rows and columns. Each row and column represent a different district. If a cell contains a 1, then the districts corresponding to the row and column of the cell are adjacent (i.e. neighbours). Otherwise, if the cell contains a 0, the corresponding districts do not share a border.
#'
#' @details The `INLA::Leuk` dataset is used to illustrate the usage of the `scaleGMRF` package in the vignette `vignette("Leuk_case_study",package="scaleGMRF")`. `nwEngland_sp` is used to create the precision matrix of a Besag/ICAR effect for the `district` areal covariate in `INLA::Leuk`.
#'
#' @source  https://CRAN.R-project.org/package=spBayesSurv
#'
"nwEngland_adj_mat"
