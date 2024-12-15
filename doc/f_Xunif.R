## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  message=FALSE,
  collapse = TRUE,
  fig.width =7,
  fig.height=12,
  comment = "#>"
)

## ----echo=FALSE---------------------------------------------------------------
library(scaleGMRF)

## -----------------------------------------------------------------------------
K <- 20
x <- factor(1:K, ordered = TRUE, levels = c(1:K))
result <- f_Xunif(x, model = "iid",fixed=T, plot_check = TRUE) 

## -----------------------------------------------------------------------------
result <- f_Xunif(x, model = "iid",fixed=F, plot_check = TRUE) 

## -----------------------------------------------------------------------------
result <- f_Xunif(x, model = "rw1", plot_check = TRUE) 

## -----------------------------------------------------------------------------
result <- f_Xunif(x, model = "rw2", plot_check = TRUE) 

## -----------------------------------------------------------------------------
data("nwEngland_adj_mat")
K <- nrow(nwEngland_adj_mat)
x <- factor(1:K, ordered = TRUE, levels = c(1:K))
result <- f_Xunif(x, model = "besag", adj_mat=nwEngland_adj_mat, plot_check = TRUE) 

## -----------------------------------------------------------------------------
x <- runif(1000, min = 2, max = 5)
result <- f_Xunif(x, model = "linear", plot_check = TRUE) 

## -----------------------------------------------------------------------------
K <- 20 
result <- f_Xunif(x, model = "pspline1", K = K, plot_check = TRUE) 

## -----------------------------------------------------------------------------
K <- 20 
result <- f_Xunif(x, model = "pspline2", K = K, plot_check = TRUE) 

## -----------------------------------------------------------------------------
x <- as.matrix(expand.grid(
  seq(0, 1, length.out = 100),
  seq(0, 1, length.out = 100)
))
result <- f_Xunif(x,model="pspline_2D", K = c(5, 5),plot_check = TRUE)

