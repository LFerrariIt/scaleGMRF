# # Comparison

# K <- 25
# adj_mat <- diag(diag(as.matrix(precmat.IGMRFreglat(5,5))))-as.matrix(precmat.IGMRFreglat(5,5))
# x <- factor(1:K, ordered = TRUE, levels = c(1:K))
# f_iid_f <-standardize_GMRF(Q=diag(K),fixed=T)
# f_iid_r <-standardize_GMRF(Q=diag(K),fixed=F)
# f_rw1   <-standardize_GMRF(Q=as.matrix(precmat.RW1(K)),rank_def=1)
# f_rw2   <-standardize_GMRF(Q=as.matrix(precmat.RW2(K)),rank_def=2)
# f_besag <-standardize_GMRF(Q=as.matrix(precmat.IGMRFreglat(5,5)),rank_def=1)
#

# f2_iid_f <-standardize_X_unif(x,model="iid",fixed=T)
# f2_iid_r <-standardize_X_unif(x,model="iid",fixed=F)
# f2_rw1   <-standardize_X_unif(x,model="rw1")
# f2_rw2   <-standardize_X_unif(x,model="rw2")
# f2_besag <-standardize_X_unif(x,model="besag",adj_mat = adj_mat)
#
# x <- runif(100,4,10)
# f_linear <-standardize_GMRF(Q=diag(1),D = matrix(x))
# f2_linear <- standardize_X_unif(x,model="linear")
#
# #------------------------------------------------------------------------------
#
# K <- 25
# # Comparison
# adj_mat <- diag(diag(as.matrix(precmat.IGMRFreglat(5,5))))-as.matrix(precmat.IGMRFreglat(5,5))
# x <- factor(1:K, ordered = TRUE, levels = c(1:K))
#
#
# standardize_GMRF(plot_check = T,Q=diag(K),fixed=T)
# standardize_GMRF(plot_check = T,Q=diag(K),fixed=F)
# standardize_GMRF(plot_check = T,fixed=F,Q=as.matrix(precmat.RW1(K)),rank_def=1)
# standardize_GMRF(plot_check = T,fixed=F,Q=as.matrix(precmat.RW2(K)),rank_def=2)
# standardize_GMRF(plot_check = T,fixed=F,Q=as.matrix(precmat.IGMRFreglat(5,5)),rank_def=1)
#
# standardize_X_unif(plot_check = T,x,model="iid",fixed=T)
# standardize_X_unif(plot_check = T,x,model="iid",fixed=F)
# standardize_X_unif(plot_check = T,x,model="rw1",fixed=F)
# standardize_X_unif(plot_check = T,x,model="rw2",fixed=F)
# standardize_X_unif(plot_check = T,x,model="besag",adj_mat = adj_mat,fixed=F)
#
# x <- sort(runif(100,4,10))
# standardize_GMRF(plot_check = T,fixed=F,Q=diag(1),D = matrix(x))
# standardize_X_unif(plot_check = T,x,model="linear",fixed=F)


