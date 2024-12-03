#                   SURVIVAL ANALYSIS CASE STUDY

# Install the package scaleGMRF
rm(list=ls())
if(!require("scaleGMRF")) {
  library(devtools)
  install_github("LFerrariIt/scaleGMRF")
}

# Library loading -------------------------------------------
library(R2BayesX)
library(fastDummies)
library(spam)
library(ggpubr)
library(scaleGMRF)
library(INLA)

# This option is necessary for the use of inla.jp.define()
inla.setOption(num.threads = "1")

# Step 1: creation of the transformed dataset---------------
#The dataset is available in INLA as Leuk
data(Leuk)

K_T <- 27 # number of intervals for time (about half a year)

# This INLA function returns the transformed data
cph.leuk <- inla.coxph(
  inla.surv(time, cens) ~ 0 ,
  data = data.frame(Leuk),
  control.hazard = list(n.intervals = K_T))

#Transformed data with appropriate column names
transf_data <- cph.leuk$data[,c(1:3,5,12:16)]
colnames(transf_data) <- c("y","E","patient","time",
                           colnames(transf_data[,5:9]))

#Step 2: define the effects of the model ----

## Age, Wbc, Tpi --------------------------------------------
# Mean and variance assuming Uniformity on the empirical range
E_AGE <- (min(transf_data$age)+max(transf_data$age))/2
E_WBC <- (min(transf_data$wbc)+max(transf_data$wbc))/2
E_TPI <- (min(transf_data$tpi)+max(transf_data$tpi))/2

C_t_AGE <- ((max(transf_data$age)-min(transf_data$age))^2)/12
C_t_WBC <- ((max(transf_data$wbc)-min(transf_data$wbc))^2)/12
C_t_TPI <- ((max(transf_data$tpi)-min(transf_data$tpi))^2)/12

# Basis matrices of the linear effect under 0-mean constraint
D_t_AGE <- as.matrix(transf_data$age-E_AGE)
D_t_WBC <- as.matrix(transf_data$wbc-E_WBC)
D_t_TPI <- as.matrix(transf_data$tpi-E_TPI)

# Set up the number K_X of basis functions for the P-Spline effects
K_X <- 50

# Basis matrices of the non-linear effect
D_r_AGE <-
  pspline_standard(transf_data$age,K = K_X,order = 2)$basis
D_r_WBC <-
  pspline_standard(transf_data$wbc,K = K_X,order = 2)$basis
D_r_TPI <-
  pspline_standard(transf_data$tpi,K = K_X,order = 2)$basis

# Other elements for the residual effect
Pspline_elements <-
  pspline_standard(transf_data$tpi,K = K_X,order = 2)
Q_r <- Pspline_elements$precisio/Pspline_elements$scaling_constant #the function already returns the scaled precision matrix: here, we save the unscaled precision and apply scaling manually
S_r <- Pspline_elements$null_space
C_r <- Pspline_elements$scaling_constant
# check that the null space is correct:
Q_r %*% S_r

## Space, time, and sex -----------------------------------

# Basis matrices
D_S <- as.matrix(
  dummy_cols(transf_data$district,remove_selected_columns = T))
D_T <- as.matrix(
  dummy_cols(transf_data$time,remove_selected_columns = T))
D_SEX <- as.matrix(
  dummy_cols(transf_data$sex,remove_selected_columns = T))

# Precision matrices
#map of the districts in North East England
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd",
                                 package = "spBayesSurv"))
Q_S <- bnd2gra(nwengland)[(1:24),(1:24)]
Q_T <- as.matrix(precmat.RW1(K_T))
Q_SEX <- gen_inv(matrix(c(0.5,-0.5,-0.5,0.5),ncol=2),rank_def = 1)

# Scaling constants
C_S <- mean(diag(gen_inv(Q_S,rank_def=1)))
C_T <- mean(diag(gen_inv(Q_T,rank_def=1)))
C_SEX <- 1/2

# Step 3: creation of the data stack for INLA ---------------

data_stack <- inla.stack(
  data = list(Y=as.matrix(transf_data$y,ncol=1)),
  # Basis matrices
  A = list(1,
           D_t_AGE,
           D_t_WBC,
           D_t_TPI,
           D_r_AGE,
           D_r_WBC,
           D_r_TPI,
           D_SEX,
           D_S,
           D_T
  ),
  effects=list(data.frame(
    intercept=rep(1,nrow(transf_data))
  ),
  # Coefficients
  list(beta_AGE=1),
  list(beta_WBC=1),
  list(beta_TPI=1),
  list(u_AGE=1:K_X),
  list(u_WBC=1:K_X),
  list(u_TPI=1:K_X),
  list(u_SEX=1:2),
  list(u_S=1:24),
  list(u_T=1:K_T)
  ),
  remove.unused = F)

# Step 4: design of the VP prior ----------------------

VP_log_joint_prior <- function(theta, theta.desc = NULL) {

  # Hyperparameter of the symmetric Dirichlet
  alpha_omega <- 1

  #theta contains the log precisions of the parameters
  t_1  <- theta[1]
  t_2  <- theta[2]
  t_3  <- theta[3]
  t_4  <- theta[4]
  t_5  <- theta[5]
  t_6  <- theta[6]
  t_7  <- theta[7]
  t_8  <- theta[8]
  t_9  <- theta[9]

  ##Step A: compute the new parameters as expressions ---------------------------

  # Total Variance
  total_Var <- expression(
        exp(-t_1)+exp(-t_2)+exp(-t_3)+
        exp(-t_4)+exp(-t_5)+exp(-t_6)+
        exp(-t_7)+exp(-t_8)+exp(-t_9))

  #proportion 1
  omega_1 <- expression(exp(-t_1)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 2
  omega_2 <- expression(exp(-t_2)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 3
  omega_3<- expression(exp(-t_3)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 4
  omega_4 <- expression(exp(-t_4)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 5
  omega_5 <- expression(exp(-t_5)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 6
  omega_6 <- expression(exp(-t_6)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 7
  omega_7 <- expression(exp(-t_7)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))

  #proportion 8
  omega_8 <- expression(exp(-t_8)/
                          (exp(-t_1)+exp(-t_2)+exp(-t_3)+
                             exp(-t_4)+exp(-t_5)+exp(-t_6)+
                             exp(-t_7)+exp(-t_8)+exp(-t_9)))


  ##Step B: compute the derivatives in a matrix ---------------------------------

  D_matrix <- matrix(data =
                       c(
                         #
                         eval(D(total_Var,"t_1")) ,
                         eval(D(omega_1,  "t_1")) ,
                         eval(D(omega_2,  "t_1")) ,
                         eval(D(omega_3,  "t_1")) ,
                         eval(D(omega_4,  "t_1")) ,
                         eval(D(omega_5,  "t_1")),
                         eval(D(omega_6,  "t_1")) ,
                         eval(D(omega_7,  "t_1")) ,
                         eval(D(omega_8,  "t_1")),
                         #
                         eval(D(total_Var,"t_2")) ,
                         eval(D(omega_1,  "t_2")) ,
                         eval(D(omega_2,  "t_2")) ,
                         eval(D(omega_3,  "t_2")) ,
                         eval(D(omega_4,  "t_2")) ,
                         eval(D(omega_5,  "t_2")),
                         eval(D(omega_6,  "t_2")) ,
                         eval(D(omega_7,  "t_2")) ,
                         eval(D(omega_8,  "t_2")),
                         #
                         eval(D(total_Var,"t_3")) ,
                         eval(D(omega_1,  "t_3")) ,
                         eval(D(omega_2,  "t_3")) ,
                         eval(D(omega_3,  "t_3")) ,
                         eval(D(omega_4,  "t_3")) ,
                         eval(D(omega_5,  "t_3")),
                         eval(D(omega_6,  "t_3")) ,
                         eval(D(omega_7,  "t_3")) ,
                         eval(D(omega_8,  "t_3")),
                         #
                         eval(D(total_Var,"t_4")) ,
                         eval(D(omega_1,  "t_4")) ,
                         eval(D(omega_2,  "t_4")) ,
                         eval(D(omega_3,  "t_4")) ,
                         eval(D(omega_4,  "t_4")) ,
                         eval(D(omega_5,  "t_4")),
                         eval(D(omega_6,  "t_4")) ,
                         eval(D(omega_7,  "t_4")) ,
                         eval(D(omega_8,  "t_4")),
                         #
                         eval(D(total_Var,"t_5")) ,
                         eval(D(omega_1,  "t_5")) ,
                         eval(D(omega_2,  "t_5")) ,
                         eval(D(omega_3,  "t_5")) ,
                         eval(D(omega_4,  "t_5")) ,
                         eval(D(omega_5,  "t_5")),
                         eval(D(omega_6,  "t_5")) ,
                         eval(D(omega_7,  "t_5")) ,
                         eval(D(omega_8,  "t_5")),
                         #
                         eval(D(total_Var,"t_6")) ,
                         eval(D(omega_1,  "t_6")) ,
                         eval(D(omega_2,  "t_6")) ,
                         eval(D(omega_3,  "t_6")) ,
                         eval(D(omega_4,  "t_6")) ,
                         eval(D(omega_5,  "t_6")),
                         eval(D(omega_6,  "t_6")) ,
                         eval(D(omega_7,  "t_6")) ,
                         eval(D(omega_8,  "t_6")),
                         #
                         eval(D(total_Var,"t_7")) ,
                         eval(D(omega_1,  "t_7")) ,
                         eval(D(omega_2,  "t_7")) ,
                         eval(D(omega_3,  "t_7")) ,
                         eval(D(omega_4,  "t_7")) ,
                         eval(D(omega_5,  "t_7")),
                         eval(D(omega_6,  "t_7")) ,
                         eval(D(omega_7,  "t_7")) ,
                         eval(D(omega_8,  "t_7")),
                         #
                         eval(D(total_Var,"t_8")) ,
                         eval(D(omega_1,  "t_8")) ,
                         eval(D(omega_2,  "t_8")) ,
                         eval(D(omega_3,  "t_8")) ,
                         eval(D(omega_4,  "t_8")) ,
                         eval(D(omega_5,  "t_8")),
                         eval(D(omega_6,  "t_8")) ,
                         eval(D(omega_7,  "t_8")) ,
                         eval(D(omega_8,  "t_8")),
                         #
                         eval(D(total_Var,"t_9")) ,
                         eval(D(omega_1,  "t_9")) ,
                         eval(D(omega_2,  "t_9")) ,
                         eval(D(omega_3,  "t_9")) ,
                         eval(D(omega_4,  "t_9")) ,
                         eval(D(omega_5,  "t_9")),
                         eval(D(omega_6,  "t_9")) ,
                         eval(D(omega_7,  "t_9")) ,
                         eval(D(omega_8,  "t_9"))
                       ),ncol=9)

  ##Step C: add log densities of necessary distributions --------------------------

  # Dirichlet distribution

  log_ddirichlet <- function (x, alpha)   {

    x <- c(x,1-sum(x))

    return(lgamma(length(x)*alpha) -
             length(x)*lgamma(alpha) +
             (alpha - 1)*sum(log(x)))
  }

  ##Step D: compute the prior for theta through change of variable formula ------

  log_prior <-
    -log(eval(total_Var)) + # Jeffreys on V
    log_ddirichlet(x=c(   # Symmetric dirichlet for omega
      eval(omega_1),
      eval(omega_2),
      eval(omega_3),
      eval(omega_4),
      eval(omega_5),
      eval(omega_6),
      eval(omega_7),
      eval(omega_8)),
      alpha=alpha_omega)+
    log(abs(det(D_matrix))) #jacobian determinant

  return(log_prior)
}

# Step 5: fit the model with unscaled effects ----------------

unscaled_model <- inla(
  Y ~ -1+intercept+
    # Linear effects
    f(beta_AGE, model="generic0",Cmatrix = as.matrix(1))+
    f(beta_WBC, model="generic0",Cmatrix = as.matrix(1))+
    f(beta_TPI, model="generic0",Cmatrix = as.matrix(1))+
    # Non-linear effects
    f(u_AGE,model="generic0",Cmatrix = Q_r, constr = F,
      extraconstr = list(A=t(S_r),e=matrix(c(0,0))))+
    f(u_WBC,model="generic0",Cmatrix = Q_r,constr = F,
      extraconstr = list(A=t(S_r),e=matrix(c(0,0))))+
    f(u_TPI,model="generic0",Cmatrix = Q_r,constr = F,
      extraconstr = list(A=t(S_r), e=matrix(c(0,0))))+
    # Cluster effect
    f(u_SEX,model="generic0",Cmatrix = Q_SEX,constr = T)+
    # Besag effect
    f(u_S,  model="generic0",Cmatrix = Q_S,constr = T)+
    # RW1 effect
    f(u_T,  model="generic0",Cmatrix = Q_T,constr = T),
  family = 'poisson',
  data = inla.stack.data(data_stack),
  E =  transf_data$E, verbose=T,
  # VP prior
  control.expert = list(jp=inla.jp.define(VP_log_joint_prior)),
  control.predictor = list(A = inla.stack.A(data_stack)))

# Step 5: fit the model with scaled effects ----------------

scaled_model <- inla(
  Y ~ -1+intercept+
    f(beta_AGE, model="generic0",Cmatrix = C_t_AGE*as.matrix(1))+
    f(beta_WBC, model="generic0",Cmatrix = C_t_WBC*as.matrix(1))+
    f(beta_TPI, model="generic0",Cmatrix = C_t_TPI*as.matrix(1))+
    f(u_AGE,model="generic0",Cmatrix = C_r*Q_r, constr = F,
      extraconstr = list(A=t(S_r),e=matrix(c(0,0))))+
    f(u_WBC,model="generic0",Cmatrix = C_r*Q_r,constr = F,
      extraconstr = list(A=t(S_r),e=matrix(c(0,0))))+
    f(u_TPI,model="generic0",Cmatrix = C_r*Q_r,constr = F,
      extraconstr = list(A=t(S_r), e=matrix(c(0,0))))+
    f(u_SEX,model="generic0",Cmatrix = C_SEX*Q_SEX,constr = T)+
    f(u_S,  model="generic0",Cmatrix = C_S*Q_S,constr = T)+
    f(u_T,  model="generic0",Cmatrix = C_T*Q_T,constr = T),
  family = 'poisson',
  data = inla.stack.data(data_stack),
  E =  transf_data$E,verbose=T,
  control.expert = list(jp=inla.jp.define(VP_log_joint_prior)),
  control.predictor = list(A = inla.stack.A(data_stack)))

# Step 6: compare the differences between the two models

covariates_plot <-
ggarrange(
  ggplot()+
  geom_step(aes(
    col = " Scaling",
    x = c(unique(cph.leuk$data$baseline.hazard.time),
          max(Leuk)) / 365.25,
    y = c(
      scaled_model$summary.random$u_T$mean,
      scaled_model$summary.random$u_T$mean[K_T]
    )
  )) +
    geom_step(aes(
      col = "No scaling",
      x = c(unique(cph.leuk$data$baseline.hazard.time), max(Leuk)) / 365.25,
      y = c(
        unscaled_model$summary.random$u_T$mean,
        unscaled_model$summary.random$u_T$mean[K_T]
      )
    )) +
    labs(x = "Time (years)", col = "", y = "Posterior mean") +
    theme_light()+
    scale_color_manual(values=c(" Scaling"="black","No scaling"="red")),
  ggplot() +
    geom_line(aes(
      col = " Scaling",
      x = transf_data$age,
      y = D_t_AGE %*% scaled_model$summary.random$beta_AGE$mean +
        D_r_AGE %*% scaled_model$summary.random$u_AGE$mean
    )) +
    geom_line(aes(
      col = "No scaling",
      x = transf_data$age,
      y = D_t_AGE %*% unscaled_model$summary.random$beta_AGE$mean +
        D_r_AGE %*% unscaled_model$summary.random$u_AGE$mean
    )) +
    labs(x = "Age", col = "", y = "Posterior mean") +
    theme_light()+
    scale_color_manual(values=c(" Scaling"="black","No scaling"="red")),
  ggplot() +
    geom_line(aes(
      col = " Scaling",
      x = transf_data$wbc,
      y = D_t_WBC %*% scaled_model$summary.random$beta_WBC$mean +
        D_r_WBC %*% scaled_model$summary.random$u_WBC$mean
    )) +
    geom_line(aes(
      col = "No scaling",
      x = transf_data$wbc,
      y = D_t_WBC %*% unscaled_model$summary.random$beta_WBC$mean +
        D_r_WBC %*% unscaled_model$summary.random$u_WBC$mean)) +
    labs(x = "White blood cell count", col = "", title="(a)",
         y = "Posterior mean") +
    scale_color_manual(values=c(" Scaling"="black","No scaling"="red"))+
    theme_light()+ theme(plot.title = element_text(hjust = 0.5)),
  ggplot() +
    geom_line(aes(
      col = " Scaling",
      x = transf_data$tpi,
      y = D_t_TPI %*% scaled_model$summary.random$beta_TPI$mean+
        D_r_TPI %*% scaled_model$summary.random$u_TPI$mean
    )) +
    geom_line(aes(
      col = "No scaling",
      x = transf_data$tpi,
      y = D_t_TPI%*%unscaled_model$summary.random$beta_TPI$mean+
        D_r_TPI%*%unscaled_model$summary.random$u_TPI$mean))+
    labs(x = "Townsend deprivation index", title="(b)",
         col = "", y = "Posterior mean") +
    scale_color_manual(values=c(" Scaling"="black","No scaling"="red"))+
    theme_light()+
    theme(plot.title = element_text(hjust = 0.5)),common.legend = T,legend = "bottom"
)

ggsave(filename = "covariates.png",covariates_plot,
       width=1600,height=1800,units="px",bg="white")

dst_map <- fortify(bnd2sp(nwengland))
dst_scaled <- scaled_model$summary.random$u_S$mean[as.numeric(districts_map$id)]
dst_unscaled <- unscaled_model$summary.random$u_S$mean[as.numeric(districts_map$id)]

districts_plot <-
ggarrange(
  ggplot(dst_map, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill=dst_scaled))+
    labs(title="Scaling",fill="")+
    scale_fill_gradient2(
      low="#2cba00",mid="gold",high="red",
      limits=range(c(dst_unscaled,dst_scaled)))+
    theme_minimal()+
    theme(legend.position = "bottom"),
  ggplot(dst_map, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill=dst_unscaled))+
    labs(title="No scaling",fill="")+
    scale_fill_gradient2(
      low="#2cba00",mid="gold",high="red",
      limits=range(c(dst_scaled,dst_scaled)))+
    theme_minimal()+
    theme(legend.position = "bottom"))

ggsave(filename = "districts.png",districts_plot,
       width=1600,height=1200,units="px",bg="white")
