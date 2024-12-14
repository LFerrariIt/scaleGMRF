#                   SURVIVAL ANALYSIS CASE STUDY

# Install the package scaleGMRF
rm(list=ls())
if(!require("scaleGMRF")) {
  library(devtools)
  install_github("LFerrariIt/scaleGMRF")
}

# Library loading -------------------------------------------
library(scaleGMRF)
library(fastDummies)
library(spam)
library(ggpubr)
library(INLA)

# This option is necessary for the use of inla.jp.define()
inla.setOption(num.threads = "1")

# Step 1: creation of the transformed dataset---------------
#The dataset is available in INLA as Leuk
data(Leuk)
data(nwEngland_sp)
data(nwEngland_adj_mat)

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

transf_data$sex <- factor(transf_data$sex,ordered = T)
transf_data$district <- factor(transf_data$district,ordered = T)
transf_data$time <- factor(transf_data$time,ordered = T)

#Step 2: define the effects of the model ----

K_X <- 50
f_t_AGE <- f_Xunif(transf_data$age,model="linear",scale_Q=F)
f_t_WBC <- f_Xunif(transf_data$wbc,model="linear",scale_Q=F)
f_t_TPI <- f_Xunif(transf_data$tpi,model="linear",scale_Q=F)
f_n_AGE <- f_Xunif(transf_data$age,model="pspline2",K=K_X,scale_Q=F)
f_n_WBC <- f_Xunif(transf_data$wbc,model="pspline2",K=K_X,scale_Q=F)
f_n_TPI <- f_Xunif(transf_data$tpi,model="pspline2",K=K_X,scale_Q=F)
f_SEX   <- f_Xunif(transf_data$sex,model="iid",scale_Q=F)
f_DST   <- f_Xunif(transf_data$district,model="besag",scale_Q=F,
                   adj_mat = nwEngland_adj_mat)
f_TEM   <- f_Xunif(transf_data$time,transf_data,model="rw1",scale_Q=F)

f_scld_t_AGE <- f_Xunif(transf_data$age,model="linear")
f_scld_t_WBC <- f_Xunif(transf_data$wbc,model="linear")
f_scld_t_TPI <- f_Xunif(transf_data$tpi,model="linear")
f_scld_n_AGE <- f_Xunif(transf_data$age,model="pspline2",K=K_X)
f_scld_n_WBC <- f_Xunif(transf_data$wbc,model="pspline2",K=K_X)
f_scld_n_TPI <- f_Xunif(transf_data$tpi,model="pspline2",K=K_X)
f_scld_SEX   <- f_Xunif(transf_data$sex,model="iid")
f_scld_DST   <- f_Xunif(transf_data$district,model="besag",
                   adj_mat = nwEngland_adj_mat)
f_scld_TEM   <- f_Xunif(transf_data$time,transf_data,model="rw1")

# Step 3: creation of the data stack for INLA ---------------

data_stack <- inla.stack(
  data = list(Y=as.matrix(transf_data$y,ncol=1)),
  # Basis matrices
  A = list(1,
           f_t_AGE$basis,
           f_t_WBC$basis,
           f_t_TPI$basis,
           f_n_AGE$basis,
           f_n_WBC$basis,
           f_n_TPI$basis,
           f_SEX$basis,
           f_DST$basis,
           f_TEM$basis
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
  list(u_DST=1:24),
  list(u_TEM=1:K_T)
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
    f(beta_AGE, model="generic0",Cmatrix = f_t_AGE$precision)+
    f(beta_WBC, model="generic0",Cmatrix = f_t_WBC$precision)+
    f(beta_TPI, model="generic0",Cmatrix = f_t_TPI$precision)+
    # Non-linear effects
    f(u_AGE,model="generic0",Cmatrix = f_n_AGE$precision, constr = F,
      extraconstr = list(A=t(f_n_AGE$null_space),e=matrix(c(0,0))))+
    f(u_WBC,model="generic0",Cmatrix = f_n_WBC$precision,constr = F,
      extraconstr = list(A=t(f_n_WBC$null_space),e=matrix(c(0,0))))+
    f(u_TPI,model="generic0",Cmatrix = f_n_TPI$precision,constr = F,
      extraconstr = list(A=t(f_n_TPI$null_space), e=matrix(c(0,0))))+
    # Cluster effect
    f(u_SEX,model="generic0",Cmatrix =  f_SEX$precision,constr = T)+
    # Besag effect
    f(u_DST,model="generic0",Cmatrix = f_DST$precision,constr = T)+
    # RW1 effect
    f(u_TEM,model="generic0",Cmatrix = f_TEM$precision,constr = T),
  family = 'poisson',
  data = inla.stack.data(data_stack),
  E =  transf_data$E, verbose=T,
  # VP prior
  control.expert = list(jp=inla.jp.define(VP_log_joint_prior)),
  control.predictor = list(A = inla.stack.A(data_stack)))

# Step 5: fit the model with scaled effects ----------------

scaled_model <- inla(
  Y ~ -1+intercept+
    # Linear effects
    f(beta_AGE, model="generic0",Cmatrix = f_scld_t_AGE$precision)+
    f(beta_WBC, model="generic0",Cmatrix = f_scld_t_WBC$precision)+
    f(beta_TPI, model="generic0",Cmatrix = f_scld_t_TPI$precision)+
    # Non-linear effects
    f(u_AGE,model="generic0",Cmatrix = f_scld_n_AGE$precision, constr = F,
      extraconstr = list(A=t(f_scld_n_AGE$null_space),e=matrix(c(0,0))))+
    f(u_WBC,model="generic0",Cmatrix = f_scld_n_WBC$precision,constr = F,
      extraconstr = list(A=t(f_scld_n_WBC$null_space),e=matrix(c(0,0))))+
    f(u_TPI,model="generic0",Cmatrix = f_scld_n_TPI$precision,constr = F,
      extraconstr = list(A=t(f_scld_n_TPI$null_space), e=matrix(c(0,0))))+
    # Cluster effect
    f(u_SEX,model="generic0",Cmatrix = f_scld_SEX$precision,constr = T)+
    # Besag effect
    f(u_DST,model="generic0",Cmatrix = f_scld_DST$precision,constr = T)+
    # RW1 effect
    f(u_TEM,model="generic0",Cmatrix = f_scld_TEM$precision,constr = T),
  family = 'poisson',
  data = inla.stack.data(data_stack),
  E =  transf_data$E, verbose=T,
  # VP prior
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
          scaled_model$summary.random$u_TEM$mean,
          scaled_model$summary.random$u_TEM$mean[K_T]
        )
      )) +
      geom_step(aes(
        col = "No scaling",
        x = c(unique(cph.leuk$data$baseline.hazard.time), max(Leuk)) / 365.25,
        y = c(
          unscaled_model$summary.random$u_TEM$mean,
          unscaled_model$summary.random$u_TEM$mean[K_T]
        )
      )) +
      labs(x = "Time (years)", col = "", y = "Posterior mean") +
      theme_light()+
      scale_color_manual(values=c(" Scaling"="black","No scaling"="red")),
    ggplot() +
      geom_line(aes(
        col = " Scaling",
        x = transf_data$age,
        y = f_t_AGE$basis %*% scaled_model$summary.random$beta_AGE$mean +
          f_n_AGE$basis %*% scaled_model$summary.random$u_AGE$mean
      )) +
      geom_line(aes(
        col = "No scaling",
        x = transf_data$age,
        y = f_t_AGE$basis %*% unscaled_model$summary.random$beta_AGE$mean +
          f_n_AGE$basis %*% unscaled_model$summary.random$u_AGE$mean
      )) +
      labs(x = "Age", col = "", y = "Posterior mean") +
      theme_light()+
      scale_color_manual(values=c(" Scaling"="black","No scaling"="red")),
    ggplot() +
      geom_line(aes(
        col = " Scaling",
        x = transf_data$wbc,
        y = f_t_WBC$basis %*% scaled_model$summary.random$beta_WBC$mean +
          f_n_WBC$basis %*% scaled_model$summary.random$u_WBC$mean
      )) +
      geom_line(aes(
        col = "No scaling",
        x = transf_data$wbc,
        y = f_t_WBC$basis %*% unscaled_model$summary.random$beta_WBC$mean +
          f_n_WBC$basis %*% unscaled_model$summary.random$u_WBC$mean)) +
      labs(x = "White blood cell count", col = "", title="(a)",
           y = "Posterior mean") +
      scale_color_manual(values=c(" Scaling"="black","No scaling"="red"))+
      theme_light()+ theme(plot.title = element_text(hjust = 0.5)),
    ggplot() +
      geom_line(aes(
        col = " Scaling",
        x = transf_data$tpi,
        y = f_t_TPI$basis %*% scaled_model$summary.random$beta_TPI$mean+
          f_n_TPI$basis %*% scaled_model$summary.random$u_TPI$mean
      )) +
      geom_line(aes(
        col = "No scaling",
        x = transf_data$tpi,
        y = f_t_TPI$basis%*%unscaled_model$summary.random$beta_TPI$mean+
          f_n_TPI$basis%*%unscaled_model$summary.random$u_TPI$mean))+
      labs(x = "Townsend deprivation index", title="(b)",
           col = "", y = "Posterior mean") +
      scale_color_manual(values=c(" Scaling"="black","No scaling"="red"))+
      theme_light()+
      theme(plot.title = element_text(hjust = 0.5)),common.legend = T,legend = "bottom"
  )

ggsave(filename = "covariates.png",covariates_plot,
       width=1600,height=1800,units="px",bg="white")

dst_map <- fortify(nwEngland_sp)
dst_scaled <- scaled_model$summary.random$u_DST$mean[as.numeric(dst_map$id)]
dst_unscaled <- unscaled_model$summary.random$u_DST$mean[as.numeric(dst_map$id)]

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
