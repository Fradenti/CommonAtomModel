
# Loading Libraries -------------------------------------------------------
library(plyr); 
library(Rcpp); 
library(RcppArmadillo)
library(reshape); 
library(pryr); 
library(knitr); 
library(gridExtra)
library(tidyverse); 
library(MCMCpack); 
library(scales)
#
#Rcpp::sourceCpp("Functions/CAM/CAM.cpp")
#
log_post_beta <- function(beta, Lk, n_j, K.star, ab,bb){
  K.star*lgamma(beta) -  sum(lgamma(beta + n_j)) + 
    ( sum(Lk+ ab-1) ) * log(beta) - bb*beta
}


log_post_LOGbeta <- function(LOGbeta, Lk, n_j, K.star, ab,bb){
  log_post_beta(exp(LOGbeta), Lk = Lk, n_j = n_j, K.star = K.star, ab = ab,bb = bb) +
    LOGbeta
}

MH.beta <- function(beta, Lk, mlk, K.star,ab,bb,sigma.prop){
  
  old.log.beta <- log(beta)
  new.log.beta <- rnorm(1,old.log.beta,sigma.prop)
  
  
  alpha <- log_post_LOGbeta(LOGbeta = new.log.beta, Lk = Lk, n_j = mlk, K.star = K.star, ab = ab,bb = bb)-
    log_post_LOGbeta(LOGbeta = old.log.beta, Lk = Lk, n_j = mlk, K.star = K.star, ab = ab,bb = bb)
  
  if(runif(1)<exp(alpha)){
    beta <- exp(new.log.beta)
  }else{
    beta <- exp(old.log.beta)
  }
  
  return(beta)
}
# Auxiliary functions -----------------------------------------------------
reset <- function(x){
  z  <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}
#######################################################################################################
CAM <- function(y_obser,                         # Observations, organized into
                y_group,                         # Groups (numeric vector 1- First group-->J- J-th group)
                K0=10, L0=20,                    # Starting number of groups
                prior,                           # List of hyperparameters
                nsim=100, burn_in = 100,         # Number of Simulations and Burn in
                thinning = 1,                    # Thinning interval 
                verbose=1,                 # Print the iteration number every ... steps
                fixedAB=T,                       # Do you want to keep the concentration parameters of the two DPs fixed?
                cheap=F,                         # If T, the output will contain only the membership labels
                kappa=0.5,                       # Slice sampler parameter
                seed=NA,                         # Specify if you want to set the random seed
                sigma.prop.beta=.25) {                       
  
  if (!is.na(seed)) {
    set.seed(seed)
  }
  
  J       <- length(unique(y_group)) # number of groups
  nj      <- table(y_group)          # number of observations inside each group
  N       <- length(y_obser)
  ############################### hyperparameters extraction
  m0 <- prior$m0
  k0 <- prior$k0
  
  a0 <- prior$a0
  b0 <- prior$b0
  
  a_alpha <- prior$a_alpha
  b_alpha <- prior$b_alpha
  a_beta  <- prior$a_beta
  b_beta  <- prior$b_beta
  
  
  #####################################################################################
  # Initialization of the Chain:
  # Containers ##########################################################################################################################
  theta                  <- matrix(NA, L0, 2)
  cij                   <- numeric(N)
  zj                    <- numeric(J)
  BETA_DP               <- numeric(nsim); 
  ALPHA_DP              <- numeric(nsim)
  pi_star_k             <- vector("list",length = nsim)
  OMEGA                 <- vector("list",length = nsim)
  Csi_ij                <- vector("list",length = nsim)
  Z_j                   <- vector("list",length = nsim)
  THETA                 <- vector("list",length = nsim)
  #####################################################################################
  # Warmstart observational/distributional clusters
  km.out    <- kmeans(y_obser, L0)
  cij       <- km.out$cluster
  theta[, 1]  <- tapply(y_obser, cij, mean)
  theta[, 2]  <- tapply(y_obser, cij, sd)
  if (K0 >= J) {
    zj        <- 1:J
  } else{
    km.outZ   <- kmeans(tapply(y_obser, y_group, mean), K0)
    zj        <- km.outZ$cluster
  }
  zj.pg    <- rep(zj, nj)
  #####################################################################################
  g_slice <- function(j) {
    (1 - kappa) * kappa ^ (j - 1)
  }
  #####################################################################################
  # hyperparameters gamma distribution of alpha and beta of DP
  alpha   <- rgamma(1, a_alpha, b_alpha)
  beta    <- rgamma(1, a_beta,  b_beta)
  #####################################################################################
  
  if (verbose) {
    cat("MCMC progress:\n")
    flush.console()
    pbar <- txtProgressBar(min = 1,
                           max = nsim*thinning + burn_in,
                           style = 3)
    on.exit(close(pbar))
    ipbar <- 0
  }
  # Main loop
  for(sim in 2:(nsim*thinning + burn_in)){
    ################################################################
    # Stochastic Truncation - Distributional
    zjr <- reset(zj)
    Uj   <- runif(J, 0, g_slice(zj))
    L.z  <- 1 + floor((log(Uj) - log(1 - kappa)) / log(kappa))
    J.z  <- max(zj)
    NN.z <- max(c(L.z, J.z))
    xi.z <- g_slice(1:NN.z)
    ################################################################
    # Stochastic Truncation - Observational
    cijr <- reset(cij)
    Uij     <- runif(N, 0, g_slice(cij))
    L.c     <- 1 + floor((log(Uij) - log(1 - kappa)) / log(kappa))
    J.c     <- max(cij)
    NN.c    <- max(c(L.c, J.c))
    xi.c    <- g_slice(1:NN.c)
    ################################################################
    # 
    if (NN.c > J.c){
      theta <- rbind(theta, matrix(0, (NN.c - J.c), 2))
    }
    ################################################################
    # Update Distributional Weights
    v.z <- pi.z  <- numeric(NN.z)
    v.omega <- omega <- matrix(0, NN.c, NN.z)
    v.z  <-
      c(Update_Distributional_Sticks(zj = zj, NN_z = NN.z, alpha = alpha))
    ################################################################
    if (length(unique(zj)) > 1) {
      pi.z <- SB_given_u2(v.z)
    } else{
      pi.z <- v.z
    }
    ################################################################   
    # Update observational Weights
    omega <- Update_omega(
      cij = cij,
      zj_pg = zj.pg,
      NN_c = NN.c,
      NN_z = NN.z,
      beta = beta
    )
    ################################################################
    # Update Distributional Labels
    zj<- c(
      Update_Zj_v2(
        Uj = Uj,
        xi_z = xi.z,
        xi_c = xi.c,
        pi_z = pi.z,
        cij = cij,
        omega = omega,
        y_group = y_group,
        NN_z = NN.z,
        J = J
      )
    )
    ################################################################
    zj.pg <- rep(zj, nj)
    ################################################################
    cij <- Update_Cij(y_obser = y_obser,
                      Uij = Uij,
                      xi_c = xi.c,
                      omega = omega,
                      zj_pg = zj.pg,
                      theta = theta,
                      N = N,
                      NN_c = NN.c)
    theta <-
      tslp  <-
      Update_theta(
        y_obser = y_obser,
        cij = cij,
        a0 = a0,
        b0 = b0,
        k0 = k0,
        m0 = m0,
        NN_c = NN.c,
        J = J
      )
    
    ################################################################
    # Update alpha distributional DP - Escobar and West 1995
    k      <- length(unique(zj))
    logeta <- log(rbeta(1, alpha + 1, J))
    Q      <- (a_alpha + k - 1) / (J * (b_alpha - logeta))
    pi_eta <- Q / (1 + Q)
    u <- runif(1)
    if (u < pi_eta) {
      alpha  <- rgamma(n = 1,
                       shape = a_alpha + k,
                       rate = b_alpha - logeta)
    } else{
      alpha  <- rgamma(n = 1,
                       shape = a_alpha + k - 1,
                       rate = b_alpha - logeta)
    }
    ################################################################
    ################################################################
    Kz    <-  length(unique(zj.pg))
    n_j   <-  table(zj.pg)
    T_j   <-  tapply(cij, zj.pg, function(x) length(unique(x)))
    
    beta  <- MH.beta(beta = beta, Lk = T_j, mlk = n_j, K.star = Kz, 
                     ab = a_beta, bb=b_beta,sigma.prop = sigma.prop.beta)
    # ################################################################
    if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
      rr                         <- floor((sim - burn_in)/thinning)
      pi_star_k[[rr]]            <- pi.z
      OMEGA[[rr]]        <- omega
      Csi_ij[[rr]]               <- cij
      Z_j[[rr]]                  <- zj
      THETA[[rr]]                <- theta
      ALPHA_DP[rr]               <- alpha
      BETA_DP[rr]                <- beta
    }
    
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
    
  }
  if(cheap){
    out <- list(Z_j=Z_j, Csi_ij=Csi_ij,
                A_DP=ALPHA_DP, B_DP=BETA_DP,
                y_obser=y_obser, y_group=y_group, nsim=nsim, prior = prior)
  }else{
    out <- list(Z_j=Z_j, Csi_ij=Csi_ij, pi_star_k=pi_star_k, THETA=THETA, 
                A_DP=ALPHA_DP, B_DP=BETA_DP,
                OMEGA=OMEGA,  y_obser=y_obser, y_group=y_group, nsim=nsim, prior = prior)
    
    
  }
  class(out) <- "CAM_SLICE" #Code-name for our model
  return(out)
}
