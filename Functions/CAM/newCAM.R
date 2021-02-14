
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



posterior.densities <- function(m1,howfarback, normalize=T, xx){
  L    <- list()
  J    <- length(m1$Z_j[1,])
  nsim <- length(m1$A_DP)
  for(h in 1:J){
    mix <- matrix(NA,howfarback,length(xx))
    for(i in 1:howfarback){
      w <- m1$OMEGA[[nsim-i]][,m1$Z_j[nsim-i,]][,h]
      
      if(normalize){
        w <- w/sum(w)
      }
      
      dens <- t(apply(m1$THETA[[nsim-i]],1,
                      function(z) dnorm(xx,z[1],sqrt(z[2])) ))
      
      mixcomp <- w*dens
      mix[i,] <-   colSums(mixcomp)
    }
    L[[h]] <- mix
    cat(h)
  }
  return(L)
}

g_slice <- function(j, kappa) {
  (1 - kappa) * kappa ^ (j - 1)
}



# Marginal Beta -----------------------------------------------------------


log_post_beta <- function(beta, MbarS, n_j, Sbar, ab, bb) {
  Sbar * lgamma(beta) -  sum(lgamma(beta + n_j)) +
    (sum(MbarS) + ab - 1) * log(beta) - bb * beta
}


log_post_LOGbeta <- function(LOGbeta, MbarS, n_j, Sbar, ab, bb) {
  log_post_beta(
    exp(LOGbeta),
    MbarS = MbarS,
    n_j = n_j,
    Sbar = Sbar,
    ab = ab,
    bb = bb
  ) +
    LOGbeta
}

MH.beta <- function(beta, MbarS, n_j, Sbar, ab, bb, sigma.prop) {
  old.log.beta <- log(beta)
  new.log.beta <- rnorm(1, old.log.beta, sigma.prop)
  acc          <- 0
  
  alpha <-
    log_post_LOGbeta(
      LOGbeta = new.log.beta,
      MbarS = MbarS,
      n_j = n_j,
      Sbar = Sbar,
      ab = ab,
      bb = bb
    ) -
    log_post_LOGbeta(
      LOGbeta = old.log.beta,
      MbarS = MbarS,
      n_j = n_j,
      Sbar = Sbar,
      ab = ab,
      bb = bb
    )
  
  if (runif(1) < exp(alpha)) {
    beta <- exp(new.log.beta)
    acc  <- 1
  } else{
    beta <- exp(old.log.beta)
  }
  
  return(c(beta, acc))
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
                verbose=1,                       # Print the iteration progression
                fixedAB=T,                       # Do you want to keep the concentration parameters of the two DPs fixed?
                cheap=F,                         # If T, the output will contain only the membership labels
                kappa=0.5,                       # Slice sampler parameter
                seed=NA,                         # Specify if you want to set the random seed
                sigma.prop.beta=.25,
                batch = 100,
                post.dens = T,
                alpha.fixed = 1,
                beta.fixed = 1,
                conditional.beta = T){                       
  
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
  theta                 <- matrix(NA, L0, 2)
  cij                   <- numeric(N)
  zj                    <- numeric(J)
  BETA_DP               <- numeric(nsim); 
  ALPHA_DP              <- numeric(nsim)
  ACC                   <- numeric(batch)
  Sigma.beta            <- numeric(nsim)
  NN.cz                 <- matrix(NA,nsim,2)
  ind                   <-  0
  PI                    <- vector("list",length = nsim)
  OMEGA                 <- vector("list",length = nsim)
  Csi_ij                <- matrix(NA,nsim,N)
  Z_j                   <- matrix(NA,nsim,J)
  THETA                 <- vector("list",length = nsim)
  if(post.dens){
    xx <- seq(min(y_obser)-5,max(y_obser)+5,by = .05)
    nx <- length(xx)
    DENS <- array(0,c(nx,J,nsim)) 
  }else{
    DENS <- NA
    xx   <- NA
  }
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
  #####################################################################################
  # hyperparameters gamma distribution of alpha and beta of DP
  if(fixedAB){
    alpha   <- alpha.fixed
    beta    <- beta.fixed
  }else{
    alpha   <- rgamma(1, a_alpha, b_alpha)
    beta    <- rgamma(1, a_beta,  b_beta)
        }
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
    Uj   <- runif(J, 0, g_slice((zj),kappa))
    L.z  <- 1 + floor((log(Uj) - log(1 - kappa)) / log(kappa))
    J.z  <- max(zj) #length(unique(zj))#
    NN.z <- max(c(L.z, J.z))
    xi.z <- g_slice(1:NN.z,kappa)
    ################################################################
    # Stochastic Truncation - Observational
    Uij     <- runif(N, 0, g_slice(cij,kappa))
    L.c     <- 1 + floor((log(Uij) - log(1 - kappa)) / log(kappa))
    J.c     <- max(cij) #length(unique(cij))#
    NN.c    <- max(c(L.c, J.c))
    xi.c    <- g_slice(1:NN.c,kappa)
    

    ################################################################
    # Update Distributional Weights
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
    # omega <- Update_omega(
    #   cij = cij,
    #   zj_pg = zj.pg,
    #   NN_c = NN.c,
    #   NN_z = NN.z,
    #   beta = beta
    # )
    v.c <- Update_Observational_Sticks(
        cij = cij,
        zj_pg = zj.pg,
        NN_c = NN.c,
        NN_z = NN.z,
        beta = beta
      )
    if(any(v.c==1)){
    v.c[v.c==1] <- 1-1e-4  # avoid numerical issues
    }
    omega <- apply(v.c,2,SB_given_u2)
    ################################################################
    # Update Atoms
    theta <- Update_theta(y_obser = y_obser,
                          cij = cij,
                          a0 = a0,
                          b0 = b0,
                          k0 = k0,
                          m0 = m0,
                          NN_c = NN.c,
                          J = J
      )
    ################################################################
    # Update Labels
    cij <- Update_Cij(y_obser = y_obser,
                      Uij = Uij,
                      xi_c = xi.c,
                      omega = omega,
                      zj_pg = zj.pg,
                      theta = theta,
                      N = N,
                      NN_c = NN.c)
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
    if(!fixedAB){
    
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
    if(conditional.beta){
      UniqZJ      <- unique(zj)
      sUniqZj     <- sort(UniqZJ,index=T)
      Sbar        <- length(UniqZJ)
      usedV       <- v.c[,sUniqZj$x]
      MbarS       <- tapply(cij, zj.pg, max) # tapply sorts the output according to zj.pg. To not lose the correspondence with usedV, I sort the column of usedV
      one_m_usedV <- log(1-usedV)
      
      astar <- a_beta + sum(MbarS)
      bstar <- b_beta - sum(sapply(1:Sbar, function(x)  sum(one_m_usedV[1:MbarS[x],x])))
    
    # UniqZJ      <- unique(zj)
    # sUniqZj     <- sort(UniqZJ,index=T)
    # Sbar        <- length(UniqZJ)
    # usedV       <- v.c[,UniqZJ]
    # MbarS       <- tapply(cij, zj.pg, max)
    # one_m_usedV <- log(1-usedV)
    # 
    # 
    # astar <- a_beta + sum(MbarS)
    # bstar <- b_beta - sum(sapply(1:Sbar, function(x)  sum(one_m_usedV[1:MbarS[x],sUniqZj$ix[x]])))
    # 
    
      beta <- rgamma(1,astar,bstar) 
    
    }else{
        
        Sbar    <-  length(unique(zj.pg))
        n_j   <-  table(zj.pg)
        MbarS   <-  tapply(cij, zj.pg, function(x) length(unique(x)))
    
        BetaAcc  <- MH.beta(beta = beta, MbarS = MbarS, n_j = n_j, Sbar = Sbar,
                         ab = a_beta, bb=b_beta,sigma.prop = sigma.prop.beta)
    
        beta <- BetaAcc[1]
    
        ACC[sim - ind * batch] = BetaAcc[2]
        if (((sim) %% batch) == 0) {
          sigma.prop.beta       = exp(ifelse(mean(ACC) < .44, #optThresh,
                                  (log(sigma.prop.beta)       - min(
                                    0.01, 1 / sqrt(sim)
                                  )) , (log(sigma.prop.beta)       + min(
                                    0.01, 1 / sqrt(sim)
                                  ))))
          ind <- ind + 1
        }
        }
    }
    ################################################################
    if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
      rr                         <- floor((sim - burn_in)/thinning)
      PI[[rr]]                   <- pi.z
      OMEGA[[rr]]                <- omega
      Csi_ij[rr,]                <- (cij)
      Z_j[rr,]                   <- (zj)
      THETA[[rr]]                <- theta
      ALPHA_DP[rr]               <- alpha
      BETA_DP[rr]                <- beta
      Sigma.beta[rr]             <- sigma.prop.beta
      NN.cz[rr,1]                <- NN.c
      NN.cz[rr,2]                <- NN.z
      if(post.dens){
        allnorm <- apply(theta,1,function(z)dnorm(xx,z[1],sqrt(z[2]))   )
        for(jjj in 1:J){
          DENS[,jjj,rr] <- 
            rowSums(matrix(rep(omega[,zj[jjj]],nx),
                           nrow = nx,byrow = T) * allnorm)
        }

      }
    }
    if (verbose) {
      ipbar <- ipbar + 1
      setTxtProgressBar(pbar, ipbar)
    }
   # cat(length(unique(zj)),"----",length(unique(cij)),"---",NN.z,"+++",NN.c,"\n")
  }
  if(cheap){
    out <- list(Z_j=Z_j, Csi_ij=Csi_ij,
                A_DP=ALPHA_DP, B_DP=BETA_DP,NN.cz=NN.cz,
                y_obser=y_obser, y_group=y_group, nsim=nsim, prior = prior,Sigma.beta=Sigma.beta)
  }else{
    out <- list(Z_j=Z_j, Csi_ij=Csi_ij, PI=PI, THETA=THETA, 
                A_DP=ALPHA_DP, B_DP=BETA_DP,NN.cz=NN.cz,post.dens=DENS,xx=xx,
                OMEGA=OMEGA,  y_obser=y_obser, y_group=y_group, nsim=nsim,Sigma.beta=Sigma.beta, prior = prior)
    
    
  }
  class(out) <- "CAM_SLICE" #Code-name for our model
  return(out)
}
