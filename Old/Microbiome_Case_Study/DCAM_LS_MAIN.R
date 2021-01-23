library(plyr); library(Rcpp); library(RcppArmadillo)
library(reshape); library(pryr); library(knitr); library(gridExtra)
library(tidyverse); library(compiler);  library(MCMCpack); library(scales)
#######################################################################################################

# Auxiliary functions -----------------------------------------------------
reset <- function(x){
  z <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}


# Main Function -----------------------------------------------------------

DCAM_LibSize <- function(y_obser, y_group, 
                         K0=10, L0=20, 
                         prior, 
                         NSIM=1000, burn_in=1000, thinning=1, 
                         verbose=TRUE, verbose.step=15, 
                         fixedAB=T, Fa=1, Fb=1,restart=F, 
                         kappa=0.5,warm.start=T,
                         Wgamma=1, User.defined.gammas=NULL) {
  
  data    <- as.tibble(cbind(y_obser,y_group)) %>% arrange(y_group)
  J       <- length(unique(data$y_group)) # number of groups
  nj      <- table(data$y_group)          # number of observations inside each group
  N       <- length(y_obser)
  
  if(Wgamma){
    if(is.null(User.defined.gammas)){  
      gamma <- rep(tapply(y_obser,y_group,sum), as.numeric(nj) )
      gamma <- gamma/min(gamma)
    }else{
      gamma <- User.defined.gammas
    }}
  
  Cnj     <- cumsum(nj)
  y_lat   <- y_obser+runif(N)
  ############################### hyperparameters extraction
  # Prior for mu and sigma
  m0 <- prior$m0; k0 <- prior$k0; 
  a0 <- prior$a0; b0 <- prior$b0;
  # hyperparameters gamma distribution of alpha and beta of DP
  BETA_DP <- numeric(NSIM); ALPHA_DP <- numeric(NSIM)
  a_alpha <- prior$a_alpha;  b_alpha <- prior$b_alpha
  a_beta  <- prior$a_beta;   b_beta  <- prior$b_beta
  alpha   <- rgamma(1, a_alpha, b_alpha)
  beta    <- rgamma(1, a_beta,  b_beta )
  #####################################################################################
  # Initialization of Output:
  # THETA ##########################################################################################################################
  tsl0 <- matrix(NA,L0,2)
  cij  <- numeric(N)
  zj   <- numeric(J)
  #####################################################################################
  if(warm.start){
    km.out    <- kmeans(y_obser,L0)
    cij       <- km.out$cluster
    tsl0[,1]  <- tapply(y_obser, cij, mean)
    tsl0[,2]  <- tapply(y_obser, cij, sd)
    if(K0>=J){
      zj        <- 1:J
    }else{
      km.outZ   <- kmeans(tapply(y_obser, y_group, mean),K0)
      zj        <- km.outZ$cluster
    }
  }else{
    zj        <- sample(K0,J,T)
    cij       <- sample(L0,N,T)
    tsl0[,1]  <- tapply(y_obser, cij, mean)
    tsl0[,2]  <- tapply(y_obser, cij, sd)
  }
  zj.pg <- rep(zj,nj)
  tsl0.tmp <- matrix(0,N,2)
  #####################################################################################
  g <- function(j) (1-kappa)*kappa^(j-1)
  #####################################################################################
  acc1 <- acc2 <- numeric((NSIM*thinning + burn_in))
  #####################################################################################
  # Containers
  pi_star_k           <- vector("list",length = NSIM)
  omega_star_lk       <- vector("list",length = NSIM)
  Csi_ij              <- vector("list",length = NSIM)
  Z_j                 <- vector("list",length = NSIM)
  theta_star_zero     <- vector("list",length = NSIM)
  #YLAT               <- matrix(NA,NSIM,N)
  
  for(sim in 2:(NSIM*thinning + burn_in)){
    ################################################################   
    Uj  <- runif(J, 0, g(zj))
    L.z  <- 1 + floor((log(Uj) - log(1 - kappa)) / log(kappa))
    J.z  <- max(zj)
    NN.z <- max(c(L.z, J.z))
    xi.z <- g(1:NN.z)
    ################################################################  
    Uij     <- runif(N, 0, g(cij))
    L.c     <- 1 + floor((log(Uij) - log(1 - kappa)) / log(kappa))
    J.c     <- max(cij)
    NN.c    <- max(c(L.c, J.c))
    xi.c    <- g(1:NN.c)
    ################################################################
    if(NN.c > J.c) tsl0 <- rbind(tsl0, matrix( 0,(NN.c - J.c),2))
    ################################################################
    v.z <- pi.z  <- numeric(NN.z)
    v.omega <- omega <- matrix(0, NN.c, NN.z)
    v.z  <- c(UPD_Pi_v_z(zj = zj, NN_z = NN.z, alpha = alpha))
    ################################################################
    if(length(unique(zj))>1) {pi.z <- SB_given_u2(v.z)}else{pi.z <- v.z}
    ################################################################
    tsl1 <- tsl0 <- UPD_tsl0(y_LAT = y_lat/gamma, 
                             cij = cij,a0 = a0,b0 = b0,k0 = k0,m0 = m0,NN_c = NN.c,J = J) 
    ################################################################
    omega <- UPD_omega(cij = cij,zj_pg = zj.pg,NN_c = NN.c,NN_z = NN.z,beta = beta)
    oldzj <- c(UPD_Zj_collapsed_uij(Uj = Uj,xi_z = xi.z,xi_c = xi.c,
                                    pi_z = pi.z,cij = cij,omega = omega,
                                    y_group = y_group,NN_z = NN.z,J = J))
    ############################### aggiunta solo sta linea
    zj <- reset(oldzj)
    omega[,unique(zj)] = omega[,unique(oldzj)] 
    zj.pg <- rep(zj,nj)
    ######################################################
    tsl0.tmp <- UPD_tsl0_for_cij(y_obser = y_obser, 
                                 Uij = Uij, xi_c = xi.c, omega = omega,
                                 zj_pg = zj.pg, tsl0 = tsl0, N = N, NN_c = NN.c,gamma = gamma)
    tsl0   <- subset(tsl0.tmp,!duplicated(tsl0.tmp[,1]))
    cij    <- apply(outer(tsl0.tmp[,1], tsl0[,1], "-"), 1, function(x) which(x == 0));
    ################################################################
    k      <- length(unique(zj))
    eta    <- rbeta(1,alpha+1,J)  
    Q      <- (a_alpha+k-1)/(J*(b_alpha-log(eta)))
    pi_eta <- Q/(1+Q)  
    alpha  <- ifelse(runif(1)<pi_eta,  rgamma(1,a_alpha+k,   b_alpha-log(eta)), 
                     rgamma(1, a_alpha+k-1,b_alpha-log(eta))  )
    ################################################################
    ################################################################
    k2      <- length(unique(cij)) # double check this sum(tapply(cij, zj.pg, function(x) length(unique(x)))) # 
    eta2    <- rbeta(1,beta+1,N)
    Q2      <- (a_beta+k2-1)/(N*(b_beta-log(eta2)))
    pi_eta2 <- Q2/(1+Q2)
    beta    <- ifelse(runif(1)<pi_eta2,  
                      rgamma(1,a_beta+k2,  b_beta-log(eta2)),
                      rgamma(1,a_beta+k2-1,b_beta-log(eta2))  )
    ################################################################
    y_lat   <- Upd_LAT(y_obser = y_obser,cij = cij,tsl0 = tsl0,N = N,gamma = gamma)
    
    
    if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
      rr                      <- floor((sim - burn_in)/thinning)
      pi_star_k[[rr]]         <- pi.z
      omega_star_lk[[rr]]     <- omega
      Csi_ij[[rr]]            <- cij
      Z_j[[rr]]               <- zj
      theta_star_zero[[rr]]   <- tsl1############## tsl0 non va bene, non comparabile con i pi grechini
      ALPHA_DP[rr]            <- alpha
      BETA_DP[rr]             <- beta
      # YLAT[rr,]               <- y_lat
    }
    
    if (sim%%(verbose.step*thinning) == 0) {
      cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in, "---", length(unique(zj)), "\n",
                sep = "" ))}
    
  }
  
  
  
  out <- list(Z_j=Z_j, Csi_ij=Csi_ij, pi_star_k=pi_star_k, theta_star_zero=theta_star_zero, 
              A_DP=ALPHA_DP, B_DP=BETA_DP,acc1=acc1, acc2=acc2,gamma=gamma,
              #ylat=YLAT,
              omega_star_lk=omega_star_lk,  y_obser=y_obser, y_group=y_group, NSIM=NSIM)
  class(out) <- "DCAM_LS" #Code-name for our model
  return(out)
}
