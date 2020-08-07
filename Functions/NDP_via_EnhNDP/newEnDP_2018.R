
# NSIM=1000;thinning=1; burn_in=100; verbose.step=10; L=20; K=10; verbose=TRUE; monitor=T; m.step=100; fixedAB=F; EPS=.0; restart=F; restartList=NA



library(plyr); library(Rcpp); library(RcppArmadillo)
library(reshape); library(pryr); library(knitr)
library(tidyverse); library(compiler);  library(MCMCpack)
#sourceCpp("C:\\Users/f.denti2/Dropbox/IRVINE_RESEARCH/CAMERLENGHI_EnhancedNDP/eNdp_september/newEnDP.cpp")


# Superposition of random measures (in Intro Thesis Francesco Denti) --------

newEnDP_Gibbs <- function(y_obser, y_group, K=20, L=50, prior, NSIM=100,
                       burn_in, thinning, verbose=TRUE, verbose.step=15, 
                       how_many_cluster_to_start=3,
                       monitor=T, m.step=100, fixedAB=T, restart=F, 
                       restartList, RDG=F, cheap=F, seed){
  
  if(!is.na(seed)){
    set.seed(seed)
  }
  
  par(mfrow=c(1,1))
  # Preparing containers + dimensions checkng
  data    <- as_tibble(cbind(y_obser,y_group)) %>% arrange(y_group)
  nij     <- nrow(data)                   # number of observations
  J       <- length(unique(data$y_group)) # number of groups
  nj      <- table(data$y_group)          # number of observations inside each group
  max_nj  <- max(nj)
  
  ############################### hyperparameters extraction
  # Prior for mu and sigma
  m0 <- prior$m0; k0 <- prior$k0; 
  a0 <- prior$a0; b0 <- prior$b0;
  # hyperparameters gamma distribution of alpha and beta of DP
  a_alpha <- prior$a_alpha;  b_alpha <- prior$b_alpha
  a_beta  <- prior$a_beta;   b_beta  <- prior$b_beta
  # alpha   <- prior$alpha;    beta    <- prior$beta
  rho_A <- prior$A
  rho_B <- prior$B
  alpha <- rgamma(1, a_alpha, b_alpha)
  beta  <- rgamma(1, a_beta,  b_beta )
  
  # Initialization of Output:
  # THETA ##########################################################################################################################
  theta_star_lk       <- array(NA, dim = c(L,K,2,NSIM))
  tslk                <- array(NA, dim = c(L,K,2))
  theta_star_zero     <- array(NA,dim=c(L,2,NSIM))
  tsl0                <- matrix(NA,L,2)
  
  if(restart){
    tslk <- restartList$tslk
    tsl0 <- restartList$tsl0
    }else{
  tslk[,,2]   <-  rinvgamma(L*K, prior$a0, prior$b0)
  tslk[,,1]   <-  rt(L*K, 2*prior$a0)*sqrt(prior$b0/(prior$k0*prior$a0))+prior$m0
  tsl0[,2]    <-  rinvgamma(L, prior$a0, prior$b0)
  tsl0[,1]    <-  rt(L, 2*prior$a0)*sqrt(prior$b0/(prior$k0*prior$a0))+prior$m0
  }
  
  theta_star_lk[,,1,1] <- tslk[,,1]
  theta_star_lk[,,2,1] <- tslk[,,2]
  theta_star_zero[,,1] <- tsl0
  
  
  # PI_STAR ##########################################################################################################################
  # SAMPLE DISTRIBUTIONAL WEIGHTS A PRIORI
  pi_star_k            <- matrix(NA,K,NSIM)
  if(restart){
    psk                <- c(restartList$psk)
  }else{
    psk                <- SB(N = K, 1,1)$pi
  }
  pi_star_k[,1]        <- psk
  
  # OMEGA_STAR ##########################################################################################################################
  # SAMPLE COMMON OBSEVATIONAL WEIGHTS A PRIORI
  
  Omega_l0             <- matrix(NA,L,NSIM)
  if(restart){
    ols0               <- c(restartList$osl0)
  }else{
    osl0               <- SB(N = L, 1,1)$pi
  }
  Omega_l0[,1]         <- osl0
  
  ####
  # SAMPLE OBSERVATIONAL WEIGHTS A PRIORI
  
  omega_star_lk        <- array(NA, dim=c(L,K,NSIM))
  if(restart){
   oslk    <- restartList$oslk
  }else{ 
   oslk    <- replicate(K,SB(N = L, 1, 1)$pi)
  }
   omega_star_lk[,,1]            <- oslk
  
  ####
  # SAMPLE DISTRIBUTIONAL Labels A PRIORI
  Z_j          <- matrix(NA, J, NSIM)      # distributional label
  if(restart){
  zj      <- restartList$zj
  }else{
    zj      <- sample(1:K, J, replace = T, prob = psk)
  }

  ####
  # SAMPLE OBSERVATIONAL Labels A PRIORI
  # Dovrebbe essere C_ij, lo lascio lk per quieto vivere al momento
  Csi_ij       <- array(0, dim=c(nij,NSIM)) # observational label
  if(restart){
   cij <- restartList$cij 
  }else{
  cij <- c()
  for(j in 1:J){
    cij <- c(cij,sample((1:(2*(L))), nj[j], replace = T, prob = c(oslk[,1],osl0)))
  }
  }
  
  Csi_ij[,1] <- cij
  
  Rho <- matrix(NA,NSIM,K)
  rho <- rbeta(K,rho_A,rho_B)
  Rho[1,] <- rho
  
  Xi <- matrix(NA,NSIM,nij)
  xi <- ifelse(cij>L,1,0)
  Xi[1,] <- xi
  
############################## GIBBS SAMPLER LOOP ################################################
  for(sim in 2:(NSIM*thinning + burn_in)){
    
    ####################################################################################################################################################################
    # STEP 1 - sample CENTER indicators
    ####################################################################################################################################################################for(j in 1:J){
    
    zj <- c(UPD_DISTR_LAB(y_obser = y_obser,
                          y_group = y_group, 
                          tslk = tslk, 
                          tsl0 = tsl0, 
                          psk = psk, 
                          rho = rho, 
                          oslk = oslk, 
                          osl0 = osl0, 
                          J =  J, K = K, L = L))
    
    ####################################################################################################################################################################
    # STEP 2 - sample INDIVIDIUAL indicators
    ####################################################################################################################################################################
    
    MAMA <- UPD_OBSER_LAB_XI(y_obser = y_obser,
                             y_group = y_group, 
                             tslk = tslk, 
                             tsl0 = tsl0, 
                             oslk = oslk, 
                             osl0 = osl0, 
                             rho = rho, 
                             J =J, L=L, K=K, NN = nij,
                             zj=zj)
    xi   <- MAMA[,2]
    cij  <- MAMA[,1]
    
    ####################################################################################################################################################################
    # STEP 3, uptading Pi^*
    ####################################################################################################################################################################
    
    res_pistar <- updt_pistar_rcpp(zj=zj,alpha = alpha,K = K)
    psk        <- res_pistar$pistar
    uk         <- res_pistar$uk
    
    ####################################################################################################################################################################
    # INTERVAL: updating useful quantities
    ####################################################################################################################################################################
     
    ListaQuantities <- Int_cpp(y_obser = y_obser, y_group = y_group, cij = cij, zj=zj, nj=nj, K=K, L=L)
    
    ####################################################################################################################################################################
    # STEP 4, uptading w^*_lk
    ####################################################################################################################################################################
    
    res_ostar <- updt_omegastar_lk(nlk=ListaQuantities$nlk,beta = beta)
    oslk      <- res_ostar$os
    v         <- res_ostar$v
    
    res_ostar0 <- updt_omegastar_0(nl0=ListaQuantities$nl0,L = L, alpha = beta)
    osl0       <- res_ostar0$pistar
    v0         <- res_ostar0$uk
    
    ####################################################################################################################################################################
    # Step 5  --  Update Theta
    ####################################################################################################################################################################
    
    tslk <- UPD_TH_CD_cpp_T_lk(AIFD = ListaQuantities,k0 = k0, m0 = m0, b0 = b0, a0 = a0, L = L,K = K)
    tsl0 <- UPD_TH_CD_cpp_T_0( AIFD = ListaQuantities,k0 = k0, m0 = m0, b0 = b0, a0 = a0, L = L,K = K)
    
    ####################################################################################################################################################################
    # Step 6 - Simply updating alpha and beta
    ####################################################################################################################################################################
    
    if(fixedAB & sim == 2){
      alpha <- prior$alpha_DP ; 
      beta  <- prior$beta_DP
    }else{
      alpha <- 1#rgamma(1, a_alpha+(K-1),   b_alpha -  sum(log(1-uk[-K])))
      beta  <- 1#rgamma(1, a_beta +K*(L-1), b_beta  -  sum(log(1-v[-L,])))
    }
    
    ####################################################################################################################################################################
    # Step 7 - Update Rhos
    ####################################################################################################################################################################
    
    if(RDG){
      rho <- rep(0,K)
    }else{
      rho <- upd_rho_j(xi = xi, lab_g = ListaQuantities$A[, 4], A = rho_A, B = rho_B, K = K)
    }
    
    ################################################################################################################################
    # Storing the results
    ################################################################################################################################
    
    if (sim > burn_in & ((sim - burn_in) %% thinning == 0)) {
      rr                      <- floor((sim - burn_in)/thinning)
      pi_star_k[,rr]          <- psk
      omega_star_lk[,,rr]     <- oslk
      theta_star_lk[,,,rr]    <- tslk
      Csi_ij[,rr]             <- cij
      Z_j[,rr]                <- zj
      Xi[rr,]                 <- xi
      Omega_l0[,rr]           <- osl0
      theta_star_zero[,,rr]   <- tsl0
      Rho[rr,]                <- rho
    #  if(sim %% m.step==0) plot(data$y_obser, col=rep(zj,nj)) #plot(BETAREG,type="l")
      }
    # Printing the results
    if(sim==NSIM){
      RL <- list(psk  =  psk,
                          oslk =  oslk,
                          tslk =  tslk,
                          cij  =  cij,
                          zj   =  zj,
                          xi=xi, osl0=osl0, rho=rho )
    }
    
    if (verbose) {
      if (sim%%(verbose.step*thinning) == 0) {
        cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in, "---", length(unique(zj)), "\n",
                  sep = "" ))}}
  }
  if(cheap){
  out <- list(Z_j=Z_j, Csi_ij=Csi_ij, RL=RL, y_obser=y_obser, y_group=y_group, NSIM=NSIM)
  }else{
  out <- list(Z_j=Z_j, Csi_ij=Csi_ij, pi_star_k=pi_star_k,rho=Rho, theta_star_zero=theta_star_zero, Omega_l0=Omega_l0,Xi=Xi,
              omega_star_lk=omega_star_lk, theta_star_lk=theta_star_lk, RL=RL, y_obser=y_obser, y_group=y_group, NSIM=NSIM)
  }
  class(out) <- "Enhanced_ndp" #Code-name for our model
  return(out)
}

#######################################################################################################
reset <- function(x){
  z <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}

#######################################################################################################################

assessing_distrib_cluster <- function(out){
  
  cc  <- apply(out$Z_j,2,function(x) length(unique(x)))
  mat <- as.matrix(table(cc)/length(cc))
  colnames(mat) <- "Probability";
  print(kable(mat,row.names = TRUE, digits = 3))
  raster_plot_given_group_labels(out$Z_j)

  }
