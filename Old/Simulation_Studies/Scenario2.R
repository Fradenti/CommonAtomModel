############################ Simulation Study 2
  # Loading Functions and Libraries -----------------------------------------
  library(Rcpp)
  library(tidyverse)
  library(scales)
  library(mcclust)
  library(parallel)
  sourceCpp("CommonAtomModel/Functions/CAM/CAM.cpp")
  source("CommonAtomModel/Functions/CAM/CAM.R")
  ##################################################################################
  # Data generation ---------------------------------------------------------
  set.seed(1234568)
  NNN     <- 40
  YD      <- YG      <- list()
  for(i in 1:6){
    YD[[i]] <- c(replicate(i,{c( rnorm(NNN*.75,0,sqrt(.6)),rnorm(NNN*.25, 3,sqrt(.6)),
                               rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.75, 3,sqrt(.6)),
                               rnorm(NNN*.33,0,sqrt(.6)),rnorm(NNN*.34, 2,sqrt(.6)),
                               rnorm(NNN*.33,-2,sqrt(.6)),
                               rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.25, -2,sqrt(.6)),
                               rnorm(NNN*.25+1,2,sqrt(.6)),  rnorm(NNN*.25, 10, sqrt(.6)))}))
  rep     <- rep(NNN,4*i)
  YG[[i]] <- rep(1:length(rep),rep)
  }
  #Seed <- list(1234,3131,2121,4141,5151,5647)
  plot(YD[[1]],col=YG[[1]])
  plot(YD[[2]],col=YG[[2]])
  plot(YD[[4]],col=YG[[4]])
  plot(YD[[6]],col=YG[[6]])
  
  
  ##################################################################################
  # Settin up the model -----------------------------------------------------
  
  
  # Creation of parallelizable function -------------------------------------
  parallel.CAM <- function(i){
    CAM(y_obser = YD[[i]],y_group = YG[[i]],
        K0 = 10,L0 = 10,
        prior = list(
          # hyperparameters NIG
          m0=0, k0=1/var(YD[[i]]), a0=3, b0=1,
          # hyperparameters alpha and beta
          a_alpha=3, b_alpha = 3,
          a_beta =3, b_beta = 3 ),
        NSIM = 50000,burn_in = 50000,thinning = 1,
        verbose.step = 50000,fixedAB = F,
        kappa=0.5,cheap=T,seed=1324354
        )
  }
  
  R <- mclapply(1:6,parallel.CAM,mc.cores = 6)
  saveRDS(R,"Scenario_2_CAMoutput_betaok.RDS")
  
  # Distributional Clustering -----------------------------------------------
  
  DF_Z <- list()
  for(i in 1:6){
  DF_Z[[i]] <- t(as.matrix(map_dfc(R[[i]]$Z_j,~.x)))
  }
  
  PSMs <- map(DF_Z, ~comp.psm(.x))
  D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)
  D_CLs
  saveRDS(D_CLs,"Scenario_2_DistributionalClustering_betaok.RDS")
  
  plot(YD[[1]],col=YG[[1]])
  
  # Observational Clustering ------------------------------------------------
  DF_Csi <- list()
  for(i in 1:6){
    DF_Csi[[i]] <- t(as.matrix(map_dfc(R[[i]]$Csi_ij,~.x)))
  }
  
  PSMs_2 <- map(DF_Csi, ~comp.psm(.x))
  O_CLs  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
  saveRDS(O_CLs,"Scenario_2_DistributionalClustering_betaok.RDS")
  
  plot(YD[[1]],col=O_CLs[[1]])
  plot(YD[[5]],col=O_CLs[[5]])
    
  
  
  
  

# NestedDP ----------------------------------------------------------------

  set.seed(1234568)
  NNN     <- 40
  YD      <- YG      <- list()
  for(i in 1:6){
    YD[[i]] <- c(replicate(i,{c( rnorm(NNN*.75,0,sqrt(.6)),rnorm(NNN*.25, 3,sqrt(.6)),
                                 rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.75, 3,sqrt(.6)),
                                 rnorm(NNN*.33,0,sqrt(.6)),rnorm(NNN*.34, 2,sqrt(.6)),
                                 rnorm(NNN*.33,-2,sqrt(.6)),
                                 rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.25, -2,sqrt(.6)),
                                 rnorm(NNN*.25+1,2,sqrt(.6)),  rnorm(NNN*.25, 10, sqrt(.6)))}))
    rep     <- rep(NNN,4*i)
    YG[[i]] <- rep(1:length(rep),rep)
  }
  plot(YD[[6]],col=YG[[6]])
  
  ##################################################################################
  
  Rcpp::sourceCpp("CommonAtomModel functions/NDP_via_EnhNDP/newEnDP_2018.cpp")
  source("CommonAtomModel functions/NDP_via_EnhNDP/newEnDP_2018.R")
  
  NSIM=50000;thinning=1; burn_in=50000; verbose.step=10; 
  L=30; K=25; verbose=TRUE; monitor=T; 
  m.step=100; fixedAB=F; 
  restart=F; restartList=NA
  
  prior <- list(
    # hyperparameters NIG
    m0=0, k0=1/var(YD[[6]]), a0=3, b0=1,
    # hyperparameters alpha and beta
    a_alpha=3, b_alpha = 3,
    a_beta =3, b_beta = 3, A=1, B=1 )
  
  
  NDP_S2 <- newEnDP_Gibbs(y_obser = YD[[6]], y_group = YG[[6]],
                          K = K,L = L,
                          prior = prior,
                          NSIM = NSIM,
                          burn_in = burn_in,thinning = thinning,
                          verbose = T,verbose.step = 50,monitor = T,
                          m.step = 500,fixedAB = F,restart = F,
                          restartList = NULL,
                          RDG=T,cheap = T,seed = 12345)
  
  
  
  