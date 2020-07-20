############################ Simulation Study 1
# Loading Functions and Libraries -----------------------------------------
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
sourceCpp("CommonAtomModel functions/DCAM/DCAM.cpp")
source("CommonAtomModel functions/DCAM/DCAM.R")
# Creation of parallelizable function -------------------------------------
parallel.DCAM <- function(i){
  DCAM(y_obser = YD[[i]],y_group = YG[[i]],
      K0 = 10,L0 = 10,
      prior = list(
        # hyperparameters NIG
        m0=0, k0=1/var(YD[[i]]), a0=3, b0=1,
        # hyperparameters alpha and beta
        a_alpha=3, b_alpha = 3,
        a_beta =3, b_beta = 3 ),
      NSIM = 50000,burn_in = 50000,thinning = 1,
      verbose.step = 50,fixedAB = F,
      kappa=0.5,cheap=T,seed=123456
  )
}
##################################################################################
# Data generation ---------------------------------------------------------
YD <- YG <- list()
set.seed(123456)
num <- c(5,10,15,25,50,75,100)
YG <- YD <- list()
for(i in 1:6){
  YD[[i]] <- c( rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(10, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(10, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T))
  repet   <- rep(100+num[i],10)
  YG[[i]] <- c(rep(1:length(repet),repet))
}

plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])
plot(YD[[4]],col=YG[[4]])
plot(YD[[5]],col=YG[[5]])
plot(YD[[6]],col=YG[[6]])

##################################################################################
# Setting up the model -----------------------------------------------------

DCAM_Res <- mclapply(1:6,parallel.DCAM,mc.cores = 6)
saveRDS(DCAM_Res,"Scenario3_DCAMoutput.RDS")

# Distributional Clustering -----------------------------------------------

DF_Z <- list()
for(i in 1:6){
  DF_Z[[i]] <- t(as.matrix(map_dfc(DCAM_Res[[i]]$Z_j,~.x)))
}

PSMs <- map(DF_Z, ~comp.psm(.x))
D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x)$cl)
D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)
D_CLs
saveRDS(D_CLs,"Scenario_3_DistributionalClustering.RDS")



# Observational Clustering ------------------------------------------------
DF_Csi <- list()
for(i in 1:6){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(DCAM_Res[[i]]$Csi_ij,~.x)))
}

PSMs_2 <- map(DF_Csi, ~comp.psm(.x))
O_CLs  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
saveRDS(O_CLs,"Scenario_3_ObservationalClustering.RDS")

plot(YD[[1]],col=O_CLs[[1]])
plot(YD[[6]],col=O_CLs[[6]])




# Nested_DP ---------------------------------------------------------------
YD <- YG <- list()
set.seed(123456)
num <- c(5,10,15,25,50,75,100)
YG <- YD <- list()
for(i in 1:6){
  YD[[i]] <- c( rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(10, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T),
                rep(0,50),rep(1,50),sample(10, num[i],T),
                rep(0,50),rep(1,50),sample(100,num[i],T),
                rep(0,50),rep(1,50),sample(50, num[i],T))
  repet   <- rep(100+num[i],10)
  YG[[i]] <- c(rep(1:length(repet),repet))
}

plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])
plot(YD[[4]],col=YG[[4]])
plot(YD[[5]],col=YG[[5]])
plot(YD[[6]],col=YG[[6]])

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


NDP_S3 <- newEnDP_Gibbs(y_obser = YD[[6]], y_group = YG[[6]],
                      K = K,L = L,
                      prior = prior,
                      NSIM = NSIM,
                      burn_in = burn_in,thinning = thinning,
                      verbose = T,verbose.step = 50,monitor = T,
                      m.step = 500,fixedAB = F,restart = F,
                      restartList = NULL,
                      RDG=T,cheap = T,seed = 12345)

