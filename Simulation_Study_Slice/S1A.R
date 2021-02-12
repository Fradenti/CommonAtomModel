############################ Simulation Study 1
# Loading Functions and Libraries -----------------------------------------
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
sourceCpp("Functions/CAM/newCAM.cpp")
source("Functions/CAM/newCAM.R")
# Creation of parallelizable function -------------------------------------
parallel.CAM <- function(i,K,nYD,nYG){
  
  Y <- nYD[[i]][[K]]
  G <- nYG[[i]][[K]]
  t1 <- Sys.time()
  OBJ <- CAM(y_obser = Y,
      y_group = G,
      K0 = 10,
      L0 = 10,
      prior = list(
        # hyperparameters NIG
        m0=mean(Y), k0=1/var(Y), a0=3, b0=1,
        # hyperparameters alpha and beta
        a_alpha=3, b_alpha = 3,
        a_beta =3, b_beta = 3 ),
      nsim = 10000,
      burn_in = 10000,
      thinning = 1,verbose = 1,
      fixedAB = T,
      kappa=0.5,cheap=T,seed=1234*i
  )
  t2 <- Sys.time()
  return(list(model=OBJ, time=t2-t1))
}
##################################################################################
# Data generation ---------------------------------------------------------
ALL_s <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Yall_s   <- ALL_s[[1]]
Gall_s   <- ALL_s[[2]]
Oall_s   <- ALL_s[[3]]

plot(Yall_s[[100]][[1]],col=Gall_s[[1]][[1]])
set.seed(12345)
S1_A_1 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=1)
S1_A_2 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=2)
S1_A_3 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=3)

L <- list(S1_A_1,S1_A_2,S1_A_3)
saveRDS(L,"Simulation_Study_Slice/S1a.RDS")















# Result Extraction -------------------------------------------------------

L <- readRDS("Simulation_Study_Slice/S1a.RDS")



# Distributional Clusters -------------------------------------------------
gt_distr   <- rep(1:6,2)
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))

time      <- timescale <- matrix(NA,30,3)
#################################################################################
for(K in 1:3){
  mod <- L[[K]]
  time[,K] <- as.numeric(unlist(map(mod,~.x$time)))
  timescale[,K] <- unlist(map(mod,~attr(.x$time,which = "units")))
}

Time <- ifelse(timescale=="secs", time/60, time)
boxplot(Time)
saveRDS(Time,"Simulation_Study_Slice/Time1A.RDS")


nclu <- aran <- frob <- matrix(NA,30,3)
PSMs <- psm <- list()
CLUs <- clu <- list()

for(K in 1:3){
  mod <- L[[K]]
  

  for(i in 1:30){
    
    R <-  mod[[i]]$model$Z_j
    
    psm[[i]]  <- PSM(R)
    clu[[i]]  <- mcclust.ext::minVI(psm[[i]],method = "greedy")$cl
    
    nclu[i,K] <- length(unique(clu[[i]]))
    aran[i,K] <- mcclust::arandi(clu[[i]],gt_distr)
    frob[i,K] <- StatPerMeCo::Frobenius(psm[[i]],DC_GT_PSM)/maxGTD
    cat(i)
  }
  rm(mod)
  
  

}
RES <- list(nclu,aran,frob)
saveRDS(RES,paste0("Simulation_Study_Slice/Results_30_1A_Distributional.RDS"))


# Observational Clusters -------------------------------------------------
ALL_S1a  <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Yall_s1a <- ALL_S1a[[1]]
Gall_s1a <- ALL_S1a[[2]]
Oall_s1a <- ALL_S1a[[3]]

nclu <- aran <- frob <- matrix(NA,30,3)
PSMs <- psm <- list()
CLUs <- clu <- list()

#################################################################################


for(K in 1:3){
  mod <- L[[K]]
  

for(i in 1:30){
  gt_distr   <- Oall_s1a[[i]][[K]]
  DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
  maxGTD <- StatPerMeCo::Frobenius(
    matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
    matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))
  
  
  R <-  mod[[i]]$model$Csi_ij
  
  psm[[i]]  <- PSM(R)
  clu[[i]]  <- mcclust.ext::minVI(psm[[i]])$cl
  
  nclu[i,K] <- length(unique(clu[[i]]))
  aran[i,K] <- mcclust::arandi(clu[[i]],gt_distr)
  frob[i,K] <- StatPerMeCo::Frobenius(psm[[i]],DC_GT_PSM)/maxGTD
  cat(i)
  }

}


RES <- list(nclu,aran,frob)
saveRDS(RES,paste0("Simulation_Study_Slice/Results_30_1A_Observational.RDS"))



