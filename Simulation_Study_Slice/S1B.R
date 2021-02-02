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
ALL_s <- readRDS("Simulated_Data/ALL_S1B_100.RDS")
Yall_s   <- ALL_s[[1]]
Gall_s   <- ALL_s[[2]]
Oall_s   <- ALL_s[[3]]

plot(Yall_s[[100]][[1]],col=Gall_s[[1]][[1]])
set.seed(12345)
S1_B_1 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=1)
S1_B_2 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=2)
S1_B_3 <- mclapply(1:30,parallel.CAM,mc.cores = 5,
                   nYD=Yall_s,nYG=Gall_s,K=3)

L <- list(S1_B_1,S1_B_2,S1_B_3)
saveRDS(L,"Simulation_Study_Slice/S1b.RDS")

# Result Extraction -------------------------------------------------------




# Distributional Clusters -------------------------------------------------
gt_distr   <- rep(1:6,2)
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))

time <- list()
attr <- psm <- clu <- list()
perc <- nclu <- aran <- frob <- numeric(30)
#################################################################################
MOD <- S1_B_3
#################################################################################
for(i in 1:30){
  R <-  MOD[[i]]$model$Z_j
  perc[i] <- mean(apply(R,1,function(x) length(unique(x))==6))
  
  psm[[i]]  <- PSM(R)
  clu[[i]]  <- mcclust.ext::minVI(psm[[i]],method = "greedy")$cl
  nclu[i] <- length(unique(clu[[i]]))
  aran[i] <- mcclust::arandi(clu[[i]],gt_distr)
  frob[i] <- StatPerMeCo::Frobenius(psm[[i]],DC_GT_PSM)/maxGTD
}

RES <- cbind(perc,nclu/6,aran,frob)
saveRDS(RES,"Simulation_Study_Slice/results_DC_S1b1.RDS")
plot(ts(RES))
boxplot(RES)
pheatmap::pheatmap(psm[[1]],cluster_rows = F,cluster_cols = F)
plot(frob)
plot(aran)


# Observational Clusters -------------------------------------------------

time <- list()
attr <- psm <- clu <- list()
perc <- nclu <- aran <- frob <- numeric(30)
plot(Yall_s1a[[1]][[1]],col=Oall_s1a[[1]][[1]])
#################################################################################
K   <- 3
MOD <- S1_B_3
MOD <- S1_B_3
MOD <- S1_B_3
#################################################################################
for(i in 1:30){
  gt_distr   <- Oall_s1a[[i]][[K]]
  DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
  maxGTD <- StatPerMeCo::Frobenius(
    matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
    matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))
  
  
  R <-  MOD[[i]]$model$Csi_ij
  perc[i] <- mean(apply(R,1,function(x) length(unique(x))==6))
  
  psm[[i]]  <- PSM(R)
  clu[[i]]  <- mcclust.ext::minVI(psm[[i]])$cl
  nclu[i] <- length(unique(clu[[i]]))
  aran[i] <- mcclust::arandi(clu[[i]],gt_distr)
  frob[i] <- StatPerMeCo::Frobenius(psm[[i]],DC_GT_PSM)/maxGTD
  cat(i)
}


ORES <- cbind(perc,nclu/6,aran,frob)
saveRDS(ORES,"Simulation_Study_Slice/results_OC_S1b1.RDS")
saveRDS(ORES,"Simulation_Study_Slice/results_OC_S1b2.RDS")
saveRDS(ORES,"Simulation_Study_Slice/results_OC_S1b3.RDS")


plot(ts(RES))
boxplot(RES)
pheatmap::pheatmap(psm[[1]])
plot(frob)
plot(aran)




