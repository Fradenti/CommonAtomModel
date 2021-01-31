############################ Simulation Study 1
# Loading Functions and Libraries -----------------------------------------
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
sourceCpp("Functions/CAM/newCAM.cpp")
source("Functions/CAM/newCAM.R")
# Creation of parallelizable function -------------------------------------
parallel.CAM <- function(i,nYD,nYG,K){
  
  Y <- Yall_s1a[[i]][[K]]
  G <- Gall_s1a[[i]][[K]]
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
      nsim = 500,
      burn_in = 500,
      thinning = 1,verbose = 1,
      fixedAB = T,
      kappa=0.5,cheap=T,seed=1234*i
  )
  t2 <- Sys.time()
  return(list(model=OBJ, time=t2-t1))
}
##################################################################################
# Data generation ---------------------------------------------------------
ALL_s1a <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Yall_s1a   <- ALL_s1a[[1]]
Gall_s1a   <- ALL_s1a[[2]]
Oall_s1a   <- ALL_s1a[[3]]

plot(Yall_s1a[[100]][[1]],col=Gall_s1a[[1]][[1]])

S1_A_1 <- mclapply(1:30,parallel.CAM,mc.cores = 3,
                   nYD=Yall_s1a,nYG=Gall_s1a,K=1)
S1_A_2 <- mclapply(1:30,parallel.CAM,mc.cores = 3,
                   nYD=Yall_s1a,nYG=Gall_s1a,K=2)
S1_A_3 <- mclapply(1:30,parallel.CAM,mc.cores = 3,
                   nYD=Yall_s1a,nYG=Gall_s1a,K=3)




# Result Extraction -------------------------------------------------------


gt_distr   <- rep(1:6,2)
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))



time <- list()
attr <- psm <- clu <- list()
perc <- nclu <- aran <- frob <- numeric(30)


for(i in 1:30){
  R <-  S1_A_3[[i]]$model$Z_j
  perc[i] <- mean(apply(R,1,function(x) length(unique(x))==6))
  
  psm[[i]]  <- PSM(R)
  clu[[i]]  <- mcclust.ext::minVI(psm[[i]],method = "greedy")$cl
  nclu[i] <- length(unique(clu[[i]]))
  aran[i] <- mcclust::arandi(clu[[i]],gt_distr)
  frob[i] <- StatPerMeCo::Frobenius(psm[[i]],DC_GT_PSM)/maxGTD
  }

RES <- cbind(perc,nclu/6,aran,frob)
plot(ts(RES))
boxplot(RES)
pheatmap::pheatmap(psm[[1]])
plot(frob)
plot(aran)

