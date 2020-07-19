############################ Simulation Study 1
# Loading Functions and Libraries -----------------------------------------
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
sourceCpp("CommonAtomModel functions/CAM/CAM.cpp")
source("CommonAtomModel functions/CAM/CAM.R")
##################################################################################
# Data generation ---------------------------------------------------------
YD <- YG <- list()
set.seed(12345689)
med <- c(0,5,10,13,16,20)
CaseA <- function(m,s=sqrt(.6),N){
  p <- c()
  xx <- c()
  for(i in 1:6){
    x <- sample(1:i,N,T)
    xx <- c(xx,x)
    p <- c(p,rnorm(N,m[x],s)) }
  return(cbind(p,xx))
}

NUM_A <- c(25,50,100)

for(i in 1:3){
  y1 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  y2 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  YD[[i]] <- c(y1[,1],y2[,1])
  repet   <- rep(NUM_A[i],12) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG[[i]] <- y_group
  }

plot(YD[[3]],col=YG[[3]])

###############################################################
# Case B
NUM_B <- c(10,20,40)
CaseB0 <- function(m,s=sqrt(.6),N, nsim=100){
  replicate(nsim, {x <- sample(1:N,1)
  cbind(ifelse(x==1,1,0),   
        rnorm(1,m[x],s),x) },simplify = T)
}

CaseB <- function(m,s=sqrt(.6),NUM){
  sa <- c(); CC <- c(); dc <- c()
  for(ll in 1:2){
    for(e in 1:6){
      LALA <- CaseB0(med,N=e,
                    nsim = e*NUM)
      sa <-  c(sa,LALA[2,])
      dc  <- c(dc,LALA[3,])
    }
  }
   return(cbind(sa,dc))
}


for(i in 1:3){
  y1 <- CaseB(m = med,s = sqrt(.6),N = NUM_B[i])
  YD[[3+i]] <- c(y1[,1])
  repet   <- c((1:6)*NUM_B[i],(1:6)*NUM_B[i]) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG[[3+i]] <- y_group
}


##################################################################################
plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])
plot(YD[[4]],col=YG[[4]])
plot(YD[[5]],col=YG[[5]])
plot(YD[[6]],col=YG[[6]])

table(YG[[6]])
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
      verbose.step = 500,fixedAB = F,
      kappa=0.5,cheap=T,seed=123456
  )
}

S1 <- mclapply(1:6,parallel.CAM,mc.cores = 6)
saveRDS(S1,"Scenario_1_CAMoutput.RDS")

# Distributional Clustering -----------------------------------------------

DF_Z <- list()
for(i in 1:6){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S1[[i]]$Z_j,~.x)))
}

PSMs <- map(DF_Z, ~comp.psm(.x))
#D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x)$cl)
D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)
D_CLs
saveRDS(D_CLs,"Scenario_1_DistributionalClustering.RDS")



# Observational Clustering ------------------------------------------------
DF_Csi <- list()
for(i in 1:6){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S1[[i]]$Csi_ij,~.x)))
}

PSMs_2 <- map(DF_Csi, ~comp.psm(.x))
O_CLs  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
saveRDS(O_CLs,"Scenario_1_DistributionalClustering.RDS")

plot(YD[[1]],col=O_CLs[[1]])
plot(YD[[4]],col=O_CLs[[4]])
