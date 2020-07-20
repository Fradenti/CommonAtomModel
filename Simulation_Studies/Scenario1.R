############################ Simulation Study 1
# Loading Functions and Libraries -----------------------------------------
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
sourceCpp("CommonAtomModel functions/CAM/CAM.cpp")
source("CommonAtomModel functions/CAM/CAM.R")
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
      verbose.step = 50,fixedAB = F,
      kappa=0.5,cheap=T,seed=123456
  )
}
##################################################################################
# Data generation ---------------------------------------------------------
YD <- YG <- list()
set.seed(123456)
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

NUM_A <- c(25,50,75)

for(i in 1:3){
  y1 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  y2 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  YD[[i]] <- c(y1[,1],y2[,1])
  repet   <- rep(NUM_A[i],12) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG[[i]] <- y_group
  }


plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])

S1_A <- mclapply(1:3,parallel.CAM,mc.cores = 3)
saveRDS(S1_A,"Scenario_1A_CAMoutput.RDS")

# Distributional Clustering -----------------------------------------------

DF_Z <- list()
for(i in 1:3){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S1_A[[i]]$Z_j,~.x)))
}

PSMs <- map(DF_Z, ~comp.psm(.x))
#D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x)$cl)
D_CLsA  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)
D_CLsA
saveRDS(D_CLsA,"Scenario_1A_DistributionalClustering.RDS")



# Observational Clustering ------------------------------------------------
DF_Csi <- list()
for(i in 1:3){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S1_A[[i]]$Csi_ij,~.x)))
}

PSMs_2 <- map(DF_Csi, ~comp.psm(.x))
O_CLsA  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
saveRDS(O_CLsA,"Scenario_1A_ObservationalClustering.RDS")




###############################################################
# Case B
##################################################################################
# Data generation ---------------------------------------------------------
YD <- YG <- list()
set.seed(123456)
med <- c(0,5,10,13,16,20)
NUM_B <- c(5,10,20)
# CaseB0 <- function(m,s=sqrt(.6),N, nsim=100){
#   replicate(nsim, {x <- sample(1:N,1)
#   cbind(ifelse(x==1,1,0),
#         rnorm(1,m[x],s),x) },simplify = T)
# }
CaseB1 <- function(m,s=sqrt(.6), nsim=100){
  ALL <- sapply(1:6,function(xx) sapply(1:xx,
                    function(x)  rnorm(nsim,m[x],s)))
  return(unlist(ALL))
}


CaseB <- function(m,s=sqrt(.6),N){
  y1 <- c();
  for(ll in 1:2){
  y1 <- c(y1,CaseB1(med,nsim = N))
    }
   return(y1)
}

YD_B <- YG_B <- list()
for(i in 1:3){
  y1 <- CaseB(m = med,s = sqrt(.6),N = NUM_B[i])
  YD_B[[i]] <- c(y1)
  repet   <- c((1:6)*NUM_B[i],(1:6)*NUM_B[i]) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG_B[[i]] <- y_group
}

YD <- YD_B
YG <- YG_B

##################################################################################
plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])

#####################################g#############################################
# Settin up the model -----------------------------------------------------

S1_B <- mclapply(1:3,parallel.CAM,mc.cores = 3)
saveRDS(S1_B,"Scenario_1B_CAMoutput.RDS")

# Distributional Clustering -----------------------------------------------

DF_Z <- list()
for(i in 1:3){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S1_B[[i]]$Z_j,~.x)))
}

PSMs <- map(DF_Z, ~comp.psm(.x))
#D_CLs  <- map(PSMs, ~mcclust.ext::minVI(.x)$cl)
D_CLs_B  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)
D_CLs_B
saveRDS(D_CLs_B,"Scenario_1B_DistributionalClustering.RDS")



# Observational Clustering ------------------------------------------------
DF_Csi <- list()
for(i in 1:3){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S1_B[[i]]$Csi_ij,~.x)))
}

PSMs_2 <- map(DF_Csi, ~comp.psm(.x))
O_CLsB  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
saveRDS(O_CLsB,"Scenario_1B_ObservationalClustering.RDS")

plot(YD_B[[3]],col=O_CLsB[[3]])



# NestedDP ----------------------------------------------------------------

YD <- YG <- list()
set.seed(123456)
med <- c(0,5,10,13,16,20)
NUM_B <- c(5,10,20)

YD_B <- YG_B <- list()
for(i in 1:3){
  y1 <- CaseB(m = med,s = sqrt(.6),N = NUM_B[i])
  YD_B[[i]] <- c(y1)
  repet   <- c((1:6)*NUM_B[i],(1:6)*NUM_B[i]) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG_B[[i]] <- y_group
}

YD <- YD_B
YG <- YG_B

##################################################################################
plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])

Rcpp::sourceCpp("CommonAtomModel functions/NDP_via_EnhNDP/newEnDP_2018.cpp")
source("CommonAtomModel functions/NDP_via_EnhNDP/newEnDP_2018.R")

NSIM=50000;thinning=1; burn_in=50000; verbose.step=10; 
L=30; K=25; verbose=TRUE; monitor=T; 
m.step=100; fixedAB=F; 
restart=F; restartList=NA

prior <- list(
  # hyperparameters NIG
  m0=0, k0=1/var(YD[[3]]), a0=3, b0=1,
  # hyperparameters alpha and beta
  a_alpha=3, b_alpha = 3,
  a_beta =3, b_beta = 3, A=1, B=1 )


NDP_S1 <- newEnDP_Gibbs(y_obser = YD[[3]], y_group = YG[[3]],
                        K = K,L = L,
                        prior = prior,
                        NSIM = NSIM,
                        burn_in = burn_in,thinning = thinning,
                        verbose = T,verbose.step = 50,monitor = T,
                        m.step = 500,fixedAB = F,restart = F,
                        restartList = NULL,
                        RDG=T,cheap = T,seed = 12345)


saveRDS(NDP_S1,"Scenario1B_Nested_RDG.RDS")