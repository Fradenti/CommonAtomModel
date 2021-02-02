library(parallel)
library("rjags")
library(MCMCvis)
rjags::load.module("glm")

# All seeds specified under R.3.5.3 --> if updated R version is used, run first
# RNGkind(sample.kind = "Rounding")

##########################
parallel_CamJags_i_S1a_K <- function(i,Yall,path,K){
  y <- matrix(Yall[[i]][[K]],K*25,12) # obs on rows and groups on cols
  dataList <- list(y=y,J=ncol(y),n=nrow(y),
                   L=30,K=30,alpha=1,beta=1,
                   a_prior =3, b_prior=1,
                   m_prior = mean(Yall[[i]][[K]]), 
                   k_prior = 1/var(Yall[[i]][[K]]))
  model.fit_basic <- jags.model(path,
                                data=dataList,
                                n.chains = 1,
                                n.adapt = 5000)
  update(model.fit_basic, 5000)
  t1 <- Sys.time()
  out <- jags.samples(model.fit_basic,
                      c('Mij', 'zj',"theta_l"),
                      10000)
  t2 <- Sys.time()
  saveRDS(
    list(model=out,time=t2-t1),
    paste0("Competitors/Results_CAM/out_CAM_S1a_n",K*25,"_",i,".RDS"))
}
# -------------------------------------------------------------------------

path    <- "JAGS_files/CAM_alphabeta_fixed.jags"

# -------------------------------------------------------------------------

# Let us start with 30 datasets -------------------------------------------
# Scenario 1a
ALL_s1a <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Yall_s1a   <- ALL_s1a[[1]]
Gall_s1a   <- ALL_s1a[[2]]
Oall_s1a   <- ALL_s1a[[3]]
set.seed(12345)
mclapply(1:30,parallel_CamJags_i_S1a_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=1)
set.seed(679)
mclapply(1:30,parallel_CamJags_i_S1a_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=2)
set.seed(995859)
mclapply(1:30,parallel_CamJags_i_S1a_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=3)
####################################################################

# # Let us start with 30 datasets -------------------------------------------
# # Scenario 1b

# Here, you need to change the size of the dataset withing the parallel Camjags functions

# ALL_s1b    <- readRDS("Simulated_Data/ALL_S1B_100.RDS")
# Yall_s1b   <- ALL_s1b[[1]]
# Gall_s1b   <- ALL_s1b[[2]]
# Oall_s1b   <- ALL_s1b[[3]]

# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s1b,path=path,K=1)
# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s1b,path=path,K=2)
# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s1b,path=path,K=3)
####################################################################












library(parallel)
library("rjags")
library(MCMCvis)
rjags::load.module("glm")




# Let us start with 30 datasets -------------------------------------------
# Scenario 2
ALL_s2 <- readRDS("Simulated_Data/ALL_S2_100.RDS")
Yall_s2   <- ALL_s2[[1]]
Gall_s2   <- ALL_s2[[2]]
Oall_s2   <- ALL_s2[[3]]

# -------------------------------------------------------------------------

path    <- "JAGS_files/CAM_alphabeta_fixed.jags"

# -------------------------------------------------------------------------

##########################
parallel_CamJags_i_S2_K <- function(i,Yall,path,K){
  J <- c(2,4,6)
  y <- matrix(Yall[[i]][[J[K]]],40,4*J[K]) # obs on rows and groups on cols
  dataList <- list(y=y,J=ncol(y),n=nrow(y),
                   L=30,K=30,alpha=1,beta=1,
                   a_prior =3, b_prior=1,
                   m_prior = mean(Yall[[i]][[J[K]]]), 
                   k_prior = 1/var(Yall[[i]][[J[K]]]))
  model.fit_basic <- jags.model(path,
                                data=dataList,
                                n.chains = 1,
                                n.adapt = 5000)
  update(model.fit_basic, 5000)
  t1 <- Sys.time()
  out <- jags.samples(model.fit_basic,
                      c('Mij', 'zj'),
                      10000)
  t2 <- Sys.time()
  saveRDS(
    list(model=out,time=t2-t1),
    paste0("Competitors/Results_CAM/out_CAM_S2_J",J[K],"_",i,".RDS"))
}

set.seed(12345)
mclapply(1:30,parallel_CamJags_i_S2_K,
         mc.cores = 15,Yall=Yall_s2,path=path,K=1)
set.seed(54321)
mclapply(1:30,parallel_CamJags_i_S2_K,
         mc.cores = 15,Yall=Yall_s2,path=path,K=2)
set.seed(5555)
mclapply(1:30,parallel_CamJags_i_S2_K,
         mc.cores = 15,Yall=Yall_s2,path=path,K=3)
####################################################################

# Let us start with 30 datasets -------------------------------------------
# Scenario 3
# ALL_s3 <- readRDS("Simulated_Data/ALL_S3_100.RDS")
# Yall_s3   <- ALL_s3[[1]]
# Gall_s3   <- ALL_s3[[2]]
# Oall_s3   <- ALL_s3[[3]]
# 
# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s3,path=path,K=1)
# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s3,path=path,K=2)
# mclapply(1:30,parallel_CamJags_i_K,
#          mc.cores = 30,Yall=Yall_s3,path=path,K=3)
####################################################################