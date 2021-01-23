library(parallel)
library("rjags")
library(MCMCvis)
rjags::load.module("glm")

##########################
parallel_CamJags_i_K <- function(i,Yall,path,K){
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
    paste0("Competitors/Results_CAM/out_CAM_S1A_n",K*25,"_",i,".RDS"))
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

mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=1)
mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=2)
mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s1a,path=path,K=3)
####################################################################

# # Let us start with 30 datasets -------------------------------------------
# # Scenario 1b
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


# Let us start with 30 datasets -------------------------------------------
# Scenario 2
ALL_s2 <- readRDS("Simulated_Data/ALL_S2_100.RDS")
Yall_s2   <- ALL_s2[[1]]
Gall_s2   <- ALL_s2[[2]]
Oall_s2   <- ALL_s2[[3]]

mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s2,path=path,K=1)
mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s2,path=path,K=2)
mclapply(1:30,parallel_CamJags_i_K,
         mc.cores = 30,Yall=Yall_s2,path=path,K=3)
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