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
  # CAM --------------------------------------------------------------------
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

ALL <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Yall   <- ALL[[1]]
Gall   <- ALL[[2]]
Oall   <- ALL[[3]]


mclapply(1:30,parallel_CamJags_i_K,mc.cores = 30,Yall=Yall,path=path,K=1)
mclapply(1:30,parallel_CamJags_i_K,mc.cores = 30,Yall=Yall,path=path,K=2)
mclapply(1:30,parallel_CamJags_i_K,mc.cores = 30,Yall=Yall,path=path,K=3)
####################################################################





















y <- matrix(YD[[2]],50,12)
dataList <- list(y=y,J=12,n=50,L=30,K=30,R=20)
# CAM --------------------------------------------------------------------
model.fit_basic <- jags.model(path,
                              data=dataList,
                              n.chains = 1,
                              n.adapt = 1000)
update(model.fit_basic, 1000)
t1 <- Sys.time()
out <- jags.samples(model.fit_basic,
                    c('Mij', 'zj',"theta_l","alpha","beta"),
                    2000)
t2 <- Sys.time()
saveRDS(list(out,time=t2-t1),"out_v2_CAM_S1A_n50.RDS")




y <- matrix(YD[[3]],75,12)
dataList <- list(y=y,J=12,n=75,L=30,K=30,R=20)
# CAM --------------------------------------------------------------------
model.fit_basic <- jags.model(path,
                              data=dataList,
                              n.chains = 1,
                              n.adapt = 5000)
update(model.fit_basic, 5000)
t1 <- Sys.time()
out <- jags.samples(model.fit_basic,
                    c('Mij', 'zj',"theta_lk","alpha","beta"),
                    10000)
t2 <- Sys.time()
saveRDS(list(out,time=t2-t1),"out_v2_CAM_S1A_n75.RDS")
