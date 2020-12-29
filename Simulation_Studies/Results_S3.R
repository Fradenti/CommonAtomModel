library(tidyverse)
# Results from simulation studies -----------------------------------------

# Simulation S3 -----------------------------------------------------------
# Read output and extract results
S3 <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM/CAM_JUL20/Scenario3_DCAMoutput.RDS")

library(coda)
MC <- as.mcmc(cbind(S3[[1]]$A_DP,S3[[1]]$B_DP,S3[[2]]$A_DP,
                    S3[[2]]$B_DP,S3[[3]]$A_DP,S3[[3]]$B_DP))
# thinning to ease computation
MC <- as.mcmc(MC[seq(1,50000,by = 20),])
heidel.diag(MC)
nclust <- unlist(purrr::map(S2[[1]]$Z_j, ~length(unique(.x))))
heidel.diag(nclust)


R <- GT_O <- list()
for(i in 1:6){
  L <- R[[i]] <- S3[[i]]$y_obser
  GT_O[[i]] <- ifelse(L==0 | L==1,1, 
                      ifelse(L<=10,2,
                             ifelse(L<=50,3,4)  ))
}
plot(R[[6]],col=GT_O[[6]])
gt_distr <- c(1,2,3,2,1,2,1,3,2,1)
#############################################################
# Distributional Clusters
DF_Z <- list()
for(i in 1:6){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S3[[i]]$Z_j,~.x)))
}
# Estimate PSMs
PSMs <- map(DF_Z, ~mcclust::comp.psm(.x))
DC_S3  <- map(PSMs, ~mcclust.ext::minVI(.x,method = "greedy")$cl)


# Observational clusters
DF_Csi <- list()
for(i in 1:6){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S3[[i]]$Csi_ij,~.x)))
}
PSMs_2 <- map(DF_Csi, ~mcclust::comp.psm(.x))
saveRDS(PSMs_2,"/home/fra/SIMULAZIONI per articoli/CAM_JUL20/ObservationalPSMs_S3.RDS")
OC_S3  <- map(PSMs_2, ~mcclust.ext::minVI(.x)$cl)
saveRDS(O_CLs,"/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_3_ObservationalClustering.RDS")

# Compute perfomance indexes

# Adjuster RI
map(DC_S3, ~length(unique(.x)))
round(unlist(map(DC_S3,~ mcclust::arandi(.x,gt_distr))),3)
map(OC_S3, ~length(unique(.x)))
round(unlist(map2(OC_S3,GT_O,mcclust::arandi)),3)

# Compute GT PSM
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
OC_GT_PSMs <- map(GT_O, ~ mcclust::comp.psm(rbind(.x,.x)))


# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- map(PSMs,~StatPerMeCo::Frobenius(.x,DC_GT_PSM))
UnnDist_OC <- map2(PSMs_2,OC_GT_PSMs,~StatPerMeCo::Frobenius(.x,.y))


# Normalized Frobenious (pairwise) distances
max_DC <- map(PSMs,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))
max_OC <- map(PSMs_2,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))

round(unlist(map2(UnnDist_DC,max_DC,~.x/.y)),3)
round(unlist(map2(UnnDist_OC,max_OC,~.x/.y)),3)






# Nested DP ------------------------------------------------------------------
NDP_S3 <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario3_Nested_RDG.RDS")

# Partitions
L <-  NDP_S3$y_obser
GT_O <- ifelse(L==0 | L==1,1, 
                      ifelse(L<=10,2,
                             ifelse(L<=50,3,4)  ))
DC_EST_PSM <- mcclust::comp.psm(t(apply(NDP_S3$Z_j,2,reset)))
OC_EST_PSM <- mcclust::comp.psm(t(NDP_S3$Csi_ij))

DC_NDP_S3 <- mcclust.ext::minVI(DC_EST_PSM,method = "greedy")$cl
OC_NDP_S3 <- mcclust.ext::minVI(OC_EST_PSM)$cl

# Adjuster RI
length(unique(DC_NDP_S3))
mcclust::arandi(gt_distr,DC_NDP_S3)
length(unique(OC_NDP_S3))
mcclust::arandi(GT_O,OC_NDP_S3)

# Compute GT PSM
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
OC_GT_PSMs <- mcclust::comp.psm(rbind(GT_O,GT_O))


# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- StatPerMeCo::Frobenius(DC_GT_PSM,DC_EST_PSM)
UnnDist_OC <- StatPerMeCo::Frobenius(OC_GT_PSMs[[3]],OC_EST_PSM)


# Normalized Frobenious (pairwise) distances
DC_EST_PSM_NDP <- list()
DC_EST_PSM_NDP[[1]] <- DC_EST_PSM
OC_EST_PSM_NDP <- list()
OC_EST_PSM_NDP[[1]] <- OC_EST_PSM
max_DC <- map(DC_EST_PSM_NDP,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))
max_OC <- map(OC_EST_PSM_NDP,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))
(unlist(map2(UnnDist_DC,max_DC,~.x/.y)))
round(unlist(map2(UnnDist_OC,max_OC,~.x/.y)),3)




