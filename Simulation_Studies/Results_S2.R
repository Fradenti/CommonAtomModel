library(tidyverse)
# Results from simulation studies -----------------------------------------
S2      <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_2_CAMoutput.RDS")

# Ground Truth
# OC
NNN <- 40
lab <- c(rep(1,NNN*.75),  rep(2,NNN*.25),
         rep(1,NNN*.25),   rep(2,NNN*.75),
         rep(1,NNN*.33),  rep(3,NNN*.34),rep(4,NNN*.33),
         rep(1,NNN*.25+1),rep(4,NNN*.25),
         rep(3,NNN*.25),  rep(5,NNN*.25))

plot(S2[[1]]$y_obser,col=lab)
GT_O <- list()
GT_O[[1]] <- lab
GT_O[[2]] <- c(lab,ttt1)
GT_O[[3]] <- c(lab,ttt2)
GT_O[[4]] <- c(lab,ttt3)
GT_O[[5]] <- c(lab,ttt4)
GT_O[[6]] <- c(lab,ttt5)
plot(S2[[4]]$y_obser,col=GT_O[[4]])

gt_distr <- list()
for(i in 1:6){
  gt_distr[[i]] <- rep(1:4,i) 
}

#############################################################
# Load results
DC_S2 <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_2_DistributionalClustering.RDS")
OC_S2 <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_2_ObservationalClustering.RDS")

# Adjuster RI
unlist(map(DC_S2, ~length(unique(.x))))
round(unlist(map2(DC_S2,gt_distr,~ mcclust::arandi(.x,.y))),3)
map(OC_S2, ~length(unique(.x)))
round(unlist(map2(OC_S2,GT_O,mcclust::arandi)),3)

# Compute GT PSM
DC_GT_PSMs <- map(gt_distr,~ mcclust::comp.psm(rbind(.x,.x)))
OC_GT_PSMs <- map(GT_O, ~ mcclust::comp.psm(rbind(.x,.x)))


# Estimate PSMs
DF_Z <- list()
for(i in 1:6){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S2[[i]]$Z_j,~.x)))
}
DC_EST_PSMs <- map(DF_Z, ~mcclust::comp.psm(.x))

DF_Csi <- list()
for(i in 1:6){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S2[[i]]$Csi_ij,~.x)))
}

OC_EST_PSMs <- map(DF_Csi, ~mcclust::comp.psm(.x))



# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- map2(DC_EST_PSMs,DC_GT_PSMs,~StatPerMeCo::Frobenius(.x,.y))
UnnDist_OC <- map2(OC_EST_PSMs,OC_GT_PSMs,~StatPerMeCo::Frobenius(.x,.y))


# Normalized Frobenious (pairwise) distances
max_DC <- map(DC_EST_PSMs,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))
max_OC <- map(OC_EST_PSMs,
              ~ StatPerMeCo::Frobenius(
                matrix(0,nrow(.x),ncol(.x)),
                matrix(1,nrow(.x),ncol(.x))))

round(unlist(map2(UnnDist_DC,max_DC,~.x/.y)),3)
round(unlist(map2(UnnDist_OC,max_OC,~.x/.y)),3)






S2_NDP <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario2_NDP.RDS")

DC_EST_PSM <- mcclust::comp.psm(t(S2_NDP$Z_j))
OC_EST_PSM <- mcclust::comp.psm(t(S2_NDP$Csi_ij))

DC_NDP_S2 <- mcclust.ext::minVI(DC_EST_PSM,method = "greedy")$cl
OC_NDP_S2 <- mcclust.ext::minVI(OC_EST_PSM)$cl

# Adjuster RI
length(unique(DC_NDP_S2))
mcclust::arandi(gt_distr[[6]],DC_NDP_S2)
length(unique(OC_NDP_S2))
mcclust::arandi(GT_O[[6]],OC_NDP_S2)



# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- StatPerMeCo::Frobenius(DC_GT_PSMs[[6]],DC_EST_PSM)
UnnDist_OC <- StatPerMeCo::Frobenius(OC_GT_PSMs[[6]],OC_EST_PSM)


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


