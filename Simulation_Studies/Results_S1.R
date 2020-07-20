library(tidyverse)
# Results from simulation studies -----------------------------------------


# Simulation S1 -----------------------------------------------------------

# Weird thing:
# GT <- list()
# # Data
# YD <- YG <- list()
# set.seed(123456)
# med <- c(0,5,10,13,16,20)
# CaseA <- function(m,s=sqrt(.6),N){
#   p <- c()
#   xx <- c()
#   for(i in 1:6){
#     x <- sample(1:i,N,T)
#     xx <- c(xx,x)
#     p <- c(p,rnorm(N,m[x],s)) }
#   return(cbind(p,xx))
# }
# 
# NUM_A <- c(25,50,75)
# for(i in 1:3){
#   y1 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
#   y2 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
#   YD[[i]] <- c(y1[,1],y2[,1])
#   repet   <- rep(NUM_A[i],12) #rep(200,6)
#   y_group <- c(rep(1:length(repet),repet))
#   YG[[i]] <- y_group
#   GT[[i]] <- c(y1[,2],y2[,2])
# }
# 

# cosa strana: su server e su mio pc stesso seme porta a risultati diversi :(
GT_O     <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/GT_ObseCluster_S1A.RDS")
S1A      <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1A_CAMoutput.RDS")

plot(S1A[[1]]$y_obser,col=GT_O[[1]])
gt_distr <- rep(1:6,2)
plot(YD[[1]],col=YG[[1]])

#############################################################
# Load results
DC_S1A <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1A_DistributionalClustering.RDS")
OC_S1A <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1A_ObservationalClustering.RDS")

# Adjuster RI
map(DC_S1A, ~length(unique(.x)))
round(unlist(map(DC_S1A,~ mcclust::arandi(.x,gt_distr))),3)
map(OC_S1A, ~length(unique(.x)))
round(unlist(map2(OC_S1A,GT_O,mcclust::arandi)),3)

# Compute GT PSM
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
OC_GT_PSMs <- map(GT_O, ~ mcclust::comp.psm(rbind(.x,.x)))


# Estimate PSMs
DF_Z <- list()
for(i in 1:3){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S1A[[i]]$Z_j,~.x)))
}
DC_EST_PSMs <- map(DF_Z, ~mcclust::comp.psm(.x))

DF_Csi <- list()
for(i in 1:3){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S1A[[i]]$Csi_ij,~.x)))
}

OC_EST_PSMs <- map(DF_Csi, ~mcclust::comp.psm(.x))



# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- map(DC_EST_PSMs,~StatPerMeCo::Frobenius(.x,DC_GT_PSM))
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








# Simulation S1B -----------------------------------------------------------

# Create indexes for Observational GT
GT <- list()
gt <- c()
NUM_B <- c(5,10,20)
for(i in 1:6) gt <- c(gt,1:i)
gt <- c(gt,gt)
for(l in 1:3){
GT[[l]] <- rep(gt,rep(NUM_B[l],7*6))
}
GT_O <- GT
gt_distr <- rep(1:6,2)

#############################################################
# Load results
S1B    <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1B_CAMoutput.RDS")
DC_S1B <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1B_DistributionalClustering.RDS")
OC_S1B <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario_1B_ObservationalClustering.RDS")
# double check
# all(S1B[[1]]$y_obser==YD[[1]])
# all(S1B[[2]]$y_obser==YD[[2]])
# all(S1B[[3]]$y_obser==YD[[3]])



# Adjuster RI
map(DC_S1B, ~length(unique(.x)))
round(unlist(map(DC_S1B,~ mcclust::arandi(.x,gt_distr))),3)
map(OC_S1B, ~length(unique(.x)))
round(unlist(map2(OC_S1B,GT_O,mcclust::arandi)),3)

# Compute GT PSM
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
OC_GT_PSMs <- map(GT_O, ~ mcclust::comp.psm(rbind(.x,.x)))


# Estimate PSMs
DF_Z <- list()
for(i in 1:3){
  DF_Z[[i]] <- t(as.matrix(map_dfc(S1B[[i]]$Z_j,~.x)))
}
DC_EST_PSMs <- map(DF_Z, ~mcclust::comp.psm(.x))

DF_Csi <- list()
for(i in 1:3){
  DF_Csi[[i]] <- t(as.matrix(map_dfc(S1B[[i]]$Csi_ij,~.x)))
}

OC_EST_PSMs <- map(DF_Csi, ~mcclust::comp.psm(.x))



# Unnormalized Frobenious (pairwise) distances
UnnDist_DC <- map(DC_EST_PSMs,~StatPerMeCo::Frobenius(.x,DC_GT_PSM))
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




S1B_NDP <- readRDS("/home/fra/SIMULAZIONI per articoli/CAM_JUL20/Scenario1B_Nested_RDG.RDS")

DC_EST_PSM <- mcclust::comp.psm(t(apply(S1B_NDP$Z_j,2,reset)))
OC_EST_PSM <- mcclust::comp.psm(t(S1B_NDP$Csi_ij))

DC_NDP_S1B <- mcclust.ext::minVI(DC_EST_PSM,method = "greedy")$cl
OC_NDP_S1B <- mcclust.ext::minVI(OC_EST_PSM)$cl

# Adjuster RI
length(unique(DC_NDP_S1B))
mcclust::arandi(gt_distr,DC_NDP_S1B)
length(unique(OC_NDP_S2))
mcclust::arandi(GT_O[[3]],OC_NDP_S1B)

# Compute GT PSM
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
OC_GT_PSMs <- map(GT_O, ~ mcclust::comp.psm(rbind(.x,.x)))


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




