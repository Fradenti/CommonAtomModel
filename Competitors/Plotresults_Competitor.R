
library(tidyverse)
library(patchwork)
library(reshape2)




tj1a_c <- readRDS("Competitors/Results_CAM_S1a.RDS")
tj1a_n <- readRDS("Competitors/Results_NDP_S1a.RDS")
tj1a_h <- readRDS("Competitors/Results_hNDP_S1a.RDS")
tj2_c  <- readRDS("Competitors/Results_CAM_S2.RDS")
tj2_n  <- readRDS("Competitors/Results_NDP_S2.RDS")
tj2_h  <- readRDS("Competitors/Results_hNDP_S2.RDS")


# this has Nan in the first column, 
mcclust::arandi(1:10,1:10) #!! it happens when there is no double value in the data
mcclust::arandi(1:10,10:1) #!! so they are essentially 1

tj2_c[[2]]

ari <- plotter_jags(x1 = tj1a_c[[2]],
                    x2 = tj1a_n[[2]],
                    x3 = tj1a_h[[2]],
                    log=F,
                    Title = "ARI - DC - Jags - S 1a",yl = "Adjusted Rand Index")

nfd <- plotter_jags(tj1a_c[[3]],tj1a_n[[3]],tj1a_h[[3]],log=F,Title = "NFD - DC - Jags",yl = "Normalized Frobenious Distance")

ari/nfd

ggsave("Competitors/ARINFD_DC_jags_S1a.png",width = 12,height = 10)
ggsave("Competitors/ARINFD_DC_jags_S1a.pdf",width = 12,height = 10)




ari <- plotter_jags(x1 = tj2_c[[2]],
                    x2 = tj2_n[[2]],
                    x3 = tj2_h[[2]],
                    log=F,
                    Title = "ARI - DC - Jags - S 2",yl = "Adjusted Rand Index")

nfd <- plotter_jags(tj2_c[[3]],
                    tj2_n[[3]],
                    tj2_h[[3]],log=F,Title = "NFD - DC - Jags",yl = "Normalized Frobenious Distance")

ari/nfd

ggsave("Competitors/ARINFD_DC_jags_S2.png",width = 12,height = 10)
ggsave("Competitors/ARINFD_DC_jags_S2.pdf",width = 12,height = 10)
