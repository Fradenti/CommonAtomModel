library(tidyverse)
library(patchwork)
library(reshape2)

# ARI and rest Analysis ---------------------------------------------------


mcclust::arandi(1:4,c(1,2,2,1))

t1a <- readRDS("Simulation_Study_Slice/Results_30_1A_Distributional.RDS")
t1b <- readRDS("Simulation_Study_Slice/Results_30_1B.RDS")
t2  <- readRDS("Simulation_Study_Slice/Results_30_2.RDS")
t3  <- readRDS("Simulation_Study_Slice/Results_30_3.RDS")

# this has Nan in the first column, 
mcclust::arandi(1:10,1:10) #!! it happens when there is no double value in the data
mcclust::arandi(1:10,10:1) #!! so they are essentially 1

t2[[2]][is.na(t2[[2]][,1]),1] <- 1

ari <- plotter(t1a[[2]],t1b[[2]],t2[[2]],t3[[2]],log=F,Title = "ARI - DC - Slice Sampler",yl = "Adjusted Rand Index")

nfd <- plotter(t1a[[3]],t1b[[3]],t2[[3]],t3[[3]],log=F,Title = "NFD - DC - Slice Sampler",yl = "Normalized Frobenious Distance")

ari/nfd

ggsave("Simulation_Study_Slice/ARINFD_DC.png",width = 12,height = 10)
ggsave("Simulation_Study_Slice/ARINFD_DC.pdf",width = 12,height = 10)


ncl <- plotter(t1a[[1]],t1b[[1]],t2[[1]],t3[[1]],log=F,Title = "NFD - DC - Slice Sampler",yl = "Normalized Frobenious Distance")
# what to do with this guy?


# Observational clusters --------------------------------------------------

t1a <- readRDS("Simulation_Study_Slice/Results_30_1A_Observational.RDS")
t1b <- readRDS("Simulation_Study_Slice/Results_30_1B_Observational.RDS")
t2  <- readRDS("Simulation_Study_Slice/Results_30_2_Observational.RDS")
t3  <- readRDS("Simulation_Study_Slice/Results_30_3_Observational.RDS")


ari <- plotter(t1a[[2]],t1b[[2]],t2[[2]],t3[[2]],log=F,Title = "ARI - OC - Slice Sampler",yl = "Adjusted Rand Index")

nfd <- plotter(t1a[[3]],t1b[[3]],t2[[3]],t3[[3]],log=F,Title = "NFD - OC - Slice Sampler",yl = "Normalized Frobenious Distance")

ari/nfd

ggsave("Simulation_Study_Slice/ARINFD_OC.png",width = 12,height = 10)
ggsave("Simulation_Study_Slice/ARINFD_OC.pdf",width = 12,height = 10)


