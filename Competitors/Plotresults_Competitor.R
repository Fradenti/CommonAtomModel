
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
                    Title = "ARI - DC - Jags - S 1A",yl = "Adjusted Rand Index")

nfd <- plotter_jags(tj1a_c[[3]],tj1a_n[[3]],tj1a_h[[3]],log=F,
                    Title = "NFD - DC - Jags - S 1A",yl = "Normalized Frobenious Distance")

ari+nfd

ggsave("Competitors/Jags_ARINFD_DC_jags_S1a.png",width = 12,height = 6)
ggsave("Competitors/Jags_ARINFD_DC_jags_S1a.pdf",width = 12,height = 6)




ari2 <- plotter_jags(x1 = tj2_c[[2]],
                    x2 = tj2_n[[2]],
                    x3 = tj2_h[[2]],
                    log=F,
                    Title = "ARI - DC - Jags - S 2",yl = "Adjusted Rand Index")

nfd2 <- plotter_jags(tj2_c[[3]],
                    tj2_n[[3]],
                    tj2_h[[3]],log=F,Title = "NFD - DC - Jags - S 2",yl = "Normalized Frobenious Distance")

ari2+nfd2

ggsave("Competitors/Jags_ARINFD_DC_jags_S2.png",width = 12,height = 6)
ggsave("Competitors/Jags_ARINFD_DC_jags_S2.pdf",width = 12,height = 6)




(ari+ari2)/(nfd+nfd2)

ggsave("Competitors/All_ARINFD_DC_jags_S2.png",width = 12,height = 8)
ggsave("Competitors/All_ARINFD_DC_jags_S2.pdf",width = 12,height = 8)




# tabluar report ----------------------------------------------------------

c1 <- map(tj1a_c, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  )) 
n1 <- map(tj1a_n, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  )) 
h1 <- map(tj1a_h, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  )) 
c2 <- map( tj2_c, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  ))  
n2 <- map( tj2_n, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  ))  
h2 <- map( tj2_h, ~ apply(.x,2,function(x)  c(mean(x),sd(x))  ))  


ARI <- rbind(cbind(c1$ARI,
                   n1$ARI,
                   h1$ARI),
                   cbind(c2$ARI,
                   n2$ARI,
                   h2$ARI))

NFD <- rbind(cbind(c1$NFD,
                   n1$NFD,
                   h1$NFD),
             cbind(c2$NFD,
                   n2$NFD,
                   h2$NFD))

round(rbind(ARI,NFD),4)
X <- as_tibble(t(round(rbind(ARI,NFD),4))) %>% 
  mutate(V2 = paste0("(",V2,")"),
         V4 = paste0("(",V4,")"),
         V6 = paste0("(",V6,")"),
         V8 = paste0("(",V8,")"))
t(X)


CLU <- rbind(cbind(c1$numclust,
                    n1$numclust,
                    h1$numclust),
              cbind(c2$numclust,
                    n2$numclust,
                    h2$numclust))
CLU




c1 <- apply(tj1a_c[[1]],2,function(x)  names(table(x))[which.max(table(x))]  ) 
n1 <- apply(tj1a_n[[1]],2,function(x)  names(table(x))[which.max(table(x))]  ) 
h1 <- apply(tj1a_h[[1]],2,function(x)  names(table(x))[which.max(table(x))]  ) 
c2 <- apply( tj2_c[[1]],2,function(x)  names(table(x))[which.max(table(x))]  )  
n2 <- apply( tj2_n[[1]],2,function(x)  names(table(x))[which.max(table(x))]  )  
h2 <- apply( tj2_h[[1]],2,function(x)  names(table(x))[which.max(table(x))]  )  
