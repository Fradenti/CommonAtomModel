
# Collect PSMs Jags S1 ----------------------------------------------------

Rcpp::sourceCpp("Functions/CAM/newCAM.cpp")

model <- c("CAM","NDP","hNDP")

K=3; i =1
PSMs <- list()
for(i in 1:30){
  psm <- list()  
   for(k in 1:3){
   path     <-   paste0("Competitors/Results_",model[k],"/out_",model[k],"_S1a_n",K*25,"_",i,".RDS")
   obj      <-   readRDS(path)
   Z        <- t(obj$model$zj[,,1])
   psm[[k]] <- PSM(Z)
   }
  PSMs[[i]] <- psm
}

saveRDS(PSMs,"Competitors/Results_Jags_Psms_S1a_K3.RDS")


# plotting ----------------------------------------------------------------

library(tidyverse)
library(patchwork)
K1 <- Results_Jags_Psms_S1a_K3[[1]]

rr <- rbind(rep(1:6,2),rep(1:6,2))

vv0 <- reshape2::melt(PSM(rr))
PSM0 <- ggplot(data=as_tibble(vv0))+geom_tile(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Unit")+ylab("Unit")+scale_fill_viridis_c(option = "B")+
  ggtitle("Posterior Coclustering Matrix - Ground Truth") + theme(legend.position = "none")+
  theme(axis.title = element_text(size=15),title = element_text(size=15))
PSM0

vv1 <- reshape2::melt(K1[[1]])
PSM1 <- ggplot(data=as_tibble(vv1))+geom_tile(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Unit")+ylab("Unit")+scale_fill_viridis_c(option = "B")+
  ggtitle("Posterior Coclustering Matrix - CAM") + theme(legend.position = "none")+
  theme(axis.title = element_text(size=15),title = element_text(size=15))
PSM1
vv2 <- reshape2::melt(K1[[2]])
PSM2 <- ggplot(data=as_tibble(vv2))+geom_tile(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Unit")+ylab("Unit")+scale_fill_viridis_c("Coclustering Probability",option = "B")+
  ggtitle("Posterior Coclustering Matrix - NDP") + theme(legend.position = "bottom")+
  theme(axis.title = element_text(size=15),title = element_text(size=15))
PSM2
vv3 <- reshape2::melt(K1[[3]])
PSM3 <- ggplot(data=as_tibble(vv3))+geom_tile(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Unit")+ylab("Unit")+scale_fill_viridis_c(option = "B")+
  ggtitle("Posterior Coclustering Matrix - nHDP") + theme(legend.position = "none")+
  theme(axis.title = element_text(size=15),title = element_text(size=15))
PSM3


(PSM0+PSM1)/(PSM2+PSM3)
ggsave("Competitors/PSMS.pdf",width = 15,height = 8)
ggsave("Competitors/PSMS.png",width = 15,height = 8)

DD <- rbind(
tibble(vv0,name="Ground Truth"),
tibble(vv1,name="CAM"),
tibble(vv2,name="NDP"),
tibble(vv3,name="nHDP")
)
DD <- DD %>% mutate(name = factor(name,levels = c("Ground Truth","CAM","NDP","nHDP")))
PSMALL <- ggplot(data=DD)+geom_tile(aes(x=Var1, y=Var2, fill=value))+
  facet_wrap(~name)+
  theme_bw()+xlab("Unit")+ylab("Unit")+
  scale_fill_viridis_c("Posterior Coclustering Probability",option = "B")+
  theme(axis.title = element_text(size=15),title = element_text(size=15),
        strip.text = element_text(size=12))+
  theme(legend.position = "bottom")+
  ggtitle("DC - Dataset 1 - S 1A - Conf. 3")
PSMALL
ggsave("Competitors/PSMS.pdf",width = 12,height = 12)
ggsave("Competitors/PSMS.png",width = 12,height = 12)
