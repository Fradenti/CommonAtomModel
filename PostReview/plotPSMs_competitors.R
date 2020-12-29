library(mcclust)

reset <- function(x){
  z  <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}

# Scenario 1 A -------------------------------------------------------------------------
gt_distr   <- rep(1:6,2)
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))


# Load --------------------------------------------------------------------

DC_memberships_outCAM   = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_CAM_S1a.RDS")$zj
DC_memberships_outNest  = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_NDP_S1a.RDS")[[1]]$zj
DC_memberships_outHNest = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_hNDP_S1a.RDS")[[1]]$zj
DC_memberships_outSemi  = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_semiHDP_S1a.RDS")[[1]]$rho



# PSM ---------------------------------------------------------------------
A0 <- DC_GT_PSM
A1 <- comp.psm(t(apply(DC_memberships_outCAM  [,,1],2,reset)))
A2 <- comp.psm(t(apply(DC_memberships_outNest [,,1],2,reset))) 
A3 <- comp.psm(t(apply(DC_memberships_outHNest[,,1],2,reset)))
A4 <- comp.psm(t(apply(DC_memberships_outSemi [,,1],2,reset))) 

library(reshape2)
library(tidyverse)
DD <- rbind(cbind(melt(A0),lab=1),
            cbind(melt(A1),lab=2),
            cbind(melt(A2),lab=3),
            cbind(melt(A3),lab=4),
            cbind(melt(A4),lab=5))
DD <- as_tibble(DD)
DD <- DD %>% mutate(lab = case_when(lab==1~"0 - DC - Ground Truth", 
                                    lab==2~"1 - DC - CAM",
                                    lab==3~"2 - DC - Nested DP",
                                    lab==4~"3 - DC - Hierarchical Nested DP",
                                    lab==5~"4 - DC - Semi Hierarchical DP"))



ggplot(data=DD)+geom_tile(aes(x=Var1,y=Var2,fill=value))+theme_bw()+
  scale_fill_viridis_c("Coclustering\nProbability",option = "B",end = .95)+xlab("")+ylab("")+
        facet_wrap(~lab,nrow = 3)+theme(legend.position=c(.75,.15),
                                        strip.text.x = element_text(size = 12))+
    scale_x_discrete(limits=as.character(1:12))+
    scale_y_discrete(limits=as.character(1:12))
ggsave("../CAM/Analysis Review 1/figures/PSMs_S1a.png",height = 10,width = 8)
ggsave("../CAM/Analysis Review 1/figures/PSMs_S1a.pdf",height = 10,width = 8)


# no mario

ggplot(data=DD %>% filter(lab!="4 - DC - Semi Hierarchical DP"))+geom_tile(aes(x=Var1,y=Var2,fill=value))+theme_bw()+
  scale_fill_viridis_c("Coclustering\nProbability",option = "B",end = .95)+xlab("")+ylab("")+
  facet_wrap(~lab,nrow = 3)+theme(legend.position="bottom",
                                  strip.text.x = element_text(size = 12))+
  scale_x_discrete(limits=as.character(1:12))+
  scale_y_discrete(limits=as.character(1:12))
ggsave("../CAM/Analysis Review 1/figures/PSMs_S1a_noMario.png",height = 10,width = 8)
ggsave("../CAM/Analysis Review 1/figures/PSMs_S1a_noMario.pdf",height = 10,width = 8)







# Scenario 2 -------------------------------------------------------------------------
gt_distr   <- rep(1:4,6)
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))


# Load --------------------------------------------------------------------
DC_memberships_outCAM   = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_CAM_S2.RDS")[[1]]$zj
DC_memberships_outNest  = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_NDP_S2.RDS")[[1]]$zj
DC_memberships_outHNest = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_hNDP_S2.RDS")[[1]]$zj
DC_memberships_outSemi  = readRDS("~/SIMULAZIONI per articoli/CAM/REV_JASA/RESULTS/out_semiHDP_S2.RDS")[[1]]$rho



# PSM ---------------------------------------------------------------------
A0 <- DC_GT_PSM
A1 <- comp.psm(t(apply(DC_memberships_outCAM  [,,1],2,reset)))
A2 <- comp.psm(t(apply(DC_memberships_outNest [,,1],2,reset))) 
A3 <- comp.psm(t(apply(DC_memberships_outHNest[,,1],2,reset)))
A4 <- comp.psm(t(apply(DC_memberships_outSemi [,,1],2,reset))) 

library(reshape2)
library(tidyverse)
DD <- rbind(cbind(melt(A0),lab=1),
            cbind(melt(A1),lab=2),
            cbind(melt(A2),lab=3),
            cbind(melt(A3),lab=4),
            cbind(melt(A4),lab=5))
DD <- as_tibble(DD)
DD <- DD %>% mutate(lab = case_when(lab==1~"0 - DC - Ground Truth", 
                                    lab==2~"1 - DC - CAM",
                                    lab==3~"2 - DC - Nested DP",
                                    lab==4~"3 - DC - Hierarchical Nested DP",
                                    lab==5~"4 - DC - Semi Hierarchical DP"))



ggplot(data=DD)+geom_tile(aes(x=Var1,y=Var2,fill=value))+theme_bw()+
  scale_fill_viridis_c("Coclustering\nProbability",option = "B",end = .95)+xlab("")+ylab("")+
  facet_wrap(~lab,nrow = 3)+theme(legend.position=c(.75,.15),
                                  strip.text.x = element_text(size = 12))+
  scale_x_discrete(limits=as.character(1:24))+
  scale_y_discrete(limits=as.character(1:24))
ggsave("../CAM/Analysis Review 1/figures/PSMs_S2.png",height = 10,width = 8)
ggsave("../CAM/Analysis Review 1/figures/PSMs_S2.pdf",height = 10,width = 8)

# no mario

ggplot(data=DD %>% filter(lab!="4 - DC - Semi Hierarchical DP"))+geom_tile(aes(x=Var1,y=Var2,fill=value))+theme_bw()+
  scale_fill_viridis_c("Coclustering\nProbability",option = "B",end = .95)+xlab("")+ylab("")+
  facet_wrap(~lab,nrow = 3)+theme(legend.position="bottom",
                                  strip.text.x = element_text(size = 12))+
  scale_x_discrete(limits=as.character(1:24))+
  scale_y_discrete(limits=as.character(1:24))
ggsave("../CAM/Analysis Review 1/figures/PSMs_S2_noMario.png",height = 10,width = 8)
ggsave("../CAM/Analysis Review 1/figures/PSMs_S2_noMario.pdf",height = 10,width = 8)
