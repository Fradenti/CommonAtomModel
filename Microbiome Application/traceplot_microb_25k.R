m1 <- readRDS("Functions/DCAM_LS/DCAM_results_Dietswap_time1_25k_Feb2021_ab_stochastic.RDS")



par(mfrow=c(1,2))
acf(m1$A_DP,main="ACF alpha")
acf(m1$B_DP,main="ACF beta")

w1 <- qplot(x=1:25000,y=m1$A_DP,geom = "line")+theme_bw()+ggtitle(bquote("Distributional SB - Concentration Parameter"~alpha))+xlab("Iteration")+ylab("Parameter")
w2 <- qplot(x=1:25000,y=m1$B_DP,geom = "line")+theme_bw()+ggtitle(bquote("Observational SB  - Concentration Parameter"~beta))  +xlab("Iteration")+ylab("Parameter")



library(forecast)
library(tidyverse)
library(patchwork)
a1 <- ggAcf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("ACF - Concentration Parameter"~alpha))
a2 <- ggPacf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("PACF - Concentration Parameter"~alpha))
b1 <- ggAcf(m1$B_DP)+theme_bw()+
  ggtitle(bquote("ACF - Concentration Parameter"~beta))
b2 <- ggPacf(m1$B_DP)+theme_bw()+
  ggtitle(bquote("PACF - Concentration Parameter"~beta))




w1/(a1+a2)
ggsave("Functions/DCAM_LS/alpha_25k_microb.pdf",width = 10,height = 9)
ggsave("Functions/DCAM_LS/alpha_25k_microb.png",width = 10,height = 9)


w2/(b1+b2)
ggsave("Functions/DCAM_LS/beta_25k_microb.pdf",width = 10,height = 9)
ggsave("Functions/DCAM_LS/beta_25k_microb.png",width = 10,height = 9)




coda::geweke.diag(m1$A_DP)











w1 <- qplot(x=1:25000,y=m1$NN.cz[,1],geom = "line")+theme_bw()+ggtitle(bquote("Distributional SB - Distributional Cluster Threshold"))+xlab("Iteration")+ylab("Parameter")
w2 <- qplot(x=1:25000,y=m1$NN.cz[,2],geom = "line")+theme_bw()+ggtitle(bquote("Observational SB - Observational Cluster Threshold"))  +xlab("Iteration")+ylab("Parameter")



library(forecast)
library(tidyverse)
library(patchwork)
a1 <- ggAcf(m1$NN.cz[,1])+theme_bw()+
  ggtitle(bquote("ACF - Distributional Cluster Threshold"))
a2 <- ggPacf(m1$NN.cz[,1])+theme_bw()+
  ggtitle(bquote("PACF - Distributional Cluster Threshold"))
b1 <- ggAcf(m1$NN.cz[,2])+theme_bw()+
  ggtitle(bquote("ACF - Observational Cluster Threshold"))
b2 <- ggPacf(m1$NN.cz[,2])+theme_bw()+
  ggtitle(bquote("PACF - Observational Cluster Threshold"))




w1/(a1+a2)
ggsave("Functions/DCAM_LS/SogliaCluster_25k_microb.pdf",width = 10,height = 9)
ggsave("Functions/DCAM_LS/SogliaCluster_25k_microb.png",width = 10,height = 9)


w2/(b1+b2)
ggsave("Functions/DCAM_LS/SogliaDistrib_25k_microb.pdf",width = 10,height = 9)
ggsave("Functions/DCAM_LS/SogliaDistrib_25k_microb.png",width = 10,height = 9)

