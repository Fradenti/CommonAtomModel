library(tidyverse)
library(coda)
library(patchwork)
Rcpp::sourceCpp("/home/fra/Documents/GitHub/CommonAtomModel/Functions/CAM/newCAM.cpp")
source("/home/fra/Documents/GitHub/CommonAtomModel/Functions/CAM/newCAM.R")

set.seed(12345)

all <- readRDS("Simulated_Data/ALL_S1A_100.RDS")
Y <- all[[1]][[1]][[3]]
G <- all[[2]][[1]][[3]]

plot(Y,col=G)

set.seed(1234589)
m1 <- CAM(y_obser = Y,
          y_group = G,
          K0 = 10,
          L0 = 10,
          prior = list(
            # hyperparameters NIG
            m0= mean(Y), 
            k0=1/var(Y), a0=3, b0=1,
            # hyperparameters alpha and beta
            a_alpha=3, b_alpha =3,
            a_beta =3, b_beta = 3),
          nsim = 10000,
          burn_in = 10000,
          thinning = 1,
          verbose= 50,
          fixedAB = F,
          conditional.beta = T,
          alpha.fixed = 1,
          kappa=0.5,
          cheap=F,
          seed=1324,
          sigma.prop.beta = 1,
          post.dens = T
)
#saveRDS(m1,"../CAM/Analysis Review 1/figures/M1_PLOTDENS_10K.RDS")

plot(ts(m1$NN.cz))
abline(h=12)
matplot(cbind(
unlist(apply(m1$Z_j,1,   function(x)  length(unique(x)) )),
unlist(apply(m1$Csi_ij,1,function(x)  length(unique(x)) ))),type="l",lty=1 )
abline(h=6)

table(apply(m1$Z_j,1,   function(x)  length(unique(x)) ))

matplot(m1$A_DP,type="l")
matplot(m1$B_DP,type="l")
acf(m1$B_DP)

coda::cumuplot(mcmc(m1$A_DP))
coda::cumuplot(mcmc(m1$B_DP))


image(PSM(m1$Z_j))


mcclust.ext::minVI(PSM(m1$Z_j),method = "greedy")


LLL <- apply(m1$post.dens,1:2,mean)
xx <- seq(-5,25,by = .05)
densit <- matrix(NA,length(xx),1)


plots <- function(m1,k,xx,iter=500,DD=NA){
  
  xx   <- m1$xx
  LLL  <- m1$post.dens
  NN   <- LLL[,k,]   
  nsim <- dim(NN)[2]
  MEAN <- rowMeans(NN)
  NN2  <- NN[,(nsim-iter+1):nsim]
  VV   <- reshape2::melt((NN2))
  VV <- VV %>% mutate(X=rep(xx,iter))
  pp <- ggplot()+theme_bw()+
    geom_line(data=VV,aes(x=X,y=value,group=Var2),col="gray",alpha=.4)
  if(!is.na(DD[1])){
    pp <- pp + geom_line(aes(x=xx,y=DD),col=1,lwd=1)
  }
  pp <-   pp + geom_line(aes(x=xx,y=MEAN),col=2,lwd=1)+xlab("Support")+
    ylab("Density")+ggtitle(paste("Unit",k))
  
  return(pp)
  
}


med <- c(0,5,10,13,16,20)
dd0  <- sapply(med,function(t)  dnorm(m1$xx,t,sqrt(0.6)) )
dd1 <- (dd0[,1])
dd2 <- rowMeans(dd0[,1:2])
dd3 <- rowMeans(dd0[,1:3])
dd4 <- rowMeans(dd0[,1:4])
dd5 <- rowMeans(dd0[,1:5])
dd6 <- rowMeans(dd0[,1:6])



g1  <- plots(m1,k=1,DD = dd1,xx,iter = 250)
g2  <- plots(m1,k=2,DD = dd2,xx,iter = 250)
g3  <- plots(m1,k=3,DD = dd3,xx,iter = 250)
g4  <- plots(m1,k=4,DD = dd4,xx,iter = 250)
g5  <- plots(m1,k=5,DD = dd5,xx,iter = 250)
g6  <- plots(m1,k=6,DD = dd6,xx,iter = 250)
g7  <- plots(m1,k=7 ,DD = dd1,xx,iter = 250 )
g8  <- plots(m1,k=8 ,DD = dd2,xx,iter = 250 )
g9  <- plots(m1,k=9 ,DD = dd3,xx,iter = 250 )
g10 <- plots(m1,k=10,DD = dd4,xx,iter = 250)
g11 <- plots(m1,k=11,DD = dd5,xx,iter = 250)
g12 <- plots(m1,k=12,DD = dd6,xx,iter = 250)



qplot(x = 1:length(m1$y_group),y=m1$y_obser,col=factor(m1$y_group))+theme_bw()+
  scale_color_viridis_d("Units",option = "B",end = .9)+geom_vline(aes(xintercept=c(0,cumsum(table((m1$y_group))))),lty=2)+
  xlab("Observation index")+ylab("Observation value")
ggsave("Other studies/Feb2021_12units_scenario1.png",width = 12,height = 7.5)
ggsave("Other studies/Feb2021_12units_scenario1.pdf",width = 12,height = 7.5)




library(patchwork)
(g1+g2+g3)/(g4+g5+g6)#/(g7+g8+g9)/(g10+g11+g12)
ggsave("Other studies/Recovered_Estimated_10k.png",width = 15,height = 10)
ggsave("Other studies/Recovered_Estimated_10k.pdf",width = 15,height = 10)



w1 <- qplot(x=1:10000,y=m1$A_DP,geom = "line")+theme_bw()+ggtitle(bquote("Distributional DP - Concentration Parameter"~alpha))+xlab("Iteration")+ylab("Parameter")
w2 <- qplot(x=1:10000,y=m1$B_DP,geom = "line")+theme_bw()+ggtitle(bquote("Observational DP  - Concentration Parameter"~beta))  +xlab("Iteration")+ylab("Parameter")


w1/w2
#ggsave("../CAM/Analysis Review 1/figures/alphabeta_10k.pdf",width = 10,height = 7)
#ggsave("../CAM/Analysis Review 1/figures/alphabeta_10k.png",width = 10,height = 7)


par(mfrow=c(1,2))
acf(m1$A_DP,main="ACF alpha")
acf(m1$B_DP,main="ACF beta")


matplot(m1$pi_star_k)



library(forecast)
a1 <- ggAcf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("ACF - Concentration Parameter"~alpha))
a2 <- ggPacf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("PACF - Concentration Parameter"~alpha))
b1 <- ggAcf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("ACF - Concentration Parameter"~beta))
b2 <- ggPacf(m1$A_DP)+theme_bw()+
  ggtitle(bquote("PACF - Concentration Parameter"~beta))




w1/(a1+a2)
#ggsave("../CAM/Analysis Review 1/figures/alpha_10k.pdf",width = 10,height = 9)
#ggsave("../CAM/Analysis Review 1/figures/alpha_10k.png",width = 10,height = 9)


w2/(b1+b2)
#ggsave("../CAM/Analysis Review 1/figures/beta_10k.pdf",width = 10,height = 9)
#ggsave("../CAM/Analysis Review 1/figures/beta_10k.png",width = 10,height = 9)




geweke.diag(m1$A_DP)

plot(ts(unlist(map(m1$Z_j,~length(unique(.x))))))
plot(ts(unlist(map(m2$Z_j,~length(unique(.x))))))
          
     
