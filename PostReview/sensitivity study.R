library(tidyverse)
library(patchwork)


set.seed(12345)
Y1 <- c(rnorm(50,0,.5),rnorm(50,2,.5),rnorm(50,8,1))
Y2 <- c(rnorm(50,2,1),rnorm(50,3,1),rnorm(50,4,1))
Y3 <- c(rnorm(50,0,.5),rnorm(50,2,.5),rnorm(50,8,1))
Y4 <- c(rnorm(50,2,1),rnorm(50,3,1),rnorm(50,-5,1))



YG <- rep(1:4,rep(150,4))
YD <- c(Y1,Y2,Y3,Y4)#,Y5)
plot(YD)
plot(YD,col=YG)


p1 <- list(
  # hyperparameters NIG
  m0=0, k0=100, a0=10, b0=10,
  # hyperparameters alpha and beta
  a_alpha=10, b_alpha = 10,
  a_beta =10, b_beta = 10 )

p2 <- list(
  # hyperparameters NIG
  m0=0, k0=.001, a0=1, b0=1,
  # hyperparameters alpha and beta
  a_alpha=.001, b_alpha = .001,
  a_beta =.001, b_beta = .001 )

p3 <- list(
  # hyperparameters NIG
  m0=mean(YD), k0=1/var(YD), a0=3, b0=1,
  # hyperparameters alpha and beta
  a_alpha=3, b_alpha = 3,
  a_beta =3, b_beta = 3 )




m1 <- CAM(y_obser = YD,
          y_group = YG,
          K0 = 20,
          L0 = 20,
          prior = p1,
          NSIM = 3000,
          burn_in = 3000,
          thinning = 1,
          verbose.step = 50,
          fixedAB = F,
          kappa=0.5,cheap=F,
          seed=1234
)

m2 <- CAM(y_obser = YD,
          y_group = YG,
          K0 = 20,
          L0 = 20,
          prior = p2,
          NSIM = 3000,
          burn_in = 3000,
          thinning = 1,
          verbose.step = 50,
          fixedAB = F,
          kappa=0.5,cheap=F,
          seed=1234
)

m3 <- CAM(y_obser = YD,
          y_group = YG,
          K0 = 20,
          L0 = 20,
          prior = p3,
          NSIM = 3000,
          burn_in = 3000,
          thinning = 1,
          verbose.step = 50,
          fixedAB = F,
          kappa=0.5,cheap=F,
          seed=1234
)



J <- 4
xx <- seq(-8,12,by = .05)


LLL1 <- list()
for(h in 1:J){
  NN <- matrix(NA,1000,length(xx))
  for(i in 1:1000){
    NN[i,] <-   colSums((m1$omega_star_lk[[i]][,m1$Z_j[[i]]][,h])*
                          t(apply(m1$theta_star_with_prior[[i]],1,function(z) dnorm(xx,z[1],sqrt(z[2])) )))}
  LLL1[[h]] <- NN
  cat(h)
}

LLL2 <- list()
for(h in 1:J){
  NN <- matrix(NA,1000,length(xx))
  for(i in 1:1000){
    NN[i,] <-   colSums((m2$omega_star_lk[[i]][,m2$Z_j[[i]]][,h])*
                          t(apply(m2$theta_star_with_prior[[i]],1,function(z) dnorm(xx,z[1],sqrt(z[2])) )))}
  LLL2[[h]] <- NN
  cat(h)
}

LLL3 <- list()
for(h in 1:J){
  NN <- matrix(NA,1000,length(xx))
  for(i in 1:1000){
    NN[i,] <-   colSums((m3$omega_star_lk[[i]][,m3$Z_j[[i]]][,h])*
                          t(apply(m3$theta_star_with_prior[[i]],1,function(z) dnorm(xx,z[1],sqrt(z[2])) )))}
  LLL3[[h]] <- NN
  cat(h)
}

# matplot(xx,t(NN),type="l",col="gray")
# matplot(xx,colMeans((NN)),type="l",col=2,add=T,lwd=2)
# matplot(xx,dd1,type="l",col=1,add=T,lwd=2)
# matplot(density(m1$y_obser[m1$y_group==6],bw = .5)$x,
#         density(m1$y_obser[m1$y_group==6],bw = .5)$y,type="l",
#         col=4,add=T,lwd=2)
# 


dY1 <- 1/3*dnorm(xx,0,.5)+1/3*dnorm(xx,2,.5)+1/3*dnorm(xx,8,1)
dY2 <- 1/3*dnorm(xx,2,1) +1/3*dnorm(xx,3,1) +1/3*dnorm(xx,4,1)
dY3 <- 1/3*dnorm(xx,0,.5)+1/3*dnorm(xx,2,.5)+1/3*dnorm(xx,8,1)
dY4 <- 1/3*dnorm(xx,2,1) +1/3*dnorm(xx,3,1) +1/3*dnorm(xx,-5,1)



# med <- c(0,5,10,13,16,20)
# dd0  <- sapply(med,function(t)  dnorm(xx,t,sqrt(0.6)) )
# dd1 <- (dd0[,1])
# dd2 <- rowMeans(dd0[,1:2])
# dd3 <- rowMeans(dd0[,1:3])
# dd4 <- rowMeans(dd0[,1:4])


plots <- function(k,LLL,DD){
  
  NN <- LLL[[k]]   
  VV <- reshape2::melt(t(NN))
  VV <- VV %>% mutate(X=rep(xx,1000))
  pp <- ggplot()+theme_bw()+
    geom_line(data=VV,aes(x=X,y=value,group=Var2),col="gray",alpha=.4)+
    geom_line(aes(x=xx,y=DD),col=1,lwd=1)+
    geom_line(aes(x=xx,y=colMeans(NN)),col=2,lwd=1)+xlab("Support")+
    ylab("Density")+ggtitle(paste("Unit",k))
  
  return(pp)
  
}

g1.1  <- plots(1,LLL = LLL1, DD=dY1)+ggtitle("Unit 1 - Prior Set 1")
g1.2  <- plots(1,LLL = LLL2, DD=dY1)+ggtitle("Unit 1 - Prior Set 2")
g1.3  <- plots(1,LLL = LLL3, DD=dY1)+ggtitle("Unit 1 - Prior Set 3")
(g1.1/g1.2)/g1.3


g2.1  <- plots(2,LLL = LLL1, DD=dY2)+ggtitle("Unit 2 - Prior Set 1")
g2.2  <- plots(2,LLL = LLL2, DD=dY2)+ggtitle("Unit 2 - Prior Set 2")
g2.3  <- plots(2,LLL = LLL3, DD=dY2)+ggtitle("Unit 2 - Prior Set 3")
g2.1/g2.2/g2.3


g3.1  <- plots(3,LLL = LLL1, DD=dY3)+ggtitle("Unit 3 - Prior Set 1")
g3.2  <- plots(3,LLL = LLL2, DD=dY3)+ggtitle("Unit 3 - Prior Set 2")
g3.3  <- plots(3,LLL = LLL3, DD=dY3)+ggtitle("Unit 3 - Prior Set 3")
g3.1/g3.2/g3.3


g4.1  <- plots(4,LLL = LLL1, DD=dY4)+ggtitle("Unit 4 - Prior Set 1") 
g4.2  <- plots(4,LLL = LLL2, DD=dY4)+ggtitle("Unit 4 - Prior Set 2")
g4.3  <- plots(4,LLL = LLL3, DD=dY4)+ggtitle("Unit 4 - Prior Set 3")
g4.1/g4.2/g4.3




(g1.1+g1.2+g1.3)/
(g2.1+g2.2+g2.3)/
(g3.1+g3.2+g3.3)/
(g4.1+g4.2+g4.3)
ggsave("/home/fra/Documents/GitHub/CommonAtomModel/PostReview/Sensivity.png",width = 10,height = 12)
ggsave("/home/fra/Documents/GitHub/CommonAtomModel/PostReview/Sensivity.pdf",width = 10,height = 12)

MM <- list(m1,m2,m3)
saveRDS(MM,"/home/fra/Documents/GitHub/CommonAtomModel/PostReview/Sensivity_datasets.RDS")
MM <- readRDS("/home/fra/Documents/GitHub/CommonAtomModel/PostReview/Sensivity_datasets.RDS")
m1 <- MM[[1]]
m2 <- MM[[2]]
m3 <- MM[[3]]



mcclust::comp.psm(matrix(unlist(m1$Z_j),3000,4,byrow = T))
mcclust::comp.psm(matrix(unlist(m2$Z_j),3000,4,byrow = T))
mcclust::comp.psm(matrix(unlist(m3$Z_j),3000,4,byrow = T))

for_psm1 <- matrix(unlist(m1$Csi_ij),3000,600,byrow = T)
for_psm2 <- matrix(unlist(m2$Csi_ij),3000,600,byrow = T)
for_psm3 <- matrix(unlist(m3$Csi_ij),3000,600,byrow = T)
psm1 <- mcclust::comp.psm(for_psm1)
psm2 <- mcclust::comp.psm(for_psm2)
psm3 <- mcclust::comp.psm(for_psm3)
image(psm1)
image(psm2)
image(psm3)
cl1 <- mcclust.ext::minVI(psm1)$cl
cl2 <- mcclust.ext::minVI(psm2)$cl
cl3 <- mcclust.ext::minVI(psm3)$cl


cl1 <- GreedyEPL::MinimiseEPL(for_psm1[2000:3000,])
plot(YD,col=cl1$decision)
cl2 <- GreedyEPL::MinimiseEPL(for_psm2[2000:3000,])
plot(YD,col=cl2$decision)
cl3 <- GreedyEPL::MinimiseEPL(for_psm3[2000:3000,])
plot(YD,col=cl3$decision)




# -------------------------------------------------------------------------


psm0 <- mcclust::comp.psm(rbind(rep(c(1,2,3,4,5,6,1,2,3,4,5,7),rep(50,12)),rep(c(1,2,3,4,5,6,1,2,3,4,5,7),rep(50,12))))
vv0 <- reshape2::melt(psm0)
PSM0 <- ggplot(data=as_tibble(vv0))+geom_raster(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Observation")+ylab("Observation")+scale_fill_viridis_c(option = "B")+
  ggtitle("Ground Truth")

vv1 <- reshape2::melt(psm1)
PSM1 <- ggplot(data=as_tibble(vv1))+geom_raster(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Observation")+ylab("Observation")+scale_fill_viridis_c(option = "B")+
  ggtitle("Prior Set 1")

vv2 <- reshape2::melt(psm2)
PSM2 <- ggplot(data=as_tibble(vv2))+geom_raster(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Observation")+ylab("Observation")+scale_fill_viridis_c(option = "B")+
  ggtitle("Prior Set 2")


vv3 <- reshape2::melt(psm3)
PSM3 <- ggplot(data=as_tibble(vv3))+geom_raster(aes(x=Var1, y=Var2, fill=value))+
  theme_bw()+xlab("Observation")+ylab("Observation")+scale_fill_viridis_c(option = "B")+
  ggtitle("Prior Set 3")
PSM3
(PSM0+PSM1)/(PSM2+PSM3)
ggsave("PostReview/Coclustering.png",width = 10,height = 10)
ggsave("PostReview/Coclustering.pdf",width = 10,height = 10)



# -------------------------------------------------------------------------


gt <- rep(c(1,2,3,4,5,6,1,2,3,4,5,7),rep(50,12))
DDf <- reshape2::melt(cbind(Y=YD,cl1,cl2,cl3,gt),"Y")
DDf <- as_tibble(DDf)


gg0 <- ggplot()+  ggtitle("Ground Truth")+
  geom_point(aes(x=1:600,y=YD,col=factor(gt)))+scale_color_viridis_d("OC",option = "B",end = .9)+theme_bw()+xlab("Observation")+ylab("")
gg1 <- ggplot()+  ggtitle("Prior Set 1")+
  geom_point(aes(x=1:600,y=YD,col=factor(cl1)))+scale_color_viridis_d("Estimated\nOC",option = "B",end = .9)+theme_bw()+xlab("Observation")+ylab("")
gg2 <- ggplot()+  ggtitle("Prior Set 2")+
  geom_point(aes(x=1:600,y=YD,col=factor(cl2)))+scale_color_viridis_d("Estimated\nOC",option = "B",end = .9)+theme_bw()+xlab("Observation")+ylab("")
gg3 <- ggplot()+  ggtitle("Prior Set 3")+
  geom_point(aes(x=1:600,y=YD,col=factor(cl3)))+scale_color_viridis_d("Estimated\nOC",option = "B",end = .9)+theme_bw()+xlab("Observation")+ylab("")


(gg0+gg1)/(gg2+gg3)
ggsave("PostReview/points.png",width = 10,height = 10)
ggsave("PostReview/points.pdf",width = 10,height = 10)
