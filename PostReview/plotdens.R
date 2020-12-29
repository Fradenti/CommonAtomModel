YD <- YG <- list()
set.seed(123456)
med <- c(0,5,10,13,16,20)
CaseA <- function(m,s=sqrt(.6),N){
  p <- c()
  xx <- c()
  for(i in 1:6){
    x <- sample(1:i,N,T)
    xx <- c(xx,x)
    p <- c(p,rnorm(N,m[x],s)) }
  return(cbind(p,xx))
}

NUM_A <- c(25,50,75)

for(i in 1:3){
  y1 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  y2 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
  YD[[i]] <- c(y1[,1],y2[,1])
  repet   <- rep(NUM_A[i],12) #rep(200,6)
  y_group <- c(rep(1:length(repet),repet))
  YG[[i]] <- y_group
}


plot(YD[[1]],col=YG[[1]])
plot(YD[[2]],col=YG[[2]])
plot(YD[[3]],col=YG[[3]])



m1 <- CAM(y_obser = YD[[3]],
          y_group = YG[[3]],
    K0 = 10,
    L0 = 10,
    prior = list(
      # hyperparameters NIG
      m0=0, k0=1/var(YD[[3]]), a0=3, b0=1,
      # hyperparameters alpha and beta
      a_alpha=3, b_alpha = 3,
      a_beta =3, b_beta = 3 ),
    NSIM = 1000,
    burn_in = 5000,
    thinning = 1,
    verbose.step = 50,
    fixedAB = F,
    kappa=0.5,cheap=F,
    seed=123456
)

# m1 <- CAM(y_obser = YD[[3]],
#           y_group = YG[[3]],
#           K0 = 10,
#           L0 = 10,
#           prior = list(
#             # hyperparameters NIG
#             m0=0, k0=1, a0=1, b0=1,
#             # hyperparameters alpha and beta
#             a_alpha=1, b_alpha = 1,
#             a_beta =1, b_beta = 1 ),
#           NSIM = 5000,
#           burn_in = 10000,
#           thinning = 1,
#           verbose.step = 50,
#           fixedAB = F,
#           kappa=0.5,cheap=F,seed=123456
# )


ZZ <- matrix(unlist(m1$Z_j),1000,12,byrow = T)
image(mcclust::comp.psm(ZZ))



# m1$theta_star_zero[[1000]]
# m1$omega_star_lk[[1000]]
# m1$pi_star_k[[1000]]
# 
# m1$y_obser[600]
# m1$Csi_ij[[1000]][600]
# m1$pi_star_k[[2]]

dims <- matrix(NA,1000,5)
re <- numeric(1000)
for(i in 1:1000){
dims[i,1] <- dim(m1$theta_star_zero[[i]])[1]
dims[i,2] <- dim(m1$pi_star_k[[i]])[1]
dims[i,3] <- dim(m1$omega_star_lk[[i]])[1]
dims[i,4] <- dim(m1$omega_star_lk[[i]])[2]
dims[i,5] <- dim(m1$theta_star_with_prior[[i]])[1]

re[i] <- sum(1-colSums(m1$omega_star_lk[[i]][1:dims[i,1],]))
}


plot(dims[,1],type="l")
plot(dims[,2],type="l")
lines(dims[,3],type="l")
lines(dims[,4],type="l")
lines(dims[,5],type="l")


plot(dims[,1]<=dims[,3],type="l")
plot(re)
m1$theta_star_with_prior



xx <- seq(-5,25,by = .05)
densit <- matrix(NA,length(xx),1)


med <- c(0,5,10,13,16,20)
dd0  <- sapply(med,function(t)  dnorm(xx,t,sqrt(0.6)) )
dd1 <- (dd0[,1])
dd2 <- rowMeans(dd0[,1:2])
dd3 <- rowMeans(dd0[,1:3])
dd4 <- rowMeans(dd0[,1:4])
dd5 <- rowMeans(dd0[,1:5])
dd6 <- rowMeans(dd0[,1:6])


LLL <- list()
for(h in 1:12){
NN <- matrix(NA,1000,length(xx))
for(i in 1:1000){
  NN[i,] <-   colSums((m1$omega_star_lk[[i]][,m1$Z_j[[i]]][,h])*
  t(apply(m1$theta_star_with_prior[[i]],1,function(z) dnorm(xx,z[1],sqrt(z[2])) )))}
LLL[[h]] <- NN
cat(h)
}


# matplot(xx,t(NN),type="l",col="gray")
# matplot(xx,colMeans((NN)),type="l",col=2,add=T,lwd=2)
# matplot(xx,dd1,type="l",col=1,add=T,lwd=2)
# matplot(density(m1$y_obser[m1$y_group==6],bw = .5)$x,
#         density(m1$y_obser[m1$y_group==6],bw = .5)$y,type="l",
#         col=4,add=T,lwd=2)
# 


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

g1  <- plots(1,LLL = LLL,DD = dd1)
g2  <- plots(2,LLL = LLL,DD = dd2)
g3  <- plots(3,LLL = LLL,DD = dd3)
g4  <- plots(4,LLL = LLL,DD = dd4)
g5  <- plots(5,LLL = LLL,DD = dd5)
g6  <- plots(6,LLL = LLL,DD = dd6)
g7  <- plots(7 ,LLL = LLL,DD = dd1 )
g8  <- plots(8 ,LLL = LLL,DD = dd2 )
g9  <- plots(9 ,LLL = LLL,DD = dd3 )
g10 <- plots(10,LLL = LLL,DD = dd4)
g11 <- plots(11,LLL = LLL,DD = dd5)
g12 <- plots(12,LLL = LLL,DD = dd6)



qplot(x = 1:length(m1$y_group),y=m1$y_obser,col=factor(m1$y_group))+theme_bw()+
  scale_color_viridis_d("Units",option = "B",end = .9)+geom_vline(aes(xintercept=c(0,cumsum(table((m1$y_group))))),lty=2)+
  xlab("Observation index")+ylab("Observation value")
ggsave("../CAM/Analysis Review 1/figures/12units_scenario1.png",width = 12,height = 7.5)
ggsave("../CAM/Analysis Review 1/figures/12units_scenario1.pdf",width = 12,height = 7.5)




library(patchwork)
(g1+g2+g3)/(g4+g5+g6)#/(g7+g8+g9)/(g10+g11+g12)
ggsave("../CAM/Analysis Review 1/figures/Estiamted_.png",width = 15,height = 10)
ggsave("../CAM/Analysis Review 1/figures/Estiamted_.pdf",width = 15,height = 10)



w1 <- qplot(x=1:1000,y=m1$A_DP,geom = "line")+theme_bw()+ggtitle(bquote("Distributional DP - Concentration Parameter"~alpha))+xlab("Iteration")+ylab("Parameter")
w2 <- qplot(x=1:1000,y=m1$B_DP,geom = "line")+theme_bw()+ggtitle(bquote("Observational DP  -Concentration Parameter"~beta))  +xlab("Iteration")+ylab("Parameter")


w1/w2
ggsave("../CAM/Analysis Review 1/figures/alphabeta.pdf",width = 10,height = 7)
ggsave("../CAM/Analysis Review 1/figures/alphabeta.png",width = 10,height = 7)


par(mfrow=c(1,2))
acf(m1$A_DP,main="ACF alpha")
acf(m1$B_DP,main="ACF beta")


matplot(m1$pi_star_k)
