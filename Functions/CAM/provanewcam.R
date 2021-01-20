x0 <- replicate(1, c(rnorm(50,3),rnorm( 150,10)))
x1 <- replicate(2, c(rnorm(100,3),rnorm( 100,10)))
x2 <- replicate(3, c(rnorm(100,-5),rnorm(100,10)))


yd <- c(x0,x1,x2)
yg <- rep(1:6,rep(200,6))
  
plot(yd,col=yg)

  
  m1 <- CAM(y_obser = yd,
            y_group = yg,
            K0 = 10,
            L0 = 10,
            prior = list(
              # hyperparameters NIG
              m0= mean(YD[[2]]), 
              k0=1/var(YD[[2]]), 
              a0=3, b0=1,
              # hyperparameters alpha and beta
              a_alpha=3, b_alpha =3,
              a_beta =3, b_beta = 3),
            nsim = 5000,
            burn_in = 5000,
            thinning = 1,
            verbose = 50,
            fixedAB = F,
            kappa=0.5,cheap=F,
            seed=123456,sigma.prop.beta = .25
  )

  
plot(ts(m1$A_DP))  
plot(ts(m1$B_DP))  

as.matrix(t(map_dfc(m1$Z_j, ~.x)  ))


r=mcclust::comp.psm(as.matrix(t(map_dfc(m1$Csi_ij, ~.x)  )))
image(r)
a=mcclust.ext::minVI(r)


plot(yd,col=a$cl)






xx <- seq(-10,25,by = .05)
densit <- matrix(NA,length(xx),1)


nsim <- 5000
LLL <- list()
for(h in 1:6){
  NN <- matrix(NA,1000,length(xx))
  for(i in 1:1000){
    
    XXXX <- 
      m1$OMEGA[[nsim-i]][,m1$Z_j[[nsim-i]]][,h] *
      t(apply(m1$THETA[[nsim-i]],1,
              function(z) dnorm(xx,z[1],sqrt(z[2])) ))
    NN[i,] <-   colSums(XXXX)
  }
  LLL[[h]] <- NN
  cat(h)
}


hist(yd[yg==1],breaks = 30,freq=F)
lines(colMeans(LLL[[1]])~xx)

hist(yd[yg==2],breaks = 30,freq=F)
lines(colMeans(LLL[[2]])~xx)


hist(yd[yg==3],breaks = 30,freq=F)
lines(colMeans(LLL[[3]])~xx)


hist(yd[yg==4],breaks = 30,freq=F)
lines(colMeans(LLL[[4]])~xx)


hist(yd[yg==5],breaks = 30,freq=F)
lines(colMeans(LLL[[5]])~xx)




