x0 <- replicate(1, c(rnorm(50,3),rnorm( 150,10,.5)))
x1 <- replicate(1, c(rnorm(100,3),rnorm( 100,10,.5)))
x2 <- replicate(3, c(rnorm(100,-5),rnorm(100,10)))
x3 <- c(rnorm(100,3),rnorm( 100,10,.5),rnorm( 100,-5,.5),rnorm( 100,0,.5))

yd <- c(x0,x1,x2,x3)
yg <- rep(1:6,c(rep(200,5),400))
  
plot(yd,col=yg)

  
  m1 <- CAM(y_obser = yd,
            y_group = yg,
            K0 = 10,
            L0 = 10,
            prior = list(
              # hyperparameters NIG
              m0= mean(yd), 
              k0=1/var(yd), 
              a0=3, b0=1,
              # hyperparameters alpha and beta
              a_alpha=3, b_alpha =3,
              a_beta =3, b_beta = 3),
            nsim = 10000,
            burn_in = 1,
            thinning = 1,
            verbose = T,
            fixedAB = F,
            kappa=0.5,cheap=F,
            seed=123456,sigma.prop.beta = 1
  )

  
  
  
  
plot(ts(m1$A_DP)) 
acf(m1$A_DP)
pacf(m1$A_DP)
plot(ts(m1$Sigma))
plot(ts(m1$B_DP))  
acf(m1$B_DP)
pacf(m1$B_DP)


as.matrix(t(map_dfc(m1$Z_j, ~.x)  ))


r=mcclust::comp.psm(as.matrix(t(map_dfc(m1$Csi_ij, ~.x)  )))
image(r)
a=mcclust.ext::minVI(r)


plot(yd,col=a$cl)





xx <- seq(-10,15,by = .05)
densit <- matrix(NA,length(xx),1)


nsim <- 5000
LLL <- list()


LLL <- posterior.densities(m1,howfarback = 500,T,xx = xx)

hh <- 6
hist(yd[yg==hh],breaks = 30,freq=F)

for(i in 1:500)
lines(LLL[[hh]][i,]~xx,col=1)

lines(colMeans(LLL[[hh]])~xx,col=2,lwd=3)

hh <- 6
hist(yd[yg==hh],breaks = 30,freq=F)
lines(colMeans(LLL[[hh]])~xx,col=2,lwd=3)



sim <- 500
ind <- 500
ygind <- yg[[ind]]

plot(yd)
points(yd[ind]~c(ind),col=2,pch=21,bg=2)

mat <- matrix(NA,nsim,2)
for(sim in 2:nsim){

mat[sim,1] <- m1$OMEGA[[sim]][ m1$Csi_ij[[sim]][ind] 
                 ,
                  m1$Z_j[[sim]]  [ygind]]
mat[sim,2] <- m1$THETA[[sim]][ m1$Csi_ij[[sim]][ind] 
                 ,
                 1]
}

plot(ts(mat))




