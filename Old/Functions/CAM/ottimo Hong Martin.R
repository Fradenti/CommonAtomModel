
##-----------------------------------------------------------------------------------------
## TITLE:      dpmslice
##
## AUTHORS:    L. Hong and R. Martin
##
## OUTLINE:    R code to implement the "independent slice-efficient" algorithm in Kalli,
##             Griffin, and Walker (Stat Comp 2011) for fitting a Dirichlet process mixture
##             of normals; see Section 2 (p 95) of their paper.
##
## REFERENCE:  Hong and Martin (2016) "A flexible Bayesian nonparametric model for 
##             predicting future insurance claims," North Amer Actuarial Journal.
##
##-----------------------------------------------------------------------------------------

M <- 10000
Y <- yd


dpmslice <- function(Y, M, alpha=1, kappa=0.5) {
  
  n <- length(Y)
  B <- round(0.2 * M)
  data.range <- range(Y)
  mu0 <- mean(data.range)
  sigma20 <- sum(c(-1, 1) * data.range)
  gamma <- 2
  beta <- 0.02 * sigma20
  Ynew <- ncomp <- numeric(M)
  y <- seq(data.range[1] - sd(Y), data.range[2] + sd(Y), len=250)
  fy <- 0 * y
  km.out <- kmeans(Y, 7)
  d <- km.out$cluster
  mu <- as.numeric(km.out$centers)
  mu.tmp <- numeric(n)
  g <- function(j) (1 - kappa) * kappa**(j - 1)
  for(r in 1:(M + B)) {
    
    u <- runif(n, 0, g(d))
    L <- 1 + floor((log(u) - log(1 - kappa)) / log(kappa))
    J <- max(d)
    NN <- max(c(L, J))
    xi <- g(1:NN)
    if(NN > J) mu <- c(mu, numeric(NN - J))
    sigma2 <- v <- w <- numeric(NN)
    if(r > B) dpred <- 0 * y
    for(j in 1:NN) {
      
      YY <- Y[d == j]
      nj <- length(YY)
      if(nj > 0) ybarj <- mean(YY) else ybarj <- 0
      if(nj > 1) ssj <- (nj - 1) * var(YY) else ssj <- 0
      sigma2[j] <- 1 / rgamma(1, shape=gamma + nj / 2, rate=beta + (ssj + nj * (ybarj - mu[j])**2) / 2)
      muj <- (mu0 / sigma20 + nj * ybarj / sigma2[j]) / (1 / sigma20 + nj / sigma2[j])
      sigj <- 1 / sqrt(1 / sigma20 + nj / sigma2[j])
      mu[j] <- rnorm(1, muj, sigj)
      v[j] <- rbeta(1, 1 + sum(d == j), alpha + sum(d > j))
      if(j == 1) w[j] <- v[j] else w[j] <- v[j] * prod(1 - v[1:(j - 1)])
      if(r > B) dpred <- dpred + w[j] * dnorm(y, mu[j], sqrt(sigma2[j]))
      
    }
    if(r > B) {
      
      rr <- r - B
      K <- sample(NN, size=1, prob=w)
      Ynew[rr] <- rnorm(1, mu[K], sqrt(sigma2[K]))
      fy <- (rr - 1) * fy / rr + dpred / rr
      ncomp[rr] <- J
      
    }
    for(i in 1:n) {
      
      p <- (xi > u[i]) * w * dnorm(Y[i], mu, sqrt(sigma2)) / xi
      mu.tmp[i] <- mu[sample(NN, size=1, prob=p)]
      
    }
    mu <- unique(mu.tmp)
    d <- apply(outer(mu.tmp, mu, "-"), 1, function(x) which(x == 0));
    
  }
  return(list(Ynew=Ynew, J=ncomp, y=y, fy=fy))
  
}


Y <- yd
# EXAMPLE: Three-point gamma mixture model in Section 5.1 of Jeon & Kim (NAAJ, 2013)

n <- 1000
set.seed(7)
u <- runif(n)
X <- (u <= 0.3) * rgamma(n, 5) + (u > 0.3 & u <= 0.7) * rgamma(n, 25) + (u > 0.7) * rgamma(n, 60)
hist(X, freq=FALSE, breaks=45, col="gray", border="white", main="")
o <- dpmslice(log(X), M, alpha, kappa)
x <- exp(o$y)
fx <- o$fy / x
lines(x, fx)
b <- 0.054
gkde.tmp <- outer(x, X, function(x, y) dgamma(x, shape=y / b + 1, scale=b))
gkde <- apply(gkde.tmp, 1, mean)
lines(x, gkde, lty=2)
J <- sort(unique(o$J))
pJ <- as.numeric(table(o$J)) / M
plot(J, pJ, type="h", xlab="Mixture order", ylab="Posterior probability")
