# R version 4.0.3

# Dataset created under three scenarios. Each object is a list containing 100 lists
# Each of the 100 is a list of three datasets with an increasing number of observations

# Useful to understand the structure:
# ALL[[a]][[b]][[c]]
# a = 1,2,3 either data, groups, or observational clusters
# b = 1,...,100 one of the 100 configurations
# c = 1,2,3 (or 6 in S2 and S3), one of the different sample sizes
# Useful functions

CaseA <- function(m,s=sqrt(.6),N){
  p <- c()
  xx <- c()
  for(i in 1:6){
    x <- sample(1:i,N,T)
    xx <- c(xx,x)
    p <- c(p,rnorm(N,m[x],s)) }
  return(cbind(p,xx))
}

CaseB1 <- function(m,s=sqrt(.6), nsim=100){
  ALL <- sapply(1:6,function(xx) sapply(1:xx,
                                        function(x)  rnorm(nsim,m[x],s)))
  return(unlist(ALL))
}

CaseB <- function(m,s=sqrt(.6),N){
  y1 <- c();
  for(ll in 1:2){
    y1 <- c(y1,CaseB1(med,nsim = N))
  }
  return(y1)
}

# Scenario 1 A

med   <- c(0,5,10,13,16,20)
NUM_A <- c(25,50,75)

Y_100 <- list()
G_100 <- list()
O_100 <- list()

set.seed(22012021)
for(k in 1:100){
  
  OCs <- YD <- YG <- list()
  
  for(i in 1:3){
    y1 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
    y2 <- CaseA(m = med,s = sqrt(.6),N = NUM_A[i])
    YD[[i]] <- c(y1[,1],y2[,1])
    OCs[[i]]     <- c(y1[,2],y2[,2])
    repet   <- rep(NUM_A[i],12) #rep(200,6)
    y_group <- c(rep(1:length(repet),repet))
    YG[[i]] <- y_group
  }
  
  Y_100[[k]] <- YD
  G_100[[k]] <- YG 
  O_100[[k]] <- OCs
}


plot(Y_100[[1]][[3]],col=G_100[[1]][[3]])
ALL_S1a <- list(Y_100,
                G_100,
                O_100)
saveRDS(ALL_S1a,"../CommonAtomModel/Simulated_Data/ALL_S1A_100.RDS")


###############################################################
# Scenario 1 B
###############################################################
med   <- c(0,5,10,13,16,20)
NUM_B <- c(5,10,20)

Y_100 <- list()
G_100 <- list()
O_100 <- list()


set.seed(123456789)
for(k in 1:100){
  OCs <- YD <- YG <- list()
  
  for(i in 1:3){
    y1 <- CaseB(m = med,s = sqrt(.6),N = NUM_B[i])
    OCs[[i]] <- unlist(sapply(1:6, function(x) rep(1:x,rep(NUM_B[i],x)  )))
    YD[[i]] <- c(y1)
    repet   <- c((1:6)*NUM_B[i],(1:6)*NUM_B[i]) #rep(200,6)
    y_group <- c(rep(1:length(repet),repet))
    YG[[i]] <- y_group
  }
  
  Y_100[[k]] <- YD
  G_100[[k]] <- YG 
  O_100[[k]] <- OCs
}

plot(Y_100[[1]][[3]],col=G_100[[1]][[3]])
plot(Y_100[[1]][[3]],col=O_100[[1]][[3]])

ALL_S1b <- list(Y_100,
                G_100,
                O_100)
saveRDS(ALL_S1b,"../CommonAtomModel/Simulated_Data/ALL_S1B_100.RDS")





###############################################################
# Scenario 2
###############################################################
NNN     <- 40
lab <- c(rep(1,NNN*.75),  rep(2,NNN*.25),
         rep(1,NNN*.25),  rep(2,NNN*.75),
         rep(1,NNN*.33),  rep(3,NNN*.34),rep(4,NNN*.33),
         rep(1,NNN*.25),rep(4,NNN*.25),
         rep(3,NNN*.25+1),  rep(5,NNN*.25))
GT_O <- list()
GT_O[[1]] <- lab
GT_O[[2]] <- c(lab,GT_O[[1]])
GT_O[[3]] <- c(lab,GT_O[[2]])
GT_O[[4]] <- c(lab,GT_O[[3]])
GT_O[[5]] <- c(lab,GT_O[[4]])
GT_O[[6]] <- c(lab,GT_O[[5]])

Y_100 <- list()
G_100 <- list()
O_100 <- list()

set.seed(22012021)
for(k in 1:100){
  
  OCs <- YD <- YG <- list()
  for(i in 1:6){
    YD[[i]] <- c(replicate(i,{c( rnorm(NNN*.75,0,sqrt(.6)),rnorm(NNN*.25, 3,sqrt(.6)),
                                 rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.75, 3,sqrt(.6)),
                                 rnorm(NNN*.33,0,sqrt(.6)),rnorm(NNN*.34, 2,sqrt(.6)),
                                 rnorm(NNN*.33,-2,sqrt(.6)),
                                 rnorm(NNN*.25,0,sqrt(.6)), rnorm(NNN*.25, -2,sqrt(.6)),
                                 rnorm(NNN*.25+1,2,sqrt(.6)),  rnorm(NNN*.25, 10, sqrt(.6)))}))
    rep     <- rep(NNN,4*i)
    YG[[i]] <- rep(1:length(rep),rep)
    OCs[[i]]<- GT_O[[i]]
    
  }
  Y_100[[k]] <- YD
  G_100[[k]] <- YG 
  O_100[[k]] <- OCs
}


plot(Y_100[[1]][[3]],col=G_100[[1]][[3]])
plot(Y_100[[1]][[3]],col=O_100[[1]][[3]])

ALL_S2 <- list(Y_100,
                G_100,
                O_100)
saveRDS(ALL_S2,"../CommonAtomModel/Simulated_Data/ALL_S2_100.RDS")







##################################################################################
# Data generation ---------------------------------------------------------
num <- c(5,10,15,25,50,75,100)

Y_100 <- list()
G_100 <- list()
O_100 <- list()

set.seed(22012021)
for(k in 1:100){
  
  OCs <- YD <- YG <- list()
  for(i in 1:6){
    YD[[i]] <- c( rep(0,50),rep(1,50),sample(50, num[i],T),
                  rep(0,50),rep(1,50),sample(100,num[i],T),
                  rep(0,50),rep(1,50),sample(10, num[i],T),
                  rep(0,50),rep(1,50),sample(100,num[i],T),
                  rep(0,50),rep(1,50),sample(50, num[i],T),
                  rep(0,50),rep(1,50),sample(100,num[i],T),
                  rep(0,50),rep(1,50),sample(50, num[i],T),
                  rep(0,50),rep(1,50),sample(10, num[i],T),
                  rep(0,50),rep(1,50),sample(100,num[i],T),
                  rep(0,50),rep(1,50),sample(50, num[i],T))
    L <- YD[[i]]
    OCs[[i]] <- ifelse(L==0 | L==1,1, 
                      ifelse(L<=10,2,
                             ifelse(L<=50,3,4)  ))
    repet   <- rep(100+num[i],10)
    YG[[i]] <- c(rep(1:length(repet),repet))
  }
  Y_100[[k]] <- YD
  G_100[[k]] <- YG 
  O_100[[k]] <- OCs 
}


plot(Y_100[[1]][[3]],col=G_100[[1]][[3]])
plot(Y_100[[1]][[3]],col=O_100[[1]][[3]])

ALL_S3 <- list(Y_100,
               G_100,
               O_100)
saveRDS(ALL_S3,"../CommonAtomModel/Simulated_Data/ALL_S3_100.RDS")
