# tsl0.tmp <- Update_tsl0_for_cij(y_obser = y_obser,
#                              Uij = Uij,
#                              xi_c = xi.c,
#                              omega = omega,
#                              zj_pg = zj.pg,
#                              tsl0 = tsl0,
#                              N = N,
#                              NN_c = NN.c,possible_label = 1:NN.c)
# 
# tsl0     <- subset(tsl0.tmp,!duplicated(tsl0.tmp[,1])) # find active values
# cij      <- apply(outer(tsl0.tmp[,1], tsl0[,1], "-"), 1, function(x) which(x == 0));


# k2      <- length(unique(cij)) # double check this sum(tapply(cij, zj.pg, function(x) length(unique(x)))) # 
# eta2    <- rbeta(1,beta+1, N)  
# Q2      <- (a_beta+k2-1)/(N*(b_beta-log(eta2)))
# pi_eta2 <- Q2/(1+Q2)  
# u2      <- runif(1)
# if(u2<pi_eta2){
#   beta    <-  rgamma(1,a_beta+k2,  b_beta-log(eta2))
#           } else{ 
#   beta    <-  rgamma(1,a_beta+k2-1,b_beta-log(eta2)) 
#   }


#    beta=3


# Kz    <-  length(unique(zj.pg))
# n_j   <-  table(zj.pg)
# T_j   <-  tapply(cij, zj.pg, function(x) length(unique(x)))
# O_j   <-  rbeta(Kz,shape1 = beta+1,shape2 = n_j)
# up1   <-  n_j/beta
# p1    <-  up1/(up1+1)
# s_j   <-  rbinom(n = Kz,size = 1,prob = p1)
# beta  <-  rgamma(1,a_beta + sum(T_j-s_j),  
#                    b_beta - sum(log(O_j)))
# 



# Re
#zj <- reset(oldzj)
#omega[, unique(zj)] = omega[, unique(oldzj)]
if (sim%%(verbose.step*thinning) == 0) {
  cat(paste("Sampling iteration: ", sim, " out of ",NSIM*thinning + burn_in, "---", length(unique(zj)), "\n",
            sep = "" ))}