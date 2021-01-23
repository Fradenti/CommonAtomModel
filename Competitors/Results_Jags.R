K      <- 1
psms   <- array(NA,c(12,12,30))
clusts <- matrix(NA,30,12)
nclu   <- numeric(30)
prop   <- numeric(30)
arands <- numeric(30)
frob   <- numeric(30)
# Results -----------------------------------------------------------------
reset <- function(x){
  z  <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}

K <- 1
ground <- rep(1:6,2)

Eval <- function(K, ground, type = "CAM"){
  right      <- length(unique(ground))
  DC_GT_PSM  <- mcclust::comp.psm(rbind(ground,ground))
  maxGTD     <- StatPerMeCo::Frobenius(
    matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
    matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))
  
for(i in 1:30){
  model     <- readRDS(paste0("Competitors/Results_CAM/out_",type,"_S1A_n",K*25,"_",i,".RDS"))
  Zeta      <- t(apply(model$model$zj[,,1],2,reset))
  prop[i]   <- mean(c(apply(Zeta,1,function(x) length(unique(x))))==right)
  
  psms[,,i]   <- mcclust::comp.psm(Zeta)
  clusts[i,]  <- mcclust.ext::minVI(psms[,,i],method = "greedy")$cl
  nclu[i]     <- length(unique(clusts[i,]))
  arands[i]   <- mcclust::arandi(clusts[i,],ground)
  frob[i]     <- StatPerMeCo::Frobenius(psms[,,i],DC_GT_PSM)/maxGTD

  cat(paste(i,"\n"))
  }

  return(list( Prop_correct = prop,
               PSM = psms,
               clustering_VI = clusts,
               nclu_VI = nclu,
               Arand = arands,
               Frob = frob))
  }  
  
  
  
res <- Eval(2,ground)
boxplot(res$Prop_correct)
boxplot(res$nclu_VI)
boxplot(res$Arand)
boxplot(res$Frob)
