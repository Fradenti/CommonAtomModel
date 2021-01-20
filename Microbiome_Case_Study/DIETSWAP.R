# Source Functions --------------------------------------------------------
Rcpp::sourceCpp("Microbiome_Case_Study/DCAM_LS_MAIN.cpp")
source("Microbiome_Case_Study/DCAM_LS_MAIN.R")

# Load libraries ----------------------------------------------------------
library(phyloseq)
library(GGally)
library(patchwork)
# Data citation doi: 10.1038/ncomms7342
# remotes::install_github("microbiome/microbiome")
library(microbiome)
library(tidyverse)
library(plot.matrix)
library(ComplexHeatmap)
library(mcclust.ext)
library(network)

# Load Data + Preprocessing -----------------------------------------------
data(dietswap)
J_all         <- nsamples(dietswap)
N_all         <- ntaxa(dietswap)
dim(dietswap@otu_table)
plot(rowSums(dietswap@otu_table))
plot(colSums(dietswap@otu_table))
# Some taxa are null!
idx_null_taxas <- which(rowSums(dietswap@otu_table)==0)
# sono 7 Genus
idx_null_taxas
tax_names <- rownames(dietswap@otu_table)
sn        <- sample_names(dietswap)
tn        <- taxa_names(dietswap)
all(table(rownames(tax_names))==1)
# they all have different names:  ✔

# Abundance table in matrix
dietswap@otu_table
OT <- matrix(c(dietswap@otu_table),N_all,J_all,byrow = F)
image(log(OT+1))
# Remove the taxa that have never been observed:
OT1 <- OT[-idx_null_taxas,]
dim(OT1)
rowSums(OT1)
# plot_heatmap(dietswap)
id1 <- as_tibble(dietswap@sam_data) %>% filter(timepoint==1)  %>% mutate(sample=as.numeric(stringr::str_extract(sample,pattern = "\\d+"))) %>% dplyr::select(sample)
# Collect ID of subjects at time 1
id1 <- as.numeric(pull(id1))

# Reduction of the main table to smaller dataset
OT_r <- OT1
N_r  <- dim(OT_r)[1]
OT_r <- OT_r[,id1]
rowSums(OT_r)
# There are still null taxa in the new reduced dataset
id_null_taxa_reduced <- which(rowSums(OT_r)==0)
OT_r2 <- OT_r[-id_null_taxa_reduced,]
N_r  <- dim(OT_r2)[1]
dim(OT_r2)
OT_r2[42:44,]
heatmap(OT_r2[1:50,],
        Rowv = NA,Colv = NA,
        col=viridis::cividis(10))
Heatmap(log(OT_r2+1),col = viridis::viridis(1000),cluster_rows = F,cluster_columns = F,
        name = "Abundance Table\nlog(Z+1) scale \n",
        row_title = "OTU (Genus)",column_labels = 1:38,show_heatmap_legend = T,
        column_split = as.character(dietswap@sam_data$nationality[id1]))

Heatmap(log(OT_r2+1),col = grey.colors(100),cluster_rows = F,cluster_columns = F,
        name = "Abundance Table\nlog(Z+1) scale \n",
        row_title = "OTU (Genus)",column_labels = 1:38,show_heatmap_legend = T,
        column_split = as.character(dietswap@sam_data$nationality[id1]))


# Summarizing results from preprocessing ----------------------------------
# OTU Table time 1:
OT_r2
dim(OT_r2)
head(OT_r2)
# Covariates:
SAM_t1 <-dietswap@sam_data[id1,]
dim(SAM_t1)
# Genus tax names:
Names  <- dietswap@tax_table[,3]
Names <- Names[-idx_null_taxas]
Names <- Names[-id_null_taxa_reduced]
length(Names)


# Create the 2 histograms for figure 1 ------------------------------------------
Y <- (OT_r2)
colnames(Y) <- paste("Subject",1:38)
Y <- as_tibble(Y)
Y

S1 <- ggplot(Y)+geom_histogram(aes(x=`Subject 1`),fill=I("grey"),col=1,bins=50)+
  ggtitle("Subject 1 - Microbiome composition")+
  xlab("OTU counts")+ylab("Frequency")+theme_bw()
S1
# ggsave("mic_SUB1.pdf",width = 8,height = 6)
# ggsave("mic_SUB1.png",width = 8,height = 6)
# ggsave("mic_SUB1.tiff",width = 8,height = 6)

S14 <- ggplot(Y)+geom_histogram(aes(x=`Subject 14`),fill=I("grey"),col=1,bins=50)+
  ggtitle("Subject 14 - Microbiome composition")+
  xlab("OTU counts")+ylab("Frequency")+theme_bw()
S14
#ggsave("mic_SUB14.pdf",width = 8,height = 6)
#ggsave("mic_SUB14.png",width = 8,height = 6)
#ggsave("mic_SUB14.tiff",width = 8,height = 6)

S1+S14
#ggsave("mic_SUB1and14.pdf",width = 16,height = 6)
#ggsave("mic_SUB1and14.png",width = 16,height = 6)
#ggsave("mic_SUB1and14.tiff",width = 16,height = 6)


# Ready to run the model! ---------------------------------------------------
y_o1 <- c(OT_r2)
y_g1 <- rep(1:length(id1),rep(N_r,length(id1)))
# library size
udg1 <- rep(tapply(y_o1, y_g1, mean),rep(N_r,length(id1)))
plot(udg1)
prior1 <- list(
  # hyperparameters NIG
  m0 = mean(y_o1), 
  k0 = 1/var(y_o1), a0=3, b0=1,
  # hyperparameters alpha and beta
  a_alpha=3, b_alpha = 3,
  a_beta =3, b_beta  = 3)

NSIM     = 50000
burn_in  = 50000
thinning =1
set.seed (19922508)

R <-DCAM_LibSize(y_obser = y_o1,
                                y_group = y_g1,
                                K0 = length(id1),
                                L0 = N_r,
                                prior    = prior1,
                                NSIM     = NSIM,
                                burn_in  = burn_in,
                                thinning = thinning,
                                verbose = T,
                                verbose.step = 2,
                                warm.start = F,
                                User.defined.gammas = udg1)
#saveRDS(R,"DCAM_results_Dietswap_time1_50k_betaok.RDS")

# Output Analysis ---------------------------------------------------------
#
R <- readRDS("DCAM_results_Dietswap_time1_50k_betaok.RDS")
# Auxiliary function 
reset <- function(x){
  z <- y <- x
  sy <- unique(sort(y))
  for(i in 1:length(unique(x))){
    z[which(x==sy[i])] <- i
  }
  return(z)
}

# 1 Distributional clusters analysis -------------------------------------------------
ZM  <- t(apply(matrix(unlist(R$Z_j),NSIM,length(R$Z_j[[1]]),byrow = T),1,reset))
psm <- mcclust::comp.psm(ZM)
image(psm)
# Collect distributional clusters - best partition accordinf to VI
cl1 <- cl1V  <- mcclust.ext::minVI(psm, method = "greedy")
table(apply(ZM,1,max))
table(cl1$cl)
image(psm)
Heatmap(psm,col = viridis::viridis(100),
        cluster_rows = F,cluster_columns = F,
        name = "Posterior\nCoclustering\nProbability",
        row_title = "Subjects",column_labels = 1:38,
        show_heatmap_legend = T,
        column_split = paste0("DC-",as.character(cl1$cl)),
        row_split = as.character(cl1$cl))

Heatmap(psm,col = gray.colors(100),cluster_rows = F,cluster_columns = F,
        name = "Posterior\nCoclustering\nProbability",
        row_title = "Subjects",column_labels = 1:38,
        show_heatmap_legend = T,
        column_split = paste0("DC-",as.character(cl1$cl)),
        row_split = as.character(cl1$cl))

#saveRDS(cl1,"Distr_Clust_50k_betaok.RDS")
cl1 <- readRDS("Distr_Clust_50k_betaok.RDS")

# Descriptive results
table(cl1$cl)
mean(OT_r2[,cl1$cl==1]==0)
mean(OT_r2[,cl1$cl==2]==0)
mean(OT_r2[,cl1$cl==3]==0)
table(cl1$cl,SAM_t1$bmi_group)
t(table(cl1$cl,SAM_t1$sex))
t(table(cl1$cl,SAM_t1$nationality))
# Create boxplots
SAM_t1["cl"] <- cl1$cl
At <- as_tibble(SAM_t1)
At <- At %>% mutate("Mean"          = apply(OT_r2,2,mean),
                    "Median"        = apply(OT_r2,2,median),
                    "Std.Dev"       = apply(OT_r2,2,sd),
                    "% of zeros"    = apply(OT_r2,2,function(x) mean(x==0)),
                    "Range"         = apply(OT_r2,2,function(x) diff(range(x))),
                    "Shannon Index" = apply(OT_r2,2,vegan::diversity),
                    "Simpson Index" = apply(OT_r2,2,function(x) vegan::diversity(x,"simpson")),
                    "Skewness"      = apply(OT_r2,2,moments::skewness),
                    "Kurtosis"      = apply(OT_r2,2,moments::kurtosis))
X <- At %>% dplyr::select(Mean,Median,Range,
                   Std.Dev,Skewness,Kurtosis,
                   `% of zeros`,`Simpson Index`,`Shannon Index`,cl)
Xr <- reshape2::melt(X,"cl")
Xr <- Xr %>% mutate( cl2 = paste0("DC-",cl)  )
ggplot(Xr)+geom_boxplot(aes(x=cl2,y=value,group=cl2))+facet_wrap(~variable,scales="free")+
  theme_bw()+xlab("Distributional Clusters")+ylab("")+
  theme(strip.text.x = element_text(size = 10))
# ggsave("Distrib_boxplots_50k_betaok.png" ,width = 12,height = 15)
# ggsave("Distrib_boxplots_50k_betaok.tiff",width = 12,height = 15)
# ggsave("Distrib_boxplots_50k_betaok.pdf", width = 12,height = 15)


# Distributional clusters by OTU
cluster_per_obs <- rep(cl1$cl,rep(nrow(OT_r2),ncol(OT_r2)))
LL <- apply(OT_r2,2,function(x) cumsum(sort((x),decreasing = T)) /sum((x)) )
dim(LL)
LL2 <- reshape2::melt(LL)
LL2 <- cbind(LL2,cluster_per_obs)
LL2 <- as_tibble(LL2) %>% mutate(cluster_per_obs=paste0("DC-",cluster_per_obs))
ggplot(LL2)+
  geom_hline(yintercept = c(0,1),lty=3)+
  geom_line(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),lwd=.5,alpha=.2)+
  geom_point(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),lwd=1,alpha=.9)+
  theme_bw()+scale_color_manual("Distributional\nCluster",values = c(2,4,1))+ylim(0,1)+
  ggtitle("CRF of microbiome subpopulations")+xlab("Taxa sorted by count")+ylab("Cumulative Relative Frequncy")
#ggsave("Ecdf_Micro_50k_betaok.png" ,height = 5,width = 8)          
#ggsave("Ecdf_Micro_50k_betaok.tiff",height = 5,width = 8)          
#ggsave("Ecdf_Micro_50k_betaok.pdf" ,height = 5,width = 8)          


ggplot(LL2)+
  geom_hline(yintercept = c(0,1),lty=3)+
  geom_line(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),lwd=.5,alpha=.05)+
  geom_point(aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),lwd=1,alpha=.05)+
  theme_bw()+scale_color_manual("Distributional\nCluster",values = c(2,4,1))+ylim(0,1)+
  ggtitle("CRF of microbiome subpopulations\nSubjects 9, 21, and 30")+xlab("Taxa sorted by count")+ylab("Cumulative Relative Frequncy")+
  geom_point(data = LL2 %>% filter(Var2==9 | Var2 ==21 |Var2 ==30),
             aes(x=Var1,y=value,group=Var2,col=as.factor(cluster_per_obs)),lwd=1,alpha=1)
#ggsave("Ecdf_Micro_dubbiosiCoclustMat.png" ,height = 5,width = 8)          
#ggsave("Ecdf_Micro_dubbiosiCoclustMat.tiff",height = 5,width = 8)          
#ggsave("Ecdf_Micro_dubbiosiCoclustMat.pdf" ,height = 5,width = 8)          

# Descriptive results
G1 <- apply(OT_r2[,cl1$cl==1],2,which.max)
G2 <- apply(OT_r2[,cl1$cl==2],2,which.max)
plot(table(G1))
plot(table(G2))
plot(OT_r2[,cl1$cl==3],type="h")
Names[c(19,78)]
Names[c(86)]

plot(rowMeans(apply(OT_r2[,cl1$cl==1],2,function(x) x/sum(x))),type="h")
rowMeans(apply(OT_r2[,cl1$cl==1],2,function(x) x/sum(x)))[c(19,78)]
plot(rowMeans(apply(OT_r2[,cl1$cl==2],2,function(x) x/sum(x))),type="h")
plot(OT_r2[,cl1$cl==3]/sum(OT_r2[,cl1$cl==3]),type="h")
Names[c(58)]
(OT_r2[,cl1$cl==3]/sum(OT_r2[,cl1$cl==3]))[86]








# 2 Observational Clusters Analysis --------------------------------------------------
NSIM = 50000
CSI <- R$Csi_ij

ALL_micr_label <- array(NA,c(ncol(OT_r2),nrow(OT_r2),NSIM))
for(i in 1:NSIM){
  ALL_micr_label[,,i] <- matrix(CSI[[i]],nrow(OT_r2),ncol(OT_r2),byrow = T) 
}
ALL_micr_label[,,1]

# index1 = which(cl1V$cl==1)
# L1 <- length(index1)
# ALL_micr_label1 <- array(NA,c(L1,N_r,NSIM))
# for(i in 1:NSIM){
#   ALL_micr_label1[,,i] <- ALL_micr_label[index1,,i] 
# }
# index2 = which(cl1V$cl==2)
# L2 <- length(index2)
# ALL_micr_label2 <- array(NA,c(L2,N_r,NSIM))
J <- 38
N_r <- 119

# Compute the best partition for observational clusters
ZMic <- matrix(unlist(R$Csi_ij),NSIM,N_r*J,byrow = T)
ZMic <- t(apply(ZMic,1,reset))
psm <- mcclust::comp.psm(ZMic)
clmic  <- mcclust.ext::minVI(psm)
saveRDS(psm,"coclusteringmatrix_microb_observlevel.RDS")
saveRDS(clmic,"clustering_microb_observlevel.RDS")

psm <- readRDS("coclusteringmatrix_microb_observlevel.RDS")
clmic <- readRDS("clustering_microb_observlevel.RDS")



# Descriptive results
plot(R$y_obser~clmic$cl)
tapply(R$y_obser, clmic$cl, median)
SS <- sort(tapply(R$y_obser, clmic$cl, median),index=T)
L  <-  clmic$cl
tapply(c(y_o1), clmic$cl, median)
ordered <- as.numeric(factor(L,levels = SS$ix))
plot(L,ordered)
boxplot(R$y_obser~ordered)

MM <- matrix(ordered,N_r,38)
meds <- tapply(c(y_o1), ordered, median)
dim(OT_r2)
OBSCmat <- matrix(ordered,119,38)


Heatmap(log(OT_r2+1),col = viridis::viridis(100),cluster_rows = F,cluster_columns = F,
        name = "Posterior\nCoclustering\nProbability",
        row_title = "Subjects",column_labels = 1:38,
        show_heatmap_legend = T,column_split = paste0("DC-",as.character(cl1$cl)))
Heatmap(OBSCmat,col = viridis::viridis(8),cluster_rows = F,cluster_columns = F,
        name = "Posterior\nCoclustering\nProbability",
        row_title = "Subjects",column_labels = 1:38,
        show_heatmap_legend = T,
        column_split = paste0("DC-",as.character(cl1$cl)))

# Creation of Abundance Classes
A   <- as_tibble(cbind(c(y_o1),ordered))
to  <- table(ordered)
lab <- c(rep("Low",sum(to[1:3])),
         rep("Medium",sum(to[4:6])),
         rep("High",sum(to[7:8])  ))
P <- rep(1:8,table(ordered))
Aa <- tibble(V1=c(y_o1),OR=factor(ordered))


Aa <- Aa %>% mutate(LAB = 
                      factor(ifelse(OR==1|OR==2|OR==3,"Low",
                                    ifelse(OR==4|OR==5,"Medium","High"))
                             ,levels=c("Low","Medium","High")))

ggplot(Aa)+
  geom_boxplot(aes(x=OR,y=V1,col=LAB))+
  #  geom_point(data=data.frame(),aes(x=1:8,y=meds))
  #facet_wrap(~LAB,scales = "free")+
  theme_bw()+scale_color_manual("Abundance\nCategory",values = c(1,2,4))+
  ylab("Abundance value")+xlab("Observational Cluster")+
  geom_text(data=data.frame(), aes(x=as.numeric(names(meds)), 
                                   y=meds+c(400,400,500,600,900,1100,500,500), 
                                   label=meds), 
            col='black', size=4)+
  geom_vline(xintercept = c(3.5,5.5),lty=3)+
  ggtitle("OTU abundance stratified by Observational Cluster")
# ggsave("Boxplots_OC_50k_betao.png" ,height = 5,width = 8)          
# ggsave("Boxplots_OC_50k_betao.tiff",height = 5,width = 8)          
# ggsave("Boxplots_OC_50k_betao.pdf" ,height = 5,width = 8)          


MM <- ifelse(MM==1|MM==2|MM==3,"Low",
             ifelse(MM==4|MM==5,"Medium","High"))
dim(MM)       
MM2 <- apply(MM,1,function(x) table(factor(x,levels = c("Low","Medium","High"))))
MM2 <- apply(MM2, 2, function(x) x/38)

Names2 <- stringr::str_remove(Names, " et rel.")
Names2 <- stringr::str_remove(Names2, " at rel.")

Heatmap(MM2[,1:40],col = viridis::viridis(8),cluster_rows = F,cluster_columns = F,
        name = "Percentage of\nAbundance class",
        row_title = "Abundance class",column_labels = Names2[1:40],
        show_heatmap_legend = T)
Heatmap(MM2[,41:80],col = viridis::viridis(8),cluster_rows = F,cluster_columns = F,
        name = "Percentage of\nAbundance class",
        row_title = "Abundance class",column_labels = Names2[41:80],
        show_heatmap_legend = T)
Heatmap(MM2[,81:119],col = viridis::viridis(8),cluster_rows = F,cluster_columns = F,
        name = "Percentage of\nAbundance class",
        row_title = "Abundance class",column_labels = Names2[81:119],
        show_heatmap_legend = T)




# 3 Co-occurrence Network analysis-----------------------------------------------------------------
# If needed, load the best distributional and observational partitions
# cl1   <- readRDS("Distr_Clust.RDS")
# clmic <- readRDS("OBScluster.RDS")

# create co-occurrence clusters
SS      <- sort(tapply(y_o1, clmic$cl, median),index=T)
L       <-  clmic$cl
tapply(c(y_o1), clmic$cl, median)
ordered <- as.numeric(factor(L,levels = SS$ix))
plot(L,ordered)
boxplot(y_o1~ordered)

MM <- matrix(ordered,N_r,38)
meds <- tapply(c(y_o1), ordered, median)
dim(OT_r2)
OBSCmat <- matrix(ordered,119,38)

MM <- ifelse(MM==1|MM==2|MM==3,"Low",
             ifelse(MM==4|MM==5,"Medium","High"))
dim(MM)       
MMnum <- ifelse(MM=="Low",1,ifelse(MM=="Medium",2,3))
cooccurrence_cluster <- comp.psm(t((MMnum)))
image(cooccurrence_cluster)
Names2 <- stringr::str_remove(Names, " et rel.")


cooccurrence_cluster1 <- comp.psm(t((MMnum[,cl1$cl==1])))
image(cooccurrence_cluster1)
image(cooccurrence_cluster1>.5)
cooccurrence_cluster2 <- comp.psm(t((MMnum[,cl1$cl==2])))
image(cooccurrence_cluster2)
image(cooccurrence_cluster2>.5)



name <- Names2
MM2 <- apply(MM,1,function(x) table(factor(x,levels = c("Low","Medium","High"))))

idx_low  <- which(apply(MMnum,1,median)<=1)#which(MM2[1,]==38)
MM3      <- MM2[,-idx_low]
MMnumred <- MMnum[-idx_low,-idx_low]
COOC_red <- cooccurrence_cluster[-idx_low,-idx_low]
image(COOC_red)
image(MMnumred)
Names_red <- Names2[-idx_low]
Names_red <- stringr::str_replace(Names_red," ","\n")


net = network(COOC_red>.5, 
              matrix.type = "adjacency", directed = FALSE)


# network plot
# GGally::ggnet2(net = net,
#        alpha = 1, mode = "circle",#, layout.par = list(cell.jitter = 1),
#        color = col1,
#        palette = viridis::magma(12),node.size = 10,
#        #node.color = med_clust,
#        label.size = 5,label=Names_red,
#        node.alpha = .8,
#        edge.alpha = .8,legend.position = "bottom")
# # for more inspiration, look at
# https://kateto.net/network-visualization

idx_low  <- which(apply(MMnum,1,median)<=1)#which(MM2[1,]==38)
cooccurrence_cluster1 <- comp.psm(t((MMnum[-idx_low,cl1$cl==1])))
image(cooccurrence_cluster1)
image(cooccurrence_cluster1>.5)
net1 = network(cooccurrence_cluster1>.5, 
               matrix.type = "adjacency", directed = FALSE)

cooccurrence_cluster2 <- comp.psm(t((MMnum[-idx_low,cl1$cl==2])))
image(cooccurrence_cluster2)
image(cooccurrence_cluster2>.5)
net2 = network(cooccurrence_cluster2>.5, 
               matrix.type = "adjacency", directed = FALSE)

col1 <- apply(MM[-idx_low,cl1$cl==1],1,function(x) which.max(table(factor(x,levels = c("Low","Medium","High")))))
col2 <- apply(MM[-idx_low,cl1$cl==2],1,function(x) which.max(table(factor(x,levels = c("Low","Medium","High")))))

col1 <- ifelse(col1==1,"Low",ifelse(col1==2,"Medium","High"))
col2 <- ifelse(col2==1,"Low",ifelse(col2==2,"Medium","High"))

col1 <- factor(col1,levels = c("Low","Medium","High"))
col2 <- factor(col2,levels = c("Low","Medium","High"))

n1 <- GGally::ggnet2(net = net1,
             alpha = 1, mode = "kamadakawai",#, layout.par = list(cell.jitter = 1),
             color = col1,palette = "Set1",
             node.size = 15,
             color.legend = "DC-1 - Modal Abundance Class",
             node.color = col1,
             label.size = 3,label=Names_red,
             node.alpha = .8,legend.position = "top",
             edge.alpha = .8)
n1
ggsave("N1.png" ,height = 10,width = 15)
ggsave("N1.tiff",height = 10,width = 15)
ggsave("N1.pdf" ,height = 10,width = 15)

n2 <- GGally::ggnet2(net = net2,
             alpha = 1, mode = "kamadakawai",#mode = "circle",#, layout.par = list(cell.jitter = 1),
             color = col2,
             color.legend = "DC-2 - Modal Abundance Class",
             palette = "Set1",
             node.size = 15,
             #node.color = med_clust,
             label.size = 4,label=Names_red,
             node.alpha = .8,
             edge.alpha = .8,legend.position = "top")#+theme(panel.background = element_rect(color = "grey"))
n2
ggsave("N2_50k_betaok.png" ,height = 10,width = 15)
ggsave("N2_50k_betaok.tiff",height = 10,width = 15)
ggsave("N2_50k_betaok.pdf" ,height = 10,width = 15)

#n1+n2

gridExtra::grid.arrange(n1,n2,ncol=2)
