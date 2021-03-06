
# Elapsed Time Analysis ---------------------------------------------------


list.files("RunningTimes/")


t1a <- readRDS("RunningTimes/Time1A.RDS")
t1b <- readRDS("RunningTimes/Time1B.RDS")
t2  <- readRDS("RunningTimes/Time2.RDS")
t3  <- readRDS("RunningTimes/Time3.RDS")
t3


library(tidyverse)
library(reshape2)

T1a <- cbind(melt(as_tibble(t1a)),Scenario="S 1A")
T1b <- cbind(melt(as_tibble(t1b)),Scenario="S 1B")
T2 <- cbind(melt(as_tibble(t2)),Scenario="S 2")
T3 <- cbind(melt(as_tibble(t3)),Scenario="S 3")

TTT <- rbind(T1a,T1b,T2,T3)

ttt <- TTT %>% group_by(variable,Scenario) %>% summarise(med=median(value))

TTT <- TTT %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
ttt <- ttt%>% mutate(variable= paste("Configuration",substr(variable,2,2)))
sl <- ggplot()+
  #geom_path(data=ttt,aes(x=variable,y=(med)/2,col=Scenario,group=Scenario),lwd=.5, position=position_dodge(.75))+
  geom_boxplot(data=TTT,aes(x=variable,y=(value)/2,col=Scenario), position=position_dodge(.75,preserve = "single")) + theme_bw()+
  ylab("Minutes") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Slice Sampler - 10k iterations")+
  theme(axis.text=element_text(size=12),
              axis.title=element_text(size=12),
              legend.text = element_text(size=12),
            legend.title = element_text(size=12))
sl
ggsave("RunningTimes/SLice.png",height = 10,width = 12)
ggsave("RunningTimes/SLice.pdf",height = 10,width = 12)




ggplot()+
  geom_path(data=ttt,aes(x=variable,y=(med)/2,col=Scenario,group=Scenario), position = position_dodge(.5))+
  geom_boxplot(data=TTT,aes(x=variable,y=(value)/2,col=Scenario), position=position_dodge(.75,preserve = "single")) + theme_bw()+
  ylab("log2(Minutes)") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Slice Sampler - 10k iterations")+
  coord_trans( y="log2")+theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=12),
                               legend.text = element_text(size=12),
                               legend.title = element_text(size=12))
ggsave("RunningTimes/Log2SLice.png",height = 10,width = 12)
ggsave("RunningTimes/Log2SLice.pdf",height = 10,width = 12)








tj1a_c <- readRDS("RunningTimes/Time_CAM_S1a.RDS")[[1]]
tj1a_n <- readRDS("RunningTimes/Time_NDP_S1a.RDS")[[1]]
tj1a_h <- readRDS("RunningTimes/Time_hNDP_S1a.RDS")[[1]]
tj2_c  <- readRDS("RunningTimes/Time_CAM_S2.RDS")[[1]]
tj2_n  <- readRDS("RunningTimes/Time_NDP_S2.RDS")[[1]]
tj2_h  <- readRDS("RunningTimes/Time_hNDP_S2.RDS")[[1]]





Tj1a_c  <- cbind(melt( as_tibble(tj1a_c)),Model="CAM")
Tj1a_n  <- cbind(melt( as_tibble(tj1a_n)),Model="NDP")
Tj1a_h  <-  cbind(melt(as_tibble(tj1a_h)),Model="nHDP")
Tj2_c   <-  cbind(melt(as_tibble(tj2_c )),Model="CAM")
Tj2_n   <-  cbind(melt(as_tibble(tj2_n )),Model="NDP")
Tj2_h   <-  cbind(melt(as_tibble(tj2_h )),Model="nHDP")



TTT1 <- rbind(Tj1a_c,Tj1a_n,Tj1a_h)

ttt1 <- TTT1 %>% group_by(variable,Model) %>% summarise(med=median(value))

TTT1 <- TTT1 %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
ttt1 <- ttt1 %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
S1J <- ggplot()+
#  geom_path(data=ttt1,aes(x=variable,y=(med),col=Model,group=Model),lwd=.5, position=position_dodge(.75,preserve = "single"))+
  geom_boxplot(data=TTT1,aes(x=variable,y=(value),col=Model), position=position_dodge(.5)) + theme_bw()+
  ylab("Minutes") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Jags - 10k iterations - S 1A")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))
S1J
ggsave("RunningTimes/Jags_S1.png",height = 10,width = 12)
ggsave("RunningTimes/Jags_S1.pdf",height = 10,width = 12)




ggplot()+
#  geom_path(data=ttt1,aes(x=variable,y=(med),col=Model,group=Model), position = position_dodge(.5))+
  geom_boxplot(data=TTT1,aes(x=variable,y=(value),col=Model), position=position_dodge(.5)) + theme_bw()+
  ylab("log2(Minutes)") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Jags - 10k iterations - S 1A")+
  coord_trans( y="log2")+theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=12),
                               legend.text = element_text(size=12),
                               legend.title = element_text(size=12))
ggsave("RunningTimes/Log2Jags_S1.png",height = 10,width = 12)
ggsave("RunningTimes/Log2Jags_S1.pdf",height = 10,width = 12)

















TTT1 <- rbind(Tj2_c,Tj2_n,Tj2_h)

ttt1 <- TTT1 %>% group_by(variable,Model) %>% summarise(med=median(value))

TTT1 <- TTT1 %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
ttt1 <- ttt1 %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
S2J<- ggplot()+
  #geom_path(data=ttt1,aes(x=variable,y=(med),col=Model,group=Model),lwd=.5, position=position_dodge(.5))+
  geom_boxplot(data=TTT1,aes(x=variable,y=(value),col=Model), position=position_dodge(.5)) + theme_bw()+
  ylab("Minutes") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Jags - 10k iterations - S 2")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))
S2J
ggsave("RunningTimes/Jags_S2.png",height = 10,width = 12)
ggsave("RunningTimes/Jags_S2.pdf",height = 10,width = 12)




ggplot()+
  geom_path(data=ttt1,aes(x=variable,y=(med),col=Model,group=Model), position = position_dodge(.5))+
  geom_boxplot(data=TTT1,aes(x=variable,y=(value),col=Model), position=position_dodge(.5)) + theme_bw()+
  ylab("log2(Minutes)") + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle("Elapsed time in minutes - Jags - 10k iterations - S 2")+
  coord_trans( y="log2")+theme(axis.text=element_text(size=12),
                               axis.title=element_text(size=12),
                               legend.text = element_text(size=12),
                               legend.title = element_text(size=12))
ggsave("RunningTimes/Log2Jags_S2.png",height = 10,width = 12)
ggsave("RunningTimes/Log2Jags_S2.pdf",height = 10,width = 12)


sl
library(patchwork)
S1J/S2J
ggsave("RunningTimes/Jags_S212.png",height = 12,width = 10)
ggsave("RunningTimes/Jags_S212.pdf",height = 12,width = 10)

sl/(S1J+S2J)
ggsave("RunningTimes/SliceJags_S212.png",height = 11,width = 13)
ggsave("RunningTimes/SliceJags_S212.pdf",height = 11,width = 13)






# time table --------------------------------------------------------------


Listone <- list()
Listone[[1]]<- readRDS("RunningTimes/Time1A.RDS")
Listone[[2]]<- readRDS("RunningTimes/Time1B.RDS")
Listone[[3]]<- readRDS("RunningTimes/Time2.RDS")
Listone[[4]]<- readRDS("RunningTimes/Time3.RDS")
Listone[[5]]<- readRDS("RunningTimes/Time_CAM_S1a.RDS")[[1]]
Listone[[6]]<- readRDS("RunningTimes/Time_NDP_S1a.RDS")[[1]]
Listone[[7]]<- readRDS("RunningTimes/Time_hNDP_S1a.RDS")[[1]]
Listone[[8]]<- readRDS("RunningTimes/Time_CAM_S2.RDS")[[1]]
Listone[[9]]<- readRDS("RunningTimes/Time_NDP_S2.RDS")[[1]]
Listone[[10]] <- readRDS("RunningTimes/Time_hNDP_S2.RDS")[[1]]


library(tidyverse)
X <- map(Listone,~ apply(.x,2,function(x)  c(mean(x),sd(x))  ))

XX <- rbind(
cbind(X[[1]],X[[2]]),
X[[3]],X[[4]]
)
round(XX,4)


Y <- rbind(cbind(X[[5]],X[[8]]),
           cbind(X[[6]],X[[9]]),
           cbind(X[[7]],X[[10]]))
round(Y           ,4)
