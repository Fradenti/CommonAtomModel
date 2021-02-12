plotter <- function(x1,x2,x3,x4, Title = "???",yl="?",log=T){

T1a <- cbind(melt(as_tibble(x1)),Scenario="S 1A")
T1b <- cbind(melt(as_tibble(x2)),Scenario="S 1B")
T2 <- cbind(melt(as_tibble( x3)),Scenario="S 2")
T3 <- cbind(melt(as_tibble( x4)),Scenario="S 3")

TTT <- rbind(T1a,T1b,T2,T3)

#ttt <- TTT %>% group_by(variable,Scenario) %>% summarise(med=median(value))

TTT <- TTT %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
#ttt <- ttt%>% mutate(variable= paste("Configuration",substr(variable,2,2)))
p <- ggplot()+
#  geom_path(data=ttt,aes(x=variable,y=(med),col=Scenario,group=Scenario),lwd=.5, position=position_dodge(.5))+
  geom_boxplot(data=TTT,aes(x=variable,y=(value),col=Scenario), position=position_dodge(.75,preserve = "single"),
               width=.5) + theme_bw()+
  ylab(yl) + xlab("Sample Size")+
  scale_color_brewer(palette = "Dark2")+
  ggtitle(Title)+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12))

if(log){
  p <- p+  coord_trans( y="log2")+theme(axis.text=element_text(size=12),
                                        axis.title=element_text(size=12),
                                        legend.text = element_text(size=12),
                                        legend.title = element_text(size=12))
  }
p 
}










plotter_jags <- function(x1,x2,x3, Title = "???",yl="?",log=T){
  
  
  Tj1a_c  <- cbind(melt( as_tibble(x1)),Model="CAM")
  Tj1a_n  <- cbind(melt( as_tibble(x2)),Model="NDP")
  Tj1a_h  <-  cbind(melt(as_tibble(x3)),Model="nHDP")
  
  
  TTT1 <- rbind(Tj1a_c,Tj1a_n,Tj1a_h)
  
 # ttt <- TTT %>% group_by(variable,Model) %>% summarise(med=median(value))
  
  TTT1 <- TTT1 %>% mutate(variable= paste("Configuration",substr(variable,2,2)))
 # ttt <- ttt%>% mutate(variable= paste("Configuration",substr(variable,2,2)))
  p <- ggplot()+
    #  geom_path(data=ttt,aes(x=variable,y=(med),col=Model,group=Model),lwd=.5, position=position_dodge(.5))+
    geom_boxplot(data=TTT1,aes(x=variable,y=(value),col=Model), position=position_dodge(.75,preserve = "single"),
                 width=.5) + theme_bw()+
    ylab(yl) + xlab("Sample Size")+
    scale_color_brewer(palette = "Dark2")+
    ggtitle(Title)+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12),
          legend.text = element_text(size=12),
          legend.title = element_text(size=12))
  
  if(log){
    p <- p+  coord_trans( y="log2")+theme(axis.text=element_text(size=12),
                                          axis.title=element_text(size=12),
                                          legend.text = element_text(size=12),
                                          legend.title = element_text(size=12))
  }
  p 
}
