rm(list=ls())
library(reshape2)
library(ggplot2)
library(data.table)
load("~/Dropbox/LPI EWS/Data/EWS results.RData")
world<-map_data("world")

se<-function(x){sd(x)/sqrt(length(x))}

fin.res
ggplot(legend=FALSE) + geom_polygon(data=world, aes(x=long, y=lat,group=group))+
  geom_point(data=fin.res,aes(x= Longitude ,y= Latitude, col=Class, size=prop))+
  theme(legend.position='top')+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  theme(axis.ticks = element_blank())



overall.means<-melt(rev(sort(tapply(fin.res$prop, list(fin.res$variable), mean))))
ses<-melt(rev(sort(tapply(fin.res$prop, list(fin.res$variable), se))))
overall.means$se<-ses$value[match(overall.means$Var1,ses$Var1)]
ggplot(overall.means, aes(x=Var1, y=value))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_errorbar(aes(max=value+se, min=value-se),stat="identity", col="grey60", width=0.4)