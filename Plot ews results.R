rm(list=ls())
library(reshape2)
library(ggplot2)
library(data.table)
library(parallel)
load("~/Dropbox/LPI EWS/Data/EWS results.RData")

fin.res$TP.per.t<-fin.res$TP/fin.res$time.series.length

world<-map_data("world")

se<-function(x){sd(x)/sqrt(length(x))}

overall.means<-melt(rev(sort(tapply(fin.res$prop, list(fin.res$variable), mean))))
ses<-melt(rev(sort(tapply(fin.res$prop, list(fin.res$variable), se))))
overall.means$se<-ses$value[match(overall.means$Var1,ses$Var1)]

pdf("~/Desktop/LivingPlanet index early warnings/plots/Proportions displaying EWS.pdf", height=5, width=8)
ggplot(rbind(head(overall.means, 10), tail(overall.means, 10)), aes(x=Var1, y=value))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+geom_errorbar(aes(max=value+se, min=value-se),stat="identity", col="grey60", width=0.4)+ylab("Proporion of correct detections\n")+xlab("")+ylim(0,1)+ggtitle("Ten best and worse metrics across all time series")
dev.off()

##the best model
best<-overall.means$Var1[1]
best.res<-subset(fin.res, variable==best)
pdf("~/Desktop/LivingPlanet index early warnings/plots/Map of best metric.pdf", height=6, width=8)
ggplot(legend=FALSE) + geom_polygon(data=world, aes(x=long, y=lat,group=group))+
  geom_point(data=best.res,aes(x= Longitude ,y= Latitude, col=Class, size=TP/time.series.length), alpha=0.7)+
  ggtitle(paste("Number of true positive EWS per time series", "\n", "(Metric = ", best.res$variable[1], ")"))+
  ylab("")+xlab("")+
  theme(legend.position='bottom')+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())
  theme(axis.ticks = element_blank())
  dev.off()
  
##by red list category
cats<-c("LC", "NT", "VU", 	"EN", "CR")
best.res.cat<-subset(best.res, Red_list_category %in% cats)
best.res.cat$Redlist<-factor(best.res.cat$Red_list_category, levels=cats)
bar.by.cat<-melt(tapply(best.res.cat$TP.per.t, list(best.res.cat$Redlist), mean))
pdf("~/Desktop/LivingPlanet index early warnings/plots/EWS by category.pdf", height=4, width=6)
ggplot(bar.by.cat, aes(x= Var1, y=value))+geom_bar(stat='identity', fill="coral", alpha=0.5, col="black")+theme_classic()+geom_errorbar(aes(max=value+se(value), min=value-se(value)),stat="identity", col="black", width=0.2)+xlab("\nIUCN category")+ylab("Mean number of EWS per time step\n")+ggtitle(paste("Metric =", best.res.cat$variable[1]))
dev.off()

##time series length agains number of TP
ts.tp<-melt(tapply(best.res$TP.per.t, list(best.res$time.series.length), mean));names(ts.tp)<-c("Number of time points", "Mean number of TP per time step")
plot(ts.tp, type="l")
#ggplot(fin.res, aes(x=time.series.length, y=TP))+geom_bar(stat="identity")+theme_classic()

##

#========================================================================================================
##timie series with perfect predictions
perfects<-subset(best.res, prop==1)

perfect.score.ids<-perfects$ID

source('~/Desktop/LivingPlanet index early warnings/RCode/Github/LPI-EWS/read and sort LPI.R', chdir = TRUE)

##dd is the LPI
head(dd)

##the perfect score time series:
per.time<-dd[match(perfect.score.ids, dd$ID),]

per.time<-dd[dd$ID %in% perfect.score.ids,]

ggplot(data=per.time, aes(x=year, y=value))+geom_line()+theme_classic()+facet_wrap(~ID,scales="free")


#========================================================================================================
#distance from cities
library("ggmap")
library(maptools)
library(maps)
library(geosphere)
library(sp)

data(world.cities);head(world.cities)
world.cities<-subset(world.cities, pop>40000)
cities.mat<-data.frame(world.cities$long, world.cities$lat)

split.res<-split(fin.res, list(fin.res$ID))

res<-rbindlist(mclapply(split.res, function(x, cities.mat){
	x$distance.from.city<-sort(spDistsN1(as.matrix(cities.mat), c(x$Longitude[1], x$Latitude[1]), longlat=TRUE))[1]
	return(x)
}, cities.mat, mc.cores=6))

res$TP.plus.FN.by.t<-(res$TP+res$FN)/res$time.series.length

res<-subset(res, Red_list_category %in% cats)

res<-subset(res, Class!="Amphibia")
pdf("~/Desktop/LivingPlanet index early warnings/plots/Distance from cities.pdf", width=5, height=8)
ggplot(subset(res, distance.from.city<1000), aes(x=distance.from.city, y=TP.plus.FN.by.t, col=Red_list_category))+geom_point()+geom_smooth(method="lm", fill="white", size=1)+facet_wrap(~Class, scales="free", ncol=1)+theme_classic()+theme(legend.position="right")+scale_x_log10()+xlab("\nDistance from nearest population center\nwith >40,000 inhabitants")+ylab("Number of TP and FN early warning signals per time step")
dev.off()


res$cut<-cut(res$distance.from.city, 10)

barplot(tapply(res$TP.plus.FN.by.t, list(res$cut), mean))
#========================================================================================================
##do some collecting methods give better results?

kk<-na.omit(melt(tapply(fin.res$prop, list(fin.res$Data_type), function(l){
	
	return(mean(l))
	
})))

kk.se<-na.omit(melt(tapply(fin.res$prop, list(fin.res$Data_type), function(l){
	
	return(se(l))
	
})))

by.data.type<-as.data.frame(kk[order(kk$value),]);names(by.data.type)<-c("Data_type", "prop.correct");by.data.type$se<-kk.se[,2]

pdf("~/Desktop/LivingPlanet index early warnings/plots/Data type .pdf", width=7, height=4)
ggplot(by.data.type, aes(x=Data_type, y=prop.correct))+geom_bar(stat="identity", fill="purple", alpha=0.5, col="black")+theme_classic()+xlab("")+ylab("Proportion of correctly identified\n TP or TN per time step\n")+geom_errorbar(aes(max=prop.correct+se, min=prop.correct-se),stat="identity", col="black", width=0.2)
dev.off()


#========================================================================================================
##do some threats show more EWS than others?


