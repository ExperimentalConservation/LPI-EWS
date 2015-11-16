rm(list=ls())
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(grid)
library(pBrackets)
library(data.table)
library(parallel)
load("~/Dropbox/LPI EWS/Data/EWS results.RData")
head(fin.res)

##load map data
world<-map_data("world")

##calculate the F score across all time series for each metric
fin.res<-rbindlist(lapply(split(fin.res, fin.res$variable), function(o){
	o$Fscore<-2*(((sum(o$TP)/sum(o$TP+o$FP))*(sum(o$TP)/sum(o$TP+o$FN)))/((sum(o$TP)/sum(o$TP+o$FP))+(sum(o$TP)/sum(o$TP+o$FN))))
	return(o)
}))

##bin the LR/lc (outdated cat) in with the LC cat
fin.res$Red_list_category[which(fin.res$Red_list_category=="LR/lc")]<-"LC"

##calcualte the proportion of true positive per time step
fin.res$TP.per.t<-fin.res$TP/fin.res$time.series.length

##function to calculate standard error
se<-function(x){sd(x)/sqrt(length(x))}

##make a data frame of F scores and order then for plotting
f.scores<-unique(as.data.frame(fin.res)[,c("variable","Fscore")])
f.scores$variable<-reorder(f.scores$variable, -f.scores$Fscore)
##make some coloyrs for the graphic 
cols.blues <- colorRampPalette(tail(brewer.pal(9, "Blues"), 6))
myPal.blues <- cols.blues(length(f.scores[,1]))
##plot it out
pdf("~/Desktop/LivingPlanet index early warnings/plots/Fscores.pdf", height=5, width=10)
ggplot(f.scores, aes(x=variable, y=Fscore))+geom_bar(aes(fill=variable),stat="identity")+theme_classic()+theme(axis.text.x = element_text(angle = 90, hjust = 1, size=7))+scale_fill_manual(values = rev(myPal.blues))+theme(legend.position="none")+xlab("")+ylab("F score\n")
dev.off()

##calculate the overall means, extract the best and worst 10, used for plotting only
overall.means<-melt(rev(sort(tapply(fin.res$prop, list(fin.res$variable), mean))))
top.bottom.10<-rbind(head(overall.means, 10), tail(overall.means, 10))$Var1
##plot out
pdf("~/Desktop/LivingPlanet index early warnings/plots/Proportions displaying EWS.pdf", height=5, width=8)
p1.1<-ggplot(rbind(head(overall.means, 10), tail(overall.means, 10)), aes(x=Var1, y=value))+geom_bar(stat="identity", fill="white")+theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+ylab("Proporion of correct detections\n")+xlab("")+ylim(0,1)+ggtitle("Ten best and worse metrics across all time series")
p1.1+geom_boxplot(data=subset(fin.res, variable%in%top.bottom.10), aes(x=variable, y=prop, fill=variable), col="black", alpha=0.3, notch=T)+theme(legend.position="non")
dev.off()


##the best model
best<-f.scores$variable[which(f.scores$Fscore==max(f.scores$Fscore))]
best.res<-subset(fin.res, variable==best)

##plot out the spatial data of the best model, only for time series with specific lat/long
pdf("~/Desktop/LivingPlanet index early warnings/plots/Map of best metric.pdf", height=6, width=8)
ggplot(legend=FALSE) + geom_polygon(data=world, aes(x=long, y=lat,group=group), fill="grey80")+
  geom_point(data=subset(best.res,Specific_location==1),aes(x= Longitude ,y= Latitude, col=Class, size=TP/time.series.length), alpha=0.5)+
  ggtitle(paste("Proportion of true positive EWS per time series", "\n", "(Metric = ", best.res$variable[1], ")"))+
  ylab("")+xlab("")+
  theme(legend.position='bottom')+
  theme(panel.background = element_blank()) +
  theme(panel.grid.major = element_blank()) +
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())+
  theme(axis.ticks = element_blank())+ coord_equal(ratio=1)+
  scale_size_continuous(range = c(3,8))+theme(plot.title = element_text(size = 10))
  dev.off()
  
##calculate and plot out the distribution of time series lengths, split by data type
lengths<-unique(as.data.frame(fin.res)[,c("ID", "time.series.length", "Data_type")])
n.time.series<-length(lengths[,1])
pdf("~/Desktop/LivingPlanet index early warnings/plots/Time series length.pdf", height=4, width=6)
	ggplot(lengths, aes(x=time.series.length, fill=Data_type))+geom_histogram(col="white")+theme_classic()+scale_fill_brewer(palette="Set1")+theme(legend.position="top")+xlab("\nTime series length")+ylab("Frequency\n")+annotate("text", x=55, y=65, label=paste("Number of time series:", n.time.series), size=3.5)
dev.off()
    
##plot the mean number of EWS by red list category. Using data from the best metric only
cats<-c("NE", "DD", "LC", "NT", "VU", 	"EN", "CR")
best.res.cat<-subset(best.res, Red_list_category %in% cats)
best.res.cat$Redlist<-factor(best.res.cat$Red_list_category, levels=cats)
bar.by.cat<-na.omit(melt(tapply(best.res.cat$TP.per.t, list(best.res.cat$Redlist), mean)))
declines.by.cat<-na.omit(melt(tapply(best.res.cat$declines/best.res.cat$time.series.length, list(best.res.cat$Redlist), mean)))
bar.by.cat$dec.per.t<-declines.by.cat$value

pdf("~/Desktop/LivingPlanet index early warnings/plots/EWS by category.pdf", height=4, width=6)
p22<-ggplot(bar.by.cat, aes(x= Var1, y=value, fill=Var1))+geom_bar(stat='identity', alpha=0.5, col="black")+theme_classic()+geom_errorbar(aes(max=value+se(value), min=value-se(value)),stat="identity", col="black", width=0.2)+xlab("\nIUCN category")+ylab("Mean number of EWS per time step\n")+ggtitle(paste("Metric =", best.res.cat$variable[1], "\n"))+scale_fill_brewer(palette="Set2")+geom_vline(x=2.5, col="black")+theme(legend.position="none")
p22+annotate("text",x=bar.by.cat$Var1, label=c(round(bar.by.cat$dec.per.t, 2)), y=0.02)
dev.off()

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
	##x<-split.res[[1]]
	dists<-spDistsN1(as.matrix(cities.mat), c(x$Longitude[1], x$Latitude[1]), longlat=TRUE)
	#top.10.cities<-sort(dists)[5]
	lp<-cbind(world.cities, dists)
	lp<-subset(lp, dists<2000)
	lp$connect<-lp$pop/lp$dists^2
	x$connectance<-sum(lp$connect)
	x$distance.from.city<-sort(spDistsN1(as.matrix(cities.mat), c(x$Longitude[1], x$Latitude[1]), longlat=TRUE))[1]
	return(x)
}, cities.mat, mc.cores=6))

res$TP.plus.FN.by.t<-(res$TP+res$FN)/res$time.series.length

res<-subset(res, Red_list_category %in% cats)
res<-subset(res, Specific_location==1)
res<-subset(res, Class!="Amphibia")
res<-subset(res, Class!="Reptilia")
res<-subset(res, variable==best)

pdf("~/Desktop/LivingPlanet index early warnings/plots/Distance from cities.pdf", width=5, height=8)
ggplot(res, aes(x=(connectance+1), y=TP.per.t, col=Red_list_category))+geom_point(alpha=1, shape=1)+geom_smooth(method="lm",aes(col=Red_list_category),size=2, fill=NA)+facet_wrap(~Class, scales="free", ncol=1)+theme_classic()+theme(legend.position="right")+xlab(bquote('\nConnectance (Population size/'~Distance^2~')'))+ylab("Number of TP early warning signals per time step\n")+scale_x_log10()+ylim(0, 0.5)
dev.off()
#========================================================================================================
##do some collecting methods give better results?

kk<-na.omit(melt(tapply(fin.res$prop, list(fin.res$Data_type), function(l){
	return(mean(l))	
})))

kk.se<-na.omit(melt(tapply(fin.res$prop, list(fin.res$Data_type), function(l){
	return(se(l))	
})))


kk.lengths<-unique(as.data.frame(fin.res)[,c("ID", "Data_type", "time.series.length")])

k.lngs.mean<-na.omit(melt(tapply(kk.lengths$time.series.length, list(kk.lengths$Data_type), mean)))
k.lngs.se<-na.omit(melt(tapply(kk.lengths$time.series.length, list(kk.lengths$Data_type), se)))

by.data.type<-as.data.frame(kk[order(kk$value),]);names(by.data.type)<-c("Data_type", "prop.correct");by.data.type$se<-kk.se[,2]
by.data.type$k.lngs.mean<-k.lngs.mean$value
by.data.type$k.lngs.se<-k.lngs.se$value

pdf("~/Desktop/LivingPlanet index early warnings/plots/Data type .pdf", width=5.6, height=3)
ggplot(by.data.type, aes(x=Data_type, y=prop.correct, fill=Data_type))+geom_bar(stat="identity", alpha=0.5, col="black", width=0.8)+theme_classic()+xlab("")+ylab("Proportion of correctly identified\n TP or TN per time step\n")+geom_errorbar(aes(max=prop.correct+se, min=prop.correct-se),stat="identity", col="black", width=0.2)+scale_fill_brewer(palette="Set2")+theme(legend.position="none")+theme(axis.text.x = element_text(angle = 0))
dev.off()

#========================================================================================================
##proportional graphs
library(plyr)
head(fin.res)
fin.res<-as.data.frame(fin.res)
props<-unique(fin.res[,c("ID", "Class", "Data_type","System","Red_list_category", "Primary_threat", "Secondary_threat")])##, 
long.props<-melt(props, id.vars = 'ID')

fin.props<-na.omit(melt(tapply(long.props$ID, list(long.props$variable, long.props$value), length)))

fin.props<-rbindlist(lapply(split(fin.props, list(fin.props$Var1)), function(x){
	x$prop<-x$value/sum(x$value)
	#x<-x[order(x$prop, decreasing=T),]
	#x$prop<-factor(x$prop, levels=x$prop)
	return(x)
	}))

fin.props$Var2<-as.character(fin.props$Var2)
	
for(i in 1:length(fin.props$Var2)){
	#i=15
	fin.props$Var2[i]<-paste(strwrap(fin.props$Var2[i], width=20), collapse="\n")
}


fin.props<-rbindlist(lapply(split(fin.props, list(fin.props$Var1)), function(y){
	y<-y[order(y$prop, decreasing=T),]
	return(y)
	}))

fin.props<-ddply(fin.props, .(Var1), transform, pos=cumsum(prop)-(0.5*prop))

fin.props$Var2 <- reorder(fin.props$Var2, fin.props$prop)
fin.props$Var2 <- factor(fin.props$Var2, levels=rev(levels(fin.props$Var2)))

vars<-unique(fin.props$Var1)
titles<-c("Class", "Data type", "System", "Red list category", "Primary threat", "Secondary threat")
lss<-list()
for(k in 1:length(vars)){
#	k=4
	sub<-subset(fin.props, Var1==vars[k])
	p1<-ggplot(sub, aes(x=Var2))+geom_bar(aes(y=prop, fill=Var2),stat = "identity", col="black")+ scale_y_continuous("", breaks=NULL)+theme_bw()+xlab("")+scale_fill_discrete(name = "")+theme(legend.position="none", plot.margin=unit(c(-0.3,-2,-0.7,-2), "cm"))+theme(text = element_text(size=6))
	p1<-p1+  theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank())+
        geom_hline(yintercept = seq(0, round(max(sub$prop), 1), by = 0.1), colour = "grey90", size = 0.2)+geom_vline(xintercept = seq(length.out=length(unique(sub$Var2))), colour = "grey90", size = 0.2)+ coord_polar()+ggtitle(paste("\n",titles[k]))
	
	p1<-p1+ scale_fill_brewer(palette="Spectral")+theme(plot.title=element_text(size=8,  vjust=-0.8))
	
	gg_table<-ggplot_gtable(ggplot_build(p1))
	gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"
	
	lss[[k]]<-(gg_table)

}

#pdf("~/Desktop/LivingPlanet index early warnings/plots/Data breakdown polar.pdf", width=5, height=7)
grid.arrange(lss[[1]], lss[[2]], lss[[3]], lss[[4]],map2,  ncol=2)
#dev.off()

blankPanel<-grid.rect(gp=gpar(col="white"))


lat.long<-unique(fin.res[,c("ID", "Longitude","Latitude")])

map2<-ggplot(legend=FALSE) + geom_polygon(data=world, aes(x=long, y=lat,group=group))+
  geom_point(data= lat.long,aes(x= Longitude ,y= Latitude),col="blueviolet ", size=2, alpha=0.7)+
  ggtitle("")+
  ylab("")+xlab("")+
 	theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),
 	axis.ticks.x = element_blank(),axis.text.x = element_blank())+
 	theme(panel.background = element_blank()) +
 	theme(axis.ticks = NULL, plot.margin=unit(c(-2,-0,-5,-0), "cm"))+ coord_equal(ratio=1)



#pdf("~/Desktop/LivingPlanet index early warnings/plots/Data breakdown polar with map.pdf", width=7, height=8)
	grid.arrange(map2, lss[[1]], lss[[2]], lss[[3]], lss[[4]],lss[[5]],lss[[6]], layout_matrix = rbind(c(1,1,1),c(2,3,4), c(5,6,7)), heights=c(5, 2, 2))
#dev.off()

cols <- colorRampPalette(brewer.pal(12, "Set2"))
myPal <- cols(length(fin.props[,1]))

fin.props$uniques<-paste(as.character(fin.props$Var1), as.character(fin.props$Var2), sep=".")
fin.props$uniques.2<-reorder(fin.props$uniques, 1:length(fin.props[,1]))

p2<-ggplot(fin.props, aes(x=uniques.2))+geom_bar(aes(y=prop, fill=Var2, group=Var1),stat = "identity", col="black")+ scale_y_continuous("", breaks=NULL)+theme_bw()+xlab("")+theme(legend.position="none", plot.margin=unit(c(-0.3,-2,-0.7,-2), "cm"))+theme(text = element_text(size=10))
p2<-p2+coord_polar()+scale_x_discrete(labels=fin.props$Var2)
p2<-p2+theme(panel.border = element_blank(),
        legend.key = element_blank(),
        axis.ticks = element_blank(),
        axis.text.y = element_blank(),
        panel.grid  = element_blank())
p2<-p2+geom_hline(yintercept = seq(0, round(max(fin.props$prop), 1), by = 0.1), colour = "grey90", size = 0.3)+geom_vline(xintercept = seq(length.out=length((fin.props$Var2))), colour = "grey90", size = 0.3)
p2<-p2+geom_vline(xintercept = c(0.5, 6.5, 10.5, 13.5, 20.5, 27.5), colour = "grey50", size = 0.5)
p2<-p2+annotate("text",x=c(3.5, 8.5, 12, 17.5, 24, 30.5), y=0.5, label=c("Class", "Data type", "System", "Red list category", "Primary threat", "Secondary threat"), size=4)
p2<-p2+scale_fill_manual(values = myPal)
gg_table<-ggplot_gtable(ggplot_build(p2))
gg_table$layout$clip[gg_table$layout$name=="panel"] <- "off"

grid.arrange(gg_table)

pdf("~/Desktop/LivingPlanet index early warnings/plots/Data breakdown polar.pdf", width=8, height=7)
	grid.arrange(gg_table)
dev.off()