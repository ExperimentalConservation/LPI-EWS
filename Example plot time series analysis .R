rm(list=ls())
set.seed(3)
eps <- 1e-7 ## finite difference interval
library(ggplot2)
library(mgcv)
library(reshape2)
library(earlywarnings)
source('~/Desktop/LivingPlanet index early warnings/RCode/Github/LPI-EWS/Generalised code - composite ews.R', chdir = TRUE)
###########################################################################################
##functino to calculate the slopes, derivatives and ci's of the data, based on the Gam method from Burthe et al.:
is.there = function(x=0, L.CI, U.CI){
		pos<- ifelse(x<U.CI, 1, -1) 
		negs<-ifelse(x>L.CI, -1, 1)
		return(pos+negs)}
													##round(length(timeseries)/4)	
gam_smoothing<-function(years, timeseries,knots){
		if(length(which(timeseries<=0))==0){
		gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian(link="log"))}else{
		gam1<-gam(timeseries~s(as.vector(years), bs="cs", k=knots), family=gaussian)}
		time.series.fit<-predict(gam1, newdata=data.frame(years=years), type="response")
		
		X0<-predict(gam1, newdata=data.frame(years=years), type= "lpmatrix")
		X1<-predict(gam1, newdata=data.frame(years=years+eps), type= "lpmatrix")
		Xi<-(X1-X0)/eps
		df <- Xi%*%coef(gam1)              ## ith smooth derivative 
		df.sd <- rowSums(Xi%*%gam1$Vp*Xi)^.5 ## cheap diag(Xi%*%b$Vp%*%t(Xi))^.5
		#plot(years,df,type="l",ylim=range(c(df+2*df.sd,df-2*df.sd)))##plot 'em
		#lines(years,df+2*df.sd,lty=2);lines(years,df-2*df.sd,lty=2)
		splines<-data.frame(years=years,deriv=df,U.CI=df+2*df.sd, L.CI=df-2*df.sd)	
		splines$sign<-is.there(0, splines$L.CI, splines$U.CI)/2
		splines$fit<-time.series.fit
		return(splines)}
###########################################################################################
##make some data
y<-c(rep(100, 10), 110, 120, 140, 140, 140, 120, 100, 80, 60, 40, 20, 10)
y<-y+rnorm(length(y), 0, 10)
x<-seq(y)

##calculate EWS
	##Dakos
	ar1.gen<-generic_ews(y, winsize=30)$ar1;dev.off()
	ar1<-c(rep(NA, length(x)-length(ar1.gen)), ar1.gen)
	##composite
	mat.dat<-matrix(c(x,y),ncol=2)
	comp.ews<-composite_ews(mat.dat, indicators=c("cv", "ar1"), 2, F)
	comp.ews<-rbind(NA, comp.ews)	
	
time.series<-cbind(data.frame(x,y, ar1))

time.series<-cbind(time.series, 
				"gam"=gam_smoothing(time.series$x, time.series$y, -1)$sign,
				"fit"=gam_smoothing(time.series$x, time.series$y, -1)$fit)


pdf("~/Desktop/LivingPlanet index early warnings/plots/Example ews plots.pdf",width=4, height=9)
par(mfrow=c(4,1), mar=c(1,4,0.5,1))
##plot time series and gam smoothing
plot(time.series$x, time.series$y, type="n", xlab="", ylab="Abundance", xaxt="n", xlim=c(0, 22), ylim=c(0, 150))
polygon(x=c(10,15, 15, 10),y=c(-10, -10, 170, 170), col= adjustcolor("darkgoldenrod2 ",alpha=0.4), lty=0)
lines(time.series$x, time.series$y, type="b")
colour <- ifelse(time.series$gam<0,"red","blue")
segments(time.series$x[x-1],
		time.series$fit[x-1],
		time.series$x[x+1],
		time.series$fit[x+1],
		col=colour, lwd=2)
legend("topleft", 
		legend=c("Time series", "GAM fit - not declining", "GAM fit - declining"), 
		lty=1, 
		pch=c(1, NA, NA),
		col=c("black", "blue", "red"), bty="n")
legend("topright", legend="a.", bty="n", cex=1.5)

##plot the generic ews with gam smoothing
gam.ar1<-gam_smoothing(time.series$x[!is.na(time.series$ar1)], time.series$ar1[!is.na(time.series$ar1)], 5)
colour.ar1 <- ifelse(gam.ar1$sign>0,"red","blue")
plot(time.series$x, time.series$ar1, type="n", xlab="Time", ylab="ar1", xaxt="n", xlim=c(0, 22), ylim=c(-0.8, 1.3))
polygon(x=c(10,15, 15, 10),y=c(-10, -10, 170, 170), col= adjustcolor("darkgoldenrod2 ",alpha=0.4), lty=0)
lines(time.series$x, time.series$ar1, type="l")
segments(gam.ar1$years[x-1],
		gam.ar1$fit[x-1],
		gam.ar1$years[x+1],
		gam.ar1$fit[x+1],
		col=colour.ar1, lwd=2)
legend("topleft", legend=c("ar1", "GAM fit - not increasing", "GAM fit - increasing"), lty=1, col=c("black", "blue", "red"), bty="n")
legend("topright", legend="b.", bty="n", cex=1.5)

##plot the composite ews, threshold crossed version
plot(comp.ews$time, comp.ews$metric.score, type="n", col="skyblue3", ylim=c(0, 7), lwd=2, xaxt="n", ylab="Composite EWS", xlim=c(0, 22))
polygon(x=c(10,15, 15, 10),y=c(-10, -10, 170, 170), col= adjustcolor("darkgoldenrod2 ",alpha=0.4), lty=0)
lines(comp.ews$time, comp.ews$metric.score, type="l", col="skyblue3", lwd=2)
lines(comp.ews$time, comp.ews$rolling.mean, type="l", lwd=1.5)
lines(comp.ews$time, comp.ews$rolling.mean+(2*comp.ews$rolling.sd), type="l", lty="dashed", lwd=1.5)
lines(comp.ews$time[which(comp.ews$threshold.crossed==1)], comp.ews$metric.score[which(comp.ews$threshold.crossed==1)], type="p", pch=19, col="maroon")
legend("topleft", legend=c("Rolling mean", 
							"Composite metric score",
							"Threshold value",
							"Thewshold crossed"), 
							lty=c("solid","solid","dashed",NA),
							pch=c(NA, NA, NA, 19),
							col=c("black", "skyblue3", "black", "maroon"), bty="n")
legend("topright", legend="c.", bty="n", cex=1.5)

##plot the composite ews, gam version 
par(mar=c(4,4,0.5,1))
gam.comp<-gam_smoothing(comp.ews$time[!is.na(comp.ews$metric.score)], 
						comp.ews$metric.score[!is.na(comp.ews$metric.score)],
						6)
plot(comp.ews$time, comp.ews$metric.score, type="n", col="skyblue3", ylim=c(0, 7), lwd=2, xlab="Time", ylab="Composite EWS", xlim=c(0, 22))
polygon(x=c(10,15, 15, 10),y=c(-10, -10, 170, 170), col= adjustcolor("darkgoldenrod2 ",alpha=0.4), lty=0)
lines(comp.ews$time, comp.ews$metric.score, type="l", col="skyblue3", lwd=2)
colour.comp<-ifelse(gam.comp$sign>0,"red","blue")
segments(gam.comp $years[x-1],
		gam.comp $fit[x-1],
		gam.comp $years[x+1],
		gam.comp $fit[x+1],
		col=colour.comp, lwd=2)
legend("topleft", 
	legend=c("Composite metric score","GAM fit - not increasing", "GAM fit - increasing"), 
	lty=1,
	lwd=1.5,
	col=c("skyblue3", "blue", "red"), bty="n")
legend("topright", legend="d.", bty="n", cex=1.5)

dev.off()


