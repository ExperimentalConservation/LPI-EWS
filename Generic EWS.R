###########################################################################################
##I wonder if including "consistent" change in a leading indicator, rather than just a single change event might do better? 
## e.g., only include an early warning signal if there are several concurrent singlas, rather than just a single one as we do at the moment?





###########################################################################################
rm(list=ls())
eps <- 1e-7 ## finite difference interval
library(data.table);library(reshape2);library(ggplot2)
library(earlywarnings);library(mgcv);library(parallel)
library(pROC)
# source('~/Dropbox/LPI EWS/R code/Generalised code - composite ews.R', chdir = TRUE)
# source('~/Dropbox/LPI EWS/R code/no plot ews.R', chdir = TRUE)
source('~/Desktop/LivingPlanet index early warnings/RCode/Github/LPI-EWS/Generalised code - composite ews.R', chdir = TRUE)
source('~/Desktop/LivingPlanet index early warnings/RCode/Github/LPI-EWS/no plot ews.R', chdir = TRUE)
###########################################################################################
##read in data
### data maipulation LPI
ddd<-read.csv("~/Dropbox/LPI EWS/Data/LPI_NA.csv")
length(ddd)
head(ddd, 2)

##different types of collected data
unique(ddd$Data_type)

exclude<-c("Index", "Measure per unit effort", "Estimate", "Proxy", "Occupancy", "Unknown")
dd1<-ddd[!(ddd$Data_type%in%exclude),]

unique(dd1$Data_type)

melt.ddd<-melt(dd1, id=c(1:49, 115:125));head(melt.ddd)

chrs<-strsplit(as.character(melt.ddd$variable), split="Y")

yrs<-do.call("rbind", mclapply(chrs, function(x){
	return(as.numeric(x[2]))
}, mc.cores=6))

melt.ddd$year<-yrs
dd<-melt.ddd

rm(ddd);rm(chrs);rm(melt.ddd);rm(yrs)

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

str(dd)
##split the data by the ID of each time series
split.dd<-split(dd, list(dd$ID))

##remove any time series that isnt at least 20 time points long
long.enough<-mclapply(split.dd, function(inc.na){
		## remove the na values
		y<-inc.na[!is.na(inc.na$value),]
		##set minimum time series length
		min.time.series<-20
		##if the time series is long enough, keep it!
		if(length(y[,1])>=min.time.series){
			return(y)
		}
}, mc.cores=4)

##number of time series we are left with
length(unique(rbindlist(long.enough)$ID))

results<-NULL
full.time.series<-NULL
ccntr<-0

#for(i in 1:length(long.enough)){
	#i=4968
	#i=29
	#i=563
full.time.series <-mclapply(long.enough, function(y){
    ##check a specific time series
    ##y<-subset(rbindlist(long.enough), ID==472)
		#=======================================================================================================================
		##y<-long.enough[[i]];print((i/length(long.enough)*100))
		if(length(y[,1])>0){
		#=======================================================================================================================
		##if there is a reasonable amount of variation in the time series data
		if(!length(which(diff(y$value)==0))>(length(y[,1])/2.8)){
			
		##the number of years of data, remove any NAs
		time.series.length<-length(na.omit(y$year))
		
		#fit a gam to the data	
		time.series.tippings<-gam_smoothing(y$year, y$value, -1)#round(length(y$year)/5))#
		
		##if the time series doenst just show a single trend (e.g. inc or dec)
		if(length(which(diff(time.series.tippings$sign)!=0))){
		
		##plot the data, and the GAM fit
		
		#ggplot(y, aes(x=year, y=value))+geom_line()+theme_classic()+stat_smooth(method="gam", formula=y~s(x, k=-1))
		# plot(y$year, y$value, type="l")
		# points(time.series.tippings$years, time.series.tippings$fit, type="l", lwd=2)
		
		#=======================================================================================================================
		##check if the data needs interpolating (is there any missing data)
		if(max(diff(y$year))>1){interp<-TRUE}else{interp<-FALSE}
		
		##take the time and abundance data
		mat.dat<-matrix(c(y$year, y$value), ncol=2)

		#calculate the generic ews using the dakos package	
		ews.res<-no.plot.ews(mat.dat, detrending="gaussian", interpolate=interp, winsize=25)
		##remove cv as it shouldnt be calculated if data is detrended
		ews.res$cv<-NULL
		
		#=======================================================================================================================
		##calculate the composite EWS, based on drake and griffen (every combination of 1:5)
		#=======================================================================================================================
		
		##blank object to save into
		comps.res<-NULL
		
		##a counter
		comp.cntr<-0
		
		##the indicators to include in the composite ews
		indicators<-c("cv", "acf", "ar1", "dr", "rr")
		
		##loop through the indicators to make every combination of 1:5
		for(h in 1:length(indicators)){
		inds<-data.frame(indicators)
			if(h!=1){for(o in 2:h){
				inds<-cbind(inds, indicators)}
				##all possible combinations of single, double etc
				combs.all<-expand.grid(inds)
				combs.all$code <- apply(combs.all[ ,1:length(combs.all)], 1, function(x) paste(sort(x), collapse = " + " ))
				combs<-as.data.frame(do.call("rbind", lapply(split(combs.all, list(combs.all$code)), function(j){
					j$code<-NULL
					return(j[1,])
				})))
				rownames(combs)<-NULL
			}else{combs<-inds}##combs is then the data frame of the indicators to be included 
			
			##then loop through the data frame of composite indicators, and for each one...
			for(m in 1:length(combs[,1])){
				
				##exclude any combinatinos of indicators where one indicators is included multiple times
				if(length(unique(unlist(combs[m,])))==h){
					##add to the counter
					comp.cntr<-comp.cntr+1					
					
					##run the composite EWS. calculates when the threshold has been calculated. Plot or not?
					comp<-composite_ews(mat.dat, indicators=c(as.character(unlist(combs[m,]))), 2, F)
					
					##make any points where the threshold is not crossed a 0 rather than an NA
					comp$threshold.crossed[is.na(comp$threshold.crossed)]<-0
					##save the output
					comps.res[[comp.cntr]]<-comp
					}
			}
		}	
    
    ##using the GAM approach, calculate the changes in direction of the ews 
    comp.inds <-mclapply(comps.res, function(ll){
      if(length(ll[,1])>0){
      return(cbind(ll, gam_smoothing(ll$time, ll$metric.score, -1)))#round(length(ll$time)/4)
    }})

    #=======================================================================================================================
	##for the dakos generic results, apply the gam techbique too
	ews.sigs<-NULL
	for(ii in 2:length(ews.res)){
		##ii=2
		ews.singals<-gam_smoothing(ews.res$timeindex, ews.res[,ii], -1)
						
		##add in a blank data frame of NAs to fill in the part of the data missed by the rolling window approach
		##(but not by the drake method)

		blank.fill<-as.data.frame(matrix(NA, nrow(subset(time.series.tippings, years<min(ews.singals$years))), ncol=ncol(ews.singals)))
		colnames(blank.fill)<-colnames(ews.singals)
		blank.fill$years<-subset(time.series.tippings, years<min(ews.singals$years))$years		
		ews.singals <-rbind(blank.fill, ews.singals)
		
		##add in the name of the leading indicator
		ews.singals$leading.indicator<-names(ews.res[ii])
		
		ews.sigs[[ii]]<-ews.singals
		#rm(ews.singals)
	}
	
	##subset the time series data (that shows whether the population is increasing, decreasing, etc) to only include data that
	##time.series.tippings<-subset(time.series.tippings, years>=min(ews.sigs[[2]]$years))

	## add the results from the gam tipping point to the data describing the population dynamics
	for(ll in 2:length(ews.sigs)){
		#ll=2
		datr<-ews.sigs[[ll]]
		time.series.tippings[paste(datr$leading.indicator[1])]<-datr$sign
	}
	
	##do the same for the composite indicators, pasiting in both the threshold approach results, and the gam results into the same data set as the above
	##time.series.tippings will then include the population dynamics, and all the ews results
	for(gg in 1:length(comp.inds)){
		#gg=1
		datr.2<-comp.inds[[gg]]
		if(length(datr.2[,1])>0){
			datr.3<-subset(datr.2, time>=min(ews.sigs[[2]]$years))
				
			time.series.tippings[paste("comb",paste(datr.2$metric.code[1]), sep=".")]<-datr.3$threshold.crossed[match(time.series.tippings$years, datr.3$time)]
			time.series.tippings[paste("gam.comb",paste(datr.2$metric.code[1]), sep=".")]<-datr.3$sign[match(time.series.tippings$years, datr.3$time)]
				}
		}
		
		##chage the sign of the return rate. theory suggest it shoudl degrease prior to a bifurcation.
		##its easier here to analyse if it is increasing
		time.series.tippings$returnrate<-time.series.tippings$returnrate*-1

		##remove the first 5 time points from the time series data, as 
		#shrt.tms<-time.series.tippings[5:length(time.series.tippings[,1]),]
		shrt.tms<-time.series.tippings
		
		##if there is enough data to continue the add a column to shrt.tms which denotes the different directions that the abundance data takes (inc, dec, etc)
		##this column will allow the data frame to be split
		if(length(shrt.tms[,1])>0){
			counter<-0;shrt.tms$split<-0			
			if(length(shrt.tms[,1])>1){for(jj in 2:length(shrt.tms$sign)){
				if(shrt.tms$sign[jj]==shrt.tms$sign[jj-1]){
					shrt.tms$split[jj]<-shrt.tms$split[jj-1]}else{
						counter<-counter+1
						shrt.tms$split[jj]<-counter
					}	
				}
			}
		
		##melt the data frame		
		melt.shrt.tms<-melt(shrt.tms, id=c(1:6, 76));head(melt.shrt.tms);tail(melt.shrt.tms)
		
		##then split it by the different ews
		split.shrt.tms<-split(melt.shrt.tms, list(melt.shrt.tms$variable))
		
		ddd<-rbindlist(lapply(split.shrt.tms, function(yy, time.series.tippings){
				##yy<-split.shrt.tms[[1]]
				
				yy<-na.omit(yy)
				yy$value[which(yy$value==-1)]<-0
				
				if(length(yy[,1])>0){

					TP<-0
					FP<-0
					TN<-0
					FN<-0
					
					window.length<-10
					
					##only consider it an EWS signal if the previous signal matched it, therefor 2:length
					for(tt in 2:(length(yy[,1])-window.length)){
						##tt<-2	
						sel.dat<-yy[(tt+1):((tt+1)+window.length),]
						
						# ##is there an ews? 1=yes, 0=no
						 val<-yy$value[tt]
						
						# #previous time step value
						#tm1.val<-yy$value[tt-1]
						
						#if(val==1){
					  #	val<-ifelse(val==tm1.val, 1, 0)
						#}
						
						##check whether true positive or false positive
						if(val==1){
							if(length(which(sel.dat$sign==-1))>0){
								TP<-TP+1}else{FP<-FP+1}
						}	
						##check whether true negative or false negative
						if(val==0){
							if(length(which(sel.dat$sign==-1))>0){
								FN<-FN+1}else{TN<-TN+1}
						}
					}
				#return(data.frame(variable=yy$variable[1], TP, FP, TN, FN))
				
				if(TP+FP+TN+FN>0){				
					roc.dat<-data.frame(variable=yy$variable[1], TP, FP, TN, FN, prop=sum(TP,TN)/sum(TP, FP, TN, FN))##TP.FP=c(rep(1, TP+TN), rep(0, FP+FN))
					return(roc.dat)
				}

				}}))					
				
		ccntr<-ccntr+1
		##for looped version
		#results[[ccntr]]<-cbind(time.series.length, ddd)	
		#full.time.series[[ccntr]]<-cbind(y[1,1:(length(y)-3)], time.series.length, ddd)
		declines<-length(which(time.series.tippings$sign==-1))
		return(cbind(as.data.frame(y)[1,c(1:(length(y)-3))], time.series.length, ddd, declines))
				}
			}
		}
	}
}, mc.cores=6)
#}

##number of analysed time series
ccntr
fin.res<-rbindlist(full.time.series);fin.res

#============================================================================================================
##add in some additional information about the
##is the predictor based on the gam approach?
lll<-unlist(lapply(strsplit(as.character(fin.res$variable), split=".",fixed=T), function(x){
	if(length(which(unlist(x)=="gam"))>0){return("Yes")}else{"No"}	
}))
fin.res$inc.gam<-lll

##is it a composite approach?
mmm<-unlist(lapply(strsplit(as.character(fin.res$variable), split=".",fixed=T), function(x){
	if(length(which(unlist(x)=="comb"))>0){return("Yes")}else{"No"}	
}))
fin.res$is.comb<-mmm

## add in the number of predictors in the method
nnn<-unlist(lapply(strsplit(as.character(fin.res$variable), split="[.+]"), function(x){
	##sx<-unlist(strsplit(as.character(fin.res$variable)[142765], "[.+]"))
	return(length(x[!(x %in% c("gam", "comb"))]))
}))
fin.res$n.preds<-nnn

save(fin.res, file="~/Dropbox/LPI EWS/Data/EWS results.RData")

#============================================================================================================

vars<-unique(fin.res$variable)

areas<-NULL

for(oo in 1:length(vars)){

	##oo<-2

	t1<-subset(fin.res, variable==vars[oo])

	rocs<-roc(t1$TP.FP, t1$ID, plot=F, smooth=F)
	
	areas[[oo]]<-data.frame(auc(rocs), vars[oo])
	
	clr<-ifelse(t1$inc.gam[1]=="Yes", "Red", "Blue")
	
	plot(rocs, add=ifelse(oo==1, F, T), col=adjustcolor(clr, alpha=0.2))

}

areas<-rbindlist(areas)
areas[order(areas$auc.rocs, decreasing=T),]
#============================================================================================================
##calculate proportion sof fp to tp
fin.res
rev(sort(tapply(fin.res$TP.FP, list(fin.res$variable), function(x){sum(x)/length(x)})))

##some investigations

#============================================================================================================

melt.res<-as.data.table(melt(fin.res, id=c(1:62)))

melt.res

world<-map_data("world")

#============================================================================================================
## plot out the locations of the time series
uniques<-rbindlist(mclapply(split(fin.res, list(fin.res$ID)), function(x){
	if(length(x[,1])>0){return(x[1,])}
}))

length(uniques$ID)

p1 <- ggplot(world, aes(long,lat,group=group)) + geom_polygon(fill="white", col="black")+theme_bw()+xlab("")+ylab("")+geom_point(data=uniques, aes(x= Longitude, y= Latitude, col=Class, group=ID), size=4,alpha=1, size=2)+theme(legend.position="top", axis.ticks=element_blank(), axis.title.x=element_blank(),axis.title.y=element_blank(),axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank())##+facet_wrap(~Class, ncol=2)
pdf("~/Desktop/LivingPlanet index early warnings/plots/Maps 1.pdf", width=6, height=6)
p1
dev.off()
#============================================================================================================
##and the lengths of the time series 
p2<-ggplot(uniques, aes(x= time.series.length))+geom_histogram(fill="black", col="white")+theme_classic()+ylab("Number of time series\n")+xlab("\nLength of time series")
pdf("~/Desktop/LivingPlanet index early warnings/plots/Time series length.pdf", width=6, height=3)
p2
dev.off()
#============================================================================================================
##proportion of true post to false negs

sort(tapply(fin.res$TP, list(fin.res$variable), function(l)sum(l)/length(l))-tapply(fin.res$FP, list(fin.res$variable), function(l)sum(l)/length(l)))


p3<-ggplot(data=propotions, aes(x=Var1, y=value))+geom_bar(stat="identity")+theme_classic()+theme(axis.text.x =element_text(angle=90, vjust=0.5))+xlab("")+ylab("Proportion of time series\ndisplaying early warning signals\n")
pdf("~/Desktop/LivingPlanet index early warnings/plots/Proportions displaying EWS.pdf", width=8, height=5)
p3
dev.off()
#============================================================================================================

cats<-c("LC", "NT", "VU", 	"EN", "CR")

dd2<-subset(melt.res, Red_list_category %in% cats)

levels(dd2$Red_list_category)

dd2$Redlist<-factor(dd2$Red_list_category, levels=cats)

dd2<-subset(dd2, variable %in% c("comb.dr"))

propotions.redlist<-na.omit(melt(tapply(dd2$value, list(dd2$Redlist, dd2$variable), function(l)sum(l)/length(l))))

#dd2<-subset(dd2, variable!="cv" & variable!="kurt" & variable!="sk")

p4<-ggplot(propotions.redlist)+geom_bar(aes(x=Var1, y=value),col="black", fill="black", stat="identity")+theme_bw()+facet_wrap(~Var2)

pdf("~/Desktop/LivingPlanet index early warnings/plots/EWS by category.pdf", width=8, height=5)
p4
dev.off()


prop.redlist.only<-melt(tapply(dd2$value, list(dd2$Redlist), function(l)sum(l)/length(l)))

ggplot(prop.redlist.only)+geom_bar(aes(x=Var1, y=value),col="black", fill="black", stat="identity")+theme_classic()+xlab("Redlist category")+ylab("Proportion showing EWS")
