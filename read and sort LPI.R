##read and sort LPI data

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
