##################################################################################
### Process rainfall data from US to produce single daily rainfall time series ###
##################################################################################

#############
### SETUP ###
#############

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"
setwd(paste0(maindir,"Analyses/Venable"))

#################
### LOAD DATA ###
#################

# nc = national climate

require(data.table)

nc1 <- fread(
	paste(maindir, "Data/Venable/National_Climatic_Data/521606.csv",sep=""),
	header=T, colClasses=rep(c("character","numeric"),c(6,4))
	)
	# TSUN missing
nc2 <- fread(
	paste(maindir, "Data/Venable/National_Climatic_Data/521610.csv",sep=""),
	header=T, colClasses=rep(c("character","numeric"),c(6,5))
	)
nc3 <- fread(
	paste(maindir, "Data/Venable/National_Climatic_Data/521611.csv",sep=""),
	header=T, colClasses=rep(c("character","numeric"),c(6,5))
	)

nclist <- list(nc1,nc2,nc3)
nc <- rbindlist(nclist,fill=T)
nc <- subset(nc, select=c(STATION_NAME,LATITUDE,LONGITUDE,DATE,PRCP))
	# dropping STATION ELEVATION TSUN TMAX TMIN TOBS
	# nlevels(STATION)=213; nlevels(STATION_NAME)=200; so using STATION_NAME

remove(nc1,nc2,nc3,nclist)
gc()

### RENAME VARIABLES 
oldnames <- names(nc)
newnames <- c("station","lat","long","date","prcp")
setnames(nc,old=oldnames,new=newnames)

### FILL IN MISSING VALUES

nc$long[nc$long=="unknown"] <- NA
nc$lat[nc$lat=="unknown"] <- NA
nc$prcp[nc$prcp==-9999] <- NA

### CALCULATE DISTANCES OF STATIONS

statx <- with(nc,tapply(as.numeric(long),station,mean,na.rm=T))
staty <- with(nc,tapply(as.numeric(lat),station,mean,na.rm=T))
statlevels <- names(statx)

plot(statx,staty)
tumax <- -111.0025
tumay <- 32.22527777777778
points(tumax,tumay,pch="+",col="red") # desert laboratory
unipos <- which(statlevels=="TUCSON U OF A NUMBER 1 AZ US")
points(statx[unipos],staty[unipos],pch=16,col="blue")

require(geosphere)
statdists <- distVincentySphere(cbind(statx,staty), c(tumax,tumay))/1000
	# /1000 because dist given in metres
dd <- data.frame(station=statlevels,statdists=statdists)
	# distance data
dd <- na.omit(dd)
	# removed one station with unknown coordinates
dd[order(dd$statdists),]

### SUBSET TO ALL STATIONS WITHIN 10KM RADIUS 

stat10 <- with(dd, station[statdists<10])
nc10 <- subset(nc,station %in% stat10)

#####################
### PLOTTING DATA ###
#####################

require(plyr)
nc10$year <- as.numeric(substr(nc10$date, start=1, stop=4))
nc10$station <- as.factor(as.character(nc10$station))
runyears <- with(nc10,table(station,year))

firstyear <- with(nc10,tapply(year,station,min,na.rm=T))
nyears <- with(nc10,tapply(year,station,function(x) length(unique(x))))

nc10stats <- levels(as.factor(nc$station)) %in% levels(nc10$station)
keepstat <- names(which(firstyear < 2009 & nyears > 10))
	# Want stations that begin before 2009 and have more than 10 years of data

nstat <- length(unique(keepstat))

nc10ys <- as.data.frame(subset(nc10, station %in% keepstat))

nc10ys$year <- substr(nc10ys$date, start=1, stop=4)
	# extract year (first four characters)
nc10ys$station <- as.factor(as.character(nc10ys$station))
nc10ys <- subset(nc10ys,select=c("station","year","prcp"))

nc10ys <- ddply(nc10ys ,.(station,year),summarise,prcp=mean(prcp,na.rm=T))
nc10_ysmat <- daply(nc10ys , .(year,station), function(x) x$prcp)


### RAIN AND GEOGRAPHIC LOCATION FIGURE (OLD) 

pdf(paste0("weather_stations_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=13,height=7)

	par(mfrow=c(1,2))

	matplot(as.numeric(rownames(nc10_ysmat)),log(nc10_ysmat),type="l")
	legend("bottomleft",legend=keepstat,col=1:nstat,lty=1:nstat,bty="n")

	plot(statx[nc10stats],staty[nc10stats])

	points(tumax,tumay,pch="+",col="red") # desert laboratory
	points(statx[nc10stats][keepstat],staty[nc10stats][keepstat],pch=16,col=1:nstat)
	text(x=statx[nc10stats][keepstat],y=staty[nc10stats][keepstat],
		labels=round(statdists[nc10stats][levels(as.factor(nc$station))[nc10stats]%in% keepstat],1),
		pos=1,col=1:nstat)

dev.off()

###################
### AGGREGATION ###
###################

### AGGREGATE RAINFALL BY DATE

nc10d <- subset(nc10,station %in% keepstat)

nc10d <- nc10[,mean(prcp, na.rm=T),by=date]
	# national climate 10km averaged by date
setnames(nc10d,old=names(nc10d),new=c("date","prcp"))

### FILL IN MISSING ROWS (DATES)

# Although -9999 entered in climate data for some missing values, 
	# not done for all, so have to fill in manually

	# Missing months for University:
	# 1989: Aug, Oct
	# 2006: Nov
	# 2009: Jul, Nov, Dec
	# 2010: Jan (1/3)
	# 2012: Aug-Dec
	# ~6 important months missing overall (plus additional days)

dateseq <- seq(as.Date(min(nc10d$date),format="%Y%m%d"),
	as.Date(max(nc10d$date),format="%Y%m%d"), 
	"day")
	# creates continuous date sequence

dateseq <- format(dateseq, "%Y%m%d")
	# puts into nc date format

setkey(nc10d,date)
nc10d <- nc10d[CJ(dateseq)]

### CALCULATE TIMES

nc10d <- as.data.frame(nc10d)
nc10d$date <- strptime(nc10d$date,format="%Y%m%d")

#################################
### OUTPUT COMBINED DATA FILE ###
#################################

write.csv(nc10d,
	paste0("aggregate_rainfall_",format(Sys.Date(),"%d%b%Y"),".csv"),
	row.names=F)

#################################
### COMPARING DATA TO VENABLE ###
#################################

require(plyr)
nc10$year <- as.numeric(substr(nc10$date, start=1, stop=4))
nc10$station <- as.factor(as.character(nc10$station))
runyears <- with(nc10,table(station,year))

firstyear <- with(nc10,tapply(year,station,min,na.rm=T))
nyears <- with(nc10,tapply(year,station,function(x) length(unique(x))))

nc10stats <- levels(as.factor(nc$station)) %in% levels(nc10$station)
keepstat <- names(which(firstyear < 2009 & nyears > 10))
	# Want stations that begin before 2009 and have more than 10 years of data

nstat <- length(unique(keepstat))

nc10ys <- as.data.frame(subset(nc10, station %in% keepstat))

nc10ys$year <- substr(nc10ys$date, start=1, stop=4)
nc10ys$station <- as.factor(as.character(nc10ys$station))
nc10ys <- subset(nc10ys,select=c("station","year","prcp"))

nc10ys <- ddply(nc10ys ,.(station,year),summarise,prcp=mean(prcp,na.rm=T))
nc10_ysmat <- daply(nc10ys , .(year,station), function(x) x$prcp)


pdf(paste0("weather_stations_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=13,height=7)

par(mfrow=c(1,2))

matplot(as.numeric(rownames(nc10_ysmat)),log(nc10_ysmat),type="l")
legend("bottomleft",legend=keepstat,col=1:nstat,lty=1:nstat,bty="n")

plot(statx[nc10stats],staty[nc10stats])

points(tumax,tumay,pch="+",col="red") # desert laboratory
points(statx[nc10stats][keepstat],staty[nc10stats][keepstat],pch=16,col=1:nstat)
text(x=statx[nc10stats][keepstat],y=staty[nc10stats][keepstat],
	labels=round(statdists[nc10stats][levels(as.factor(nc$station))[nc10stats]%in% keepstat],1),
	pos=1,col=1:nstat)

dev.off()

### ANGERT DATA 
	# Doesn't name the source of this data (uni or desert laboratory)

ang <- read.csv(paste(maindir,
	"Data/Venable/angert_2009_figS1_b.csv",sep=""),
	header=T,skip=1)
	# Actually Kimball et al 2010

require("reshape2")
ang2 <- colsplit(ang[,1],",",names=c("x","y"))
ang3 <- ang2[order(ang2$x),]
ang3$year <- 1983:2007 

names(ang3)[names(ang3)=="y"] <- "precip"

### COMPARISON 

ncda <- subset(ncd,date$mon<=4,select=names(ncd)[names(ncd)!="date"])
	# a = april (4th month)

aa <- ddply(ncda, .(station,year), summarize,
	prcp = mean(prcp,na.rm=T)	# should we be omitting NAs here?
	) 
	# april aggregated 

plot(log(prcp)~year,data=aa,type="n")
for(i in 1:length(levels(aa$station))){
	lines(log(prcp)~year,data=subset(aa,station==levels(aa$station)[i]),col=i)
	}

aau <- subset(aa,station=="TUCSON U OF A NUMBER 1 AZ US")
	# u = unversity of arizona

auam <- merge(aau,ang3)
	# angert university april merge

aaa <- ddply(aa, .(year), summarize,
	prcp = mean(prcp,na.rm=T)	# should we be omitting NAs here?
	) 
	# april aggregated all

aaam <- merge(aaa,ang3)
	# angert all april merge

plot(prcp~year,data=aau,type="b")

par(mfrow=c(1,2))
plot(prcp~precip,data=auam)
with(auam,cor.test(prcp,precip))
plot(prcp~precip,data=aaam)
with(aaam,cor.test(prcp,precip))
	# tighter relationship for just university of arizona 
	# (but not much in it)
