###############################################################
### Finding best times for open and close of growing season ###
###############################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
source("venable_figure_functions_31Jul2015.R")

cd <- read.csv("cd_20Jun2015.csv",header=T)

####################
### CLIMATE DATA ###
####################

# Need to construct joint model that finds right period given likelihood of
# two above models multiplied together

nc <- read.csv("aggregate_rainfall_03May2015.csv",header=T)
	# nc = national climate

nc$date <- strptime(nc$date,format="%Y-%m-%d")
	# leave as late as possible
nc$year <- nc$date$year + 1900 # actual year, not survey year
nc$month <- nc$date$mday
nc$yday <- nc$date$yday # 01 Jan = yday 0

### FUNCTION

seasdaycalc <- function(yday,year,iswinter){

	require(lubridate)
		# for leap year calculations

	isleap <- leap_year(year)

	seasday <- vector("numeric",length=length(yday))
	seasday[iswinter==T] <- ifelse(isleap[iswinter==T], 
		yday[iswinter==T]-365-1,  
		yday[iswinter==T]-364-1
		)
		# yday starts at zero, so leap -> 365, no leap -> 364
		# Need to add 1 because previous year starts at 0
	seasday[iswinter==F] <- yday[iswinter==F]

	return(seasday)

	}

### CENSUS

gystrip <- strptime(cd$germ_date,format="%d/%m/%Y")
dystrip <- strptime(cd$death_date,format="%d/%m/%Y")

gyday <- gystrip$yday
dyday <- dystrip$yday
gyear <- gystrip$year + 1900
dyear <- dystrip$year + 1900

iswinter_g <- gyear!=cd$year
iswinter_d <- dyear!=cd$year

cd$gseasday <- seasdaycalc(gyday,gyear,iswinter_g)
cd$dseasday <- NA
ddatepres <- (is.na(dyday)|is.na(dyear))==F
cd$dseasday[ddatepres] <- seasdaycalc(dyday[ddatepres],dyear[ddatepres],iswinter_d[ddatepres])
	# ignores leap years
	# if year matches year in cd, this is year FOLLOWING germination, so add 365
cd$gseasday[cd$gseasday>cd$dseasday] <- NA
	# Sets incorrect values to missing (Ursula looking into)

### SET DAYS
startday <- quantile(cd$gseasday,0.10,na.rm=T) 
endday <- quantile(cd$dseasday,0.90,na.rm=T)

### CALCULATE T TIMES

T1start <- median(cd$gseasday,na.rm=T)
T2start <- 31 + 14 # Jan + 1/2Feb
T3start <- median(cd$dseasday[cd$seeds>0],na.rm=T)

T1 <- (T2start - T1start)/365
T2 <- (T3start - T2start)/365
T3 <- 1-(T1+T2)
	# times in years

Tvalues <- data.frame(period=c("T1","T2","T3"),duration=c(T1,T2,T3))

write.csv(Tvalues,file=paste0("Tvalues_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)

# quantile(cd$gseasday,probs=c(0.1,0.9),na.rm=T)
	# germination dates mostly lay between -59 and 32

###################
### DATE FIGURE ###
###################

### MERGE INTO ONE VARIABLE
require(reshape2)
cdm <- subset(cd,select=c("year","plot","gseasday","dseasday"))
cdm <- melt(cdm, id.var=c("year","plot"))

### PLOT FIGURE

pdf("germ_death_seasdates.pdf",width=8,height=4)

layout(rbind(rep(1:2,c(1,4))))
par(bty="l")
par(mar=c(10,5,2,2),oma=c(0,1,2,0))
boxplot(value~variable,range=0,lty=1,las=2,data=cdm,border=c("red","blue"),
	names=c("Germination","Death / final survey"),ylab="Day (relative to 1st Jan)")
abline(h=startday,lty=2)
abline(h=endday,lty=2)
lettlab2(1,myline=0.5,cex=0.9)

par(mar=c(10,2,2,2))
boxplot(value~variable+year,
	# at=c(1:3, 4:6 + 1),
	range=0,
	lty=1,
	las=2,
	data=cdm,
	border=c("red","blue"),
	xaxt="n",
	xlab="Year"
	)
axis(side=1,
	at=seq(from=2,to=length(unique((cdm$year)))*2,by=2)-0.5, 
		# -0.5 puts in middle
	sort(unique(cdm$year)),
	lty="blank",
	las=2
	)

abline(h=startday,lty=2)
abline(h=endday,lty=2)
lettlab2(2,myline=0.5,cex=0.9)

legend("topright",
	legend=c("Germination","Death / final survey"),
	lty=c(1,1),
	col=c("red","blue"),
	bty="n")

dev.off()





