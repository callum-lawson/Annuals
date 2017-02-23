###########################################################################
### Sets up seedling and seed databases for later analysis and plotting ###
###########################################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
source("venable_functions_22Apr2015.R")

library(plyr)
library(data.table)

########################
### LOAD CENSUS DATA ###
########################

cd <- read.csv("census_data_sharedspecies_13May2015.csv",header=T) # cd = census data

names(cd)[names(cd)=="plot.habitat.replicate"] <- "plot"
names(cd)[names(cd)=="Year"] <- "year"

cd[cd %in% c("","-")] <- NA 	# assume that "" AND "-" = MISSING data
							# (applied to whole database)

cd$seeds <- suppressWarnings( as.numeric(as.character(cd$seeds)) )
	# removed warning about NAs

cd$reprod <- cd$seeds!=0 	# Did plant have >0 offspring?
						# !=0 because want to include -99 as reproduction

cd <- subset(cd,species!="fica")
cd$species <- as.factor(as.character(cd$species))
	# Removes Fica

###########################
### LOAD SEED BANK DATA ###
###########################

sb <- read.csv("Seed_Bank_sharedspecies_13May2015.csv",header=T)

names(sb)[names(sb)=="viable.seeds"] <- "nlive"
names(sb)[names(sb)=="non.viable.seeds"] <- "ndead"

sb[sb==-99] <- NA 	# assume that -99 is MISSING data

sb <- subset(sb,species!="fica")
sb$species <- as.factor(as.character(sb$species))
	# Removes Fica

### INDIVIDUAL CORRECTIONS

mistakes <- with(sb, which(
	(year==2014 & plot=="E" & habitat=="open" & replicate==5 & species=="scba" & nlive==0 & ndead==1)
	| (year==2014 & plot=="I" & habitat=="shrub" & replicate==4 & species=="pehe" & nlive==0 & ndead==1)
	| (year==2014 & plot=="L" & habitat=="shrub" & replicate==3 & species=="pehe" & nlive==0 & ndead==0)
	))

sb <- sb[-mistakes,]
sb$species[sb$year==2014 & sb$plot=="L" & sb$habitat=="shrub" & sb$replicate==3 & sb$species=="pere" & sb$nlive==2 & sb$ndead==1] <- "pehe"

### CHECK FOR DUPLICATES
sb_cat <- sb[,names(sb) %in% c("nlive","ndead")==F]
dupfor <- duplicated(sb_cat,fromLast=F)
dupback <- duplicated(sb_cat,fromLast=T)
	# one plot duplicated in 2013 
# sb[dupfor|dupback,]

### DUPLICATE CORRECTIONS
sb$plot[dupback & sb$plot=="G"] <- "F"
sb <- sb[-which(dupfor & sb$plot!="F"),]
	# deleting/moving first records (dupfor) is arbitrary

### PAD OUT WITH NAS
sb_dt <- as.data.table(sb)
	# sb data table

setkey(sb_dt,species,year,plot,habitat,replicate)
sb_dt <- sb_dt[CJ(
	unique(species),
	unique(year),
	unique(plot),
	unique(habitat),
	unique(replicate)
	),
	allow.cartesian=TRUE]

sb <- as.data.frame(sb_dt)

yearvals <- unique(sb$year)
area_byyear <- vector(length=length(yearvals))
area_byyear[yearvals %in% c(1990:1992,2001:2005)] <- 0.002290221
area_byyear[yearvals %in% 1993:2000] <- 0.00194
area_byyear[yearvals %in% 2006:2014] <- 0.002156515
	# From "plot and core areas over time_01Oct2015"

# sb$area <- pi*(0.054/2)^2	# diameter is 5.4cm
sb$area <- area_byyear[as.numeric(as.factor(sb$year))]
sb$nlive[with(sb, is.na(nlive) & year %in% 2013:2014 & species %in% c("erdi","lela","siir"))] <- 0

# sb$area[is.na(sb$nlive)==T] <- NA
	# ignoring areas where ndead measured but nlive not

############################
### AGGREGATION - CENSUS ###
############################

# NB: in all cases, average over lower levels first
# Census: plots -> years
# Seeds: replicates -> plots -> years

mysum <- function(x) ifelse(F %in% is.na(x),sum(x,na.rm=T),NA)
mymean <- function(x) ifelse(F %in% is.na(x),mean(x,na.rm=T),NA)

### AGGREGATE (cd -> csyp)

csyp <- ddply(cd, .(species,year,plot), summarize,
	ngerm = length(reprod), 				
		# length of dataframe
	nfert = mysum(reprod==T),
	nbarren = mysum(reprod==F),
	nmeas = nfert + nbarren,
	pr = mymean(reprod),
		# pr = NA when all plants had reprod=NA
	rseed = mean(seeds[seeds>0])
		# n seeds given that reproduced; NaN when rseed = 0 or -99
	) 
	# csyp = census species year plot
	# seedlings germinated but reproduction unknown -> accounted for by ommiting NAs
	# nmeas = number of seedlings in which reproduction was measured
	# (differs from ngerm if seeds germinated but reproduction unknown)

### ADD NGERM=0 

pa <- read.csv("venablePlots_processed_13May2015.csv",header=T)

pa_dt <- as.data.table(subset(pa,select=c("year","plot","area")))
	# add open/shrub, aspect if needed
	# area also serves as way to test whether combo is present
	# INCLUDES 2006

setkey(pa_dt,year,plot)
pa_dt <- pa_dt[CJ(unique(year),unique(plot))]
pa_dt <- as.data.frame(pa_dt)

pa_dt$plotpres <- ifelse(is.na(pa_dt$area),F,T)

csyp_dt <- as.data.table(csyp)
	# csyp data table
setkey(csyp_dt,species,year,plot)
csyp_dt <- csyp_dt[CJ(unique(species),min(year,na.rm=T):max(year,na.rm=T),unique(plot))]
	# includes 2006 by making year sequence

csyp_pa_dt <- as.data.frame(merge(csyp_dt,pa_dt,by=c("year","plot"),all=T))
	# adds "plotpres" to each year-plot combo in csyp_dt
# nrow(csyp_dt) - nrow(csyp_pa_dt)
	# = 0 -> all year-plot combos present in both datasets

csyp_pa_dt$ngerm[is.na(csyp_pa_dt$ngerm) & csyp_pa_dt$plotpres==T] <- 0
csyp <- as.data.frame(csyp_pa_dt)
	# Formerly: as.data.frame(csyp_pa_dt[is.na(csyp_pa_dt$ngerm)==F,])
	
	# The 2006 removal used to be here (see dataprocess from 01 Aug 2015)

	# table(csyp$ngerm==0)/nrow(csyp)
	# ~80% of plots have no germinants!

### NEW VARIABLES (csyp)

csyp$pcr = with(csyp, ifelse(pr==0,0,ifelse(is.nan(rseed), NA, pr*rseed)))
	# per-capita reproduction (mean seeds per germinant)
	# = 0 if pr=0 (pr=0 -> rseed=NA -> need to fill in 0 manually)
	# = NA if(pr>0 AND rseed=NaN)
	# NaN problem = don't know whether 0 reproducers died or had 0 offspring

csyp$newseed = with(csyp, ifelse(ngerm==0,0,ngerm*pcr))
	# newseed = estimated seeds per plot 
	# = number of individuals * mean seeds per individual
	# (ngerm=0 -> pcr=NA -> need to fill in 0 manually)

	# reproduction measurements missing from 10% of plots with germinants
	# that reproduced:
		# with(subset(csyp,pr>0),table(is.na(rseed)))
	# but doesn't seem to be biased:
		# with(subset(csyp,pr>0 & is.na(rseed)),mean(pr))
		# with(subset(csyp,pr>0 & is.na(rseed)==F),mean(pr))

csyp$germd <- with(csyp, ngerm/area)	
	# germinant density - when zero, area missing, so would get NA
csyp$nsd <- with(csyp, newseed/area)	
	# new seed density 

### PREVIOUS YEARS

csyp$prevyear <- csyp$year - 1
csp <- ddply(csyp, .(species,plot), summarize,
	species = species,
	year = year,
	plot = plot,
	prevnsd = nsd[match(prevyear,year)]
	)

csyp <- merge(csyp,csp,by=c("species","year","plot"))

### AGGREGATE (csyp -> csy)

csy <- ddply(csyp, .(species,year), summarize,	
	totgerm = mysum(ngerm),
	totfert = mysum(nfert),
	totbarren = mysum(nbarren),
	area_g = mysum(area[is.na(ngerm)==F]), 
		# area of germinant counts
	area_f = mysum(area[is.na(nfert)==F & is.na(nbarren)==F]), 
		# area of fecundity counts
	area_n = mysum(area[(is.na(nfert)==F & is.na(nbarren)==F) | ngerm==0]), 
		# area of seed counts (excludes plots missing fecundity)
	nplot_n = length(area[(is.na(nfert)==F & is.na(nbarren)==F) | ngerm==0]), 

	prbar = mymean(pr),
	rseedbar = mymean(rseed),
	pcrbar = mymean(pcr),
	totnewseed = mysum(newseed),
	newseedhat = totnewseed/area_g,
	newseedbar = mymean(newseed),
	germdhat = totgerm/area_g,
	germdbar = mymean(germd),
	nsdbar = mymean(nsd),
	prevnsdbar = mymean(prevnsd)
	)

# areas for germinants and nfert/nbarren very different because 
# most germinant plots don't have seedlings (so no reproduction measured)

# with(csy,(ngerm / area_g) - germd)
	# sum(ngerm) / sum(area) different to mean(germd) because latter
	# weights different plots equally (regardless of area)

csy[csy=="NaN"] <- NA # get rid of NaNs arising from averaging rows without germinants

###############################
### AGGREGATION - SEED BANK ###
###############################

# hat -> sum then mean
# bar -> mean then mean

### AGGGREGATE (sb -> ssyp)

ssyp <- ddply(sb, .(species,year,plot), summarize,
	nplot = length(na.omit(nlive)),
	totlive = mysum(nlive),
	area_o = mysum(area),
	olsdbar = mymean(nlive/area), 
	odsdbar = mymean(ndead/area)
	) 
	# !!! Could change to olsdhat using totlive/area_o
	# ssyp = seed species year plot

### PREVIOUS YEARS

ssyp$prevyear <- ssyp$year - 1
ssp <- ddply(ssyp, .(species,plot), summarize,
	species = species,
	plot = plot,
	year = year,
	prevolsdbar = olsdbar[match(prevyear,year)],
	prevodsdbar = odsdbar[match(prevyear,year)]
	)

ssyp <- merge(ssyp,ssp,by=c("species","year","plot"))

### AGGREGATE (ssyp -> ssy)

ssy <- ddply(ssyp, .(species,year), summarize,
	area_o = mysum(area_o[is.na(totlive)==F]),
	totlive = mysum(totlive),
	olsdbar = mymean(olsdbar), 	
	odsdbar = mymean(odsdbar),
	prevolsdbar = mymean(prevolsdbar), 	
	prevodsdbar = mymean(prevodsdbar)			
	)
	# ssy = seed species year

###################################
### MERGED CENSUS AND SEED BANK ###
###################################

msy <- merge(csy,ssy,by=c("species","year"),all.x=T,all.y=T)
	# merged (census and seed) species year
	# includes all data, even when year is missing from either dataset
	# (e.g. prior to 1990/1991, when seed surveys started)

msy$csdbar <- with(msy, nsdbar + olsdbar) # combined (live) seed density
msy$prevcsdbar <- with(msy, prevnsdbar + prevolsdbar)
msy$totmeas <- with(msy, totfert + totbarren)

####################
### APPEND MEANS ###
####################

### msy -> csyp AND ssyp

msy_app_c <- subset(msy,select=c("species","year","olsdbar","nsdbar",
	"germdhat","germdbar","prevnsdbar","prevolsdbar")	
	)
msy_app_s <- subset(msy,select=c("species","year","nsdbar",
	"germdhat","germdbar")	
	)

csyp <- merge(csyp,msy_app_c,by=c("species","year"))
ssyp <- merge(ssyp,msy_app_s,by=c("species","year"))
	# olsdbar here refers to value for that specific plot 

### csyp (ngerm) -> cd

cysp_app <- subset(csyp,select=c("species","year","plot","ngerm","area","germd"))
cd <- merge(cd,cysp_app,by=c("species","year","plot"))

### msy (ngermhat) -> cd

msy_app <- subset(msy,select=c("species","year","germdhat"))
cd <- merge(cd,msy_app,by=c("species","year"))

####################
### ADD CLIMATE  ###
####################

nc <- read.csv("aggregate_rainfall_03May2015.csv",header=T)
	# nc = national climate

nc$date <- strptime(nc$date,format="%Y-%m-%d")
	# leave as late as possible
nc$year <- nc$date$year + 1900 # actual year, not survey year
nc$month <- nc$date$mday
nc$yday <- nc$date$yday # 01 Jan = yday 0

### CHOOSE WINDOW
# Based on seasonwindow script

# Choose start and end days 
startday <- 364-59 # yday has 364 when not leap year
endday <- 112
earlyendday <- 32 # 90% quantile for germination dates
 	
ncy <- ddply(nc, .(year), summarize,
	springprcp = sum(prcp[yday<=endday], na.rm=T), 	
	wintprcp = sum(prcp[yday>=startday], na.rm=T),
	espringprcp = sum(prcp[yday<=earlyendday], na.rm=T),
	ewintprcp = sum(prcp[yday>=(startday-30)], na.rm=T)	# -30 because month before first germ 			
	) 
	# IGNORING NAS - only 8 in whole nc dataset
	# national climate aggregated 

prevyear <- ncy$year - 1
ncy$prevwintprcp <- with(ncy, wintprcp[match(prevyear,year)])
ncy$eprevwintprcp <- with(ncy, ewintprcp[match(prevyear,year)])

ncy$seasprcp <- with(ncy, prevwintprcp + springprcp)/10 
ncy$germprcp <- with(ncy, eprevwintprcp + espringprcp)/10
	# accounting for different years
	# /10 = convert to mm (original values in 10ths of a mm)

### APPEND TO DATA
 
csyp$prcp <- ncy$seasprcp[match(csyp$year,ncy$year)] 
cd$prcp <- ncy$seasprcp[match(cd$year,ncy$year)]
msy$prcp <- ncy$seasprcp[match(msy$year,ncy$year)]

csyp$gprcp <- ncy$germprcp[match(csyp$year,ncy$year)] 
cd$gprcp <- ncy$germprcp[match(cd$year,ncy$year)]
msy$gprcp <- ncy$germprcp[match(msy$year,ncy$year)]

###############################
### ADD INTEGER SEEDS TO CD ###
###############################

cd$seedsint <- with(cd, as.integer(round(seeds,0)))
	# Reproduction not always whole number
	# -> in some years, seeds have not been counted independently

	# which(cd$seeds>0 & cd$seeds<1) # empty 
		# -> nothing rounded down to 0

cdpos <- subset(cd,seedsint>0)

##########################
### SUMMARISE COVERAGE ###
##########################

#with(cd,table(species,year))
#with(csyp,table(species,year))
#with(csy,table(species,year))

#with(sb,table(species,year))
#with(ssyp,table(species,year))
#with(ssy,table(species,year))

#with(msy,table(species,year))
	# fica absent from both datasets in 2006, so excluded

	# Could use xtabs

##################
### WRITE DATA ###
##################

write.csv(cd,file=paste0("cd_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(cdpos,file=paste0("cdpos_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(csy,file=paste0("csy_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(csyp,file=paste0("csyp_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(msy,file=paste0("msy_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(sb,file=paste0("sb_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(ssy,file=paste0("ssy_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(ssyp,file=paste0("ssyp_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
write.csv(ncy,file=paste0("ncy_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)

#######################
### CLEAN WORKSPACE ###
#######################

rm("csp","cysp_app","csyp_dt","csyp_pa_dt","dupback","dupfor","mistakes",
	"msy_app","pa","pa_dt","sb_cat","sb_dt","ssp")
# day variables removed

gc()