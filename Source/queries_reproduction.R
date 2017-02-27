######################################################################
### Subsets and saves census data to provide examples of potential ###
### problems for Venable to investigate                            ###
######################################################################

# CHECK DATES ON SOURCE SCRIPTS

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))

########################
### LOAD CENSUS DATA ###
########################

cd <- read.csv("census_data_sharedspecies.csv",header=T) # cd = census data

names(cd)[names(cd)=="plot.habitat.replicate"] <- "plot"
names(cd)[names(cd)=="Year"] <- "year"

cd[cd %in% c("","-")] <- NA 	# assume that "" AND "-" = MISSING data
							# (applied to whole database)

cd$seeds <- suppressWarnings( as.numeric(as.character(cd$seeds)) )
	# removed warning about NAs

cd$reprod <- cd$seeds!=0 	# Did plant have >0 offspring?
						# !=0 because want to include -99 as reproduction

### CORRECT PLOT NAMES IN CD

cd <- cd[-(cd$year==2012 & cd$plot=="3-open-east"),]
cd$plot[cd$year==2012 & cd$plot=="3-shrub-north"] <- "3-shrub-south"

# remove these later

cd$plot <- as.factor(as.character(cd$plot))

### LOAD PLOT AREAS AND ADD TO CD

pa <- read.csv("venablePlots_02May2015.csv",header=T)

spa <- subset(pa,select=c("plot","year","area"))
	# add open/shrub, aspect if needed

cd <- merge(cd,spa,by=c("year","plot"),all.x=T)
	# not all.y because the following plots are missing from spa:

plots_open <- c("30-open-north","30-open-south","9-open-east","9-open-north","9-open-south","9-open-west")
plots_shrub <- c("30-shrub-north","30-shrub-south","9-shrub-east","9-shrub-north","9-shrub-south","9-shrub-west") 
	# these plots are missing from venablePlots

cd$area[cd$plot %in% plots_open] <- 0.1
cd$area[cd$plot %in% plots_shrub] <- 0.05
	# correct plot areas in CD
cd <- cd[is.na(cd$area)==F,]	
	# Delete 2012 3-open-east (discovered as incorrect plot)

cd$seeds[cd$seeds>5000] <- NA
	# Such high values clearly wrong

###########################
### LOAD SEED BANK DATA ###
###########################

sb <- read.csv("Seed_Bank_sharedspecies.csv",header=T)

names(sb)[names(sb)=="viable.seeds"] <- "nlive"
names(sb)[names(sb)=="non.viable.seeds"] <- "ndead"

sb[sb==-99] <- NA 	# assume that -99 is MISSING data

############################
### AGGREGATION - CENSUS ###
############################

# NB: in all cases, average over lower levels first
# Census: plots -> years
# Seeds: replicates -> plots -> years

### AGGREGATE (cd -> csyp)

csyp <- ddply(cd, .(species,year,plot), summarize,
	ngerm = length(reprod), 				
		# length of dataframe
	nfert = ifelse(F %in% is.na(reprod), sum(reprod==T,na.rm=T), NA),
	nbarren = ifelse(F %in% is.na(reprod), sum(reprod==F,na.rm=T), NA),
	nmeas = nfert + nbarren,
	pr = ifelse(F %in% is.na(reprod), mean(reprod,na.rm=T), NA)	,
		# pr = NA when all plants had reprod=NA
	rseed = mean(seeds[seeds>0]), 		
		# n seeds given that reproduced; NaN when rseed = 0 or -99
	area = mean(area,na.rm=T)
	) 
	# csyp = census species year plot
	# seedlings germinated but reproduction unknown -> accounted for by ommiting NAs
	# nmeas = number of seedlings in which reproduction was measured
	# (differs from ngerm if seeds germinated but reproduction unknown)

### ADD NGERM=0

require(data.table)

pa_dt <- as.data.table(subset(pa,select=c("year","plot","area")))
	# area just serves as way to test whether combo is present
setkey(pa_dt,year,plot)
pa_dt <- pa_dt[CJ(unique(year),unique(plot))]
pa_dt <- as.data.frame(pa_dt)

names(pa_dt)[names(pa_dt)=="area"] <- "plotpres"
pa_dt$plotpres <- ifelse(is.na(pa_dt$plotpres),F,T)

csyp_dt <- as.data.table(csyp)
	# csyp data table
setkey(csyp_dt,species,year,plot)
csyp_dt <- csyp_dt[CJ(unique(species),unique(year),unique(plot))]
	# ignores 2006 - could set year to sequence to include

csyp_pa_dt <- as.data.frame(merge(csyp_dt,pa_dt,by=c("year","plot"),all=T))
	# adds "plotpres" to each year-plot combo in csyp_dt
# nrow(csyp_dt) - nrow(csyp_pa_dt)
	# = 0 -> all year-plot combos present in both datasets

plots_30 <- c("30-open-north","30-open-south","30-shrub-north","30-shrub-south")
plots_9 <- c("9-open-east","9-open-north","9-open-south","9-open-west",
	"9-shrub-east","9-shrub-north","9-shrub-south","9-shrub-west")  
csyp_pa_dt$plotpres[csyp_pa_dt$plot %in% plots_30 & csyp_pa_dt$year %in% 1991:1992] <- T
csyp_pa_dt$plotpres[csyp_pa_dt$plot %in% plots_30 & csyp_pa_dt$year %in% 1991:1992==F] <- F
csyp_pa_dt$plotpres[csyp_pa_dt$plot %in% plots_9 & csyp_pa_dt$year %in% 1990:2008] <- T
csyp_pa_dt$plotpres[csyp_pa_dt$plot %in% plots_9 & csyp_pa_dt$year %in% 1990:2008==F] <- F

csyp_pa_dt$ngerm[is.na(csyp_pa_dt$ngerm) & csyp_pa_dt$plotpres==T] <- 0
csyp <- as.data.frame(csyp_pa_dt[is.na(csyp_pa_dt$ngerm)==F,])
	# table(csyp$ngerm==0)/nrow(csyp)
	# ~80% of plots have no germinants!

### NEW VARIABLES (csyp)

csyp$pcr = with(csyp, ifelse(pr==0,0,ifelse(is.nan(rseed), NA, pr*rseed)))
	# per-capita reproduction (mean seeds per germinant)
	# = 0 if pr=0 
	# = NA if(pr>0 AND rseed=NaN)
		# NaN problem = don't know whether 0 reproducers died or had 0 offspring

csyp$newseed = with(csyp, ngerm*pcr)
	# newseed = estimated seeds per plot 
	# = number of individuals * mean seeds per individual

csyp$germd <- with(csyp, ifelse(ngerm==0,0,ngerm/area))	
	# germinant density - when zero, area missing, so would get NA
csyp$nsd <- with(csyp, newseed/area)	# new seed density 

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
	pr = mean(pr,na.rm=T),
	rseed = mean(rseed,na.rm=T),
	pcr = mean(pcr,na.rm=T),
	newseed = mean(newseed,na.rm=T),
	germd = mean(germd,na.rm=T),
	nsd = mean(nsd,na.rm=T),
	prevnsd = mean(prevnsd,na.rm=T)
	)

csy[csy=="NaN"] <- NA # get rid of NaNs arising from averaging rows without germinants

###############################
### AGGREGATION - SEED BANK ###
###############################

### AGGGREGATE (sb -> ssyp)

ssyp <- ddply(sb, .(species,year,plot), summarize,
	nplot = length(na.omit(nlive)),
	nlive = mean(nlive,na.rm=T), 
	ndead = mean(ndead,na.rm=T)
	) 
	# ssyp = seed species year plot

### NEW VARIABLES (ssyp)

corearea <- pi*(0.054/2)^2		# diameter is 5.4cm

ssyp$olsd <- with(ssyp, nlive/corearea)		# old live seed dens
ssyp$odsd <- with(ssyp, ndead/corearea)		# old dead seed dens

### PREVIOUS YEARS

ssyp$prevyear <- ssyp$year - 1
ssp <- ddply(ssyp, .(species,plot), summarize,
	species = species,
	plot = plot,
	year = year,
	prevolsd = olsd[match(prevyear,year)],
	prevodsd = odsd[match(prevyear,year)]
	)

ssyp <- merge(ssyp,ssp,by=c("species","year","plot"))

### AGGREGATE (ssy -> ssyp)

ssy <- ddply(ssyp, .(species,year), summarize,
	olsd = mean(olsd,na.rm=T), 	
	odsd = mean(odsd,na.rm=T),
	prevolsd = mean(prevolsd,na.rm=T), 	
	prevodsd = mean(prevodsd,na.rm=T)			
	)
	# ssy = seed species year

###################################
### MERGED CENSUS AND SEED BANK ###
###################################

msy <- merge(csy,ssy,by=c("species","year"),all.x=T,all.y=T)
	# merged (census and seed) species year
	# includes all data, even when year is missing from either dataset
	# (e.g. prior to 1990/1991, when seed surveys started)

msy$csd <- with(msy, nsd + olsd) # combined (live) seed density
msy$prevcsd <- with(msy, prevnsd + prevolsd)

msy$csdr <- with(msy, log(csd/prevcsd)) 
	# r = intrinsic growth rate

####################
### APPEND MEANS ###
####################

### msy -> csyp AND ssyp

msy_app <- with(msy, data.frame(
	species=species,
	year=year,
	olsdbar=olsd,
	nsdbar=nsd,
	germdbar=germd,
	prevnsdbar=prevnsd,
	prevolsdbar=prevolsd
	))

csyp <- merge(csyp,msy_app,by=c("species","year"))
ssyp <- merge(ssyp,msy_app,by=c("species","year"))

### csyp (ngerm) -> cd

cysp_app <- subset(csyp,select=c("species","year","plot","ngerm"))
cd <- merge(cd,cysp_app,by=c("species","year","plot"))
cd$germd <- with(cd, ngerm/area) 

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

##################
### CLEAN DATA ###
##################

rm("csp","cysp_app","csyp_dt","csyp_pa_dt","msy_app","pa","pa_dt", "plots_30","plots_9","plots_open","plots_shrub","ssp","spa")
# day variables removed

gc()

