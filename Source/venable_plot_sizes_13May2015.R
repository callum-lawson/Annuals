#############################################
### Check and edit Venable plot area data ###
#############################################

# CHECK DATES ON SOURCE SCRIPTS

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
require(plyr)

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

#########################################
### CHECKS FOR VENABLE PLOT AREA DATA ###
#########################################

# pa <- read.csv("plotareas_bestguess.csv",header=T)
pa <- read.csv(paste0(maindir,"Data/Venable/","venablePlots2_09May2015.csv"),header=T)
	# Note that different working directory

excludeplots <- paste0(rep(paste0("C",c(1:3,5:6)),each=3),rep(LETTERS[1:3],times=5))
pa <- subset(pa,PlotID %in% excludeplots==F)
	# no longer getting rid of 2006 (or 2013/2014)
pa$PlotID <- as.factor(as.character(pa$PlotID))

###################
### EXPLORATION ###
###################

nlevels(pa$PlotID)
	# N plots = 79; more than 72 claimed

### PLOT X YEAR MATRICES
with(cd,table(plot,year))
with(pa,table(PlotID,year))

with(pa, table(year,Plot.,Habitat,Replicate))
with(subset(pa, Habitat=="open" & Replicate=="west"),table(PlotID,year))
only35 <- pa[grepl("3-", pa$PlotID)|grepl("5-", pa$PlotID),]
only35$PlotID <- as.factor(as.character(only35$PlotID))
with(only35,table(PlotID,year))
	# 3-open-west appears for 3 years in which 5-open-west disappears

notcurrent <- subset(pa,current==F)$PlotID
subset(pa,PlotID %in% notcurrent)
	# current missasigned? 

with(pa, table(PlotID,Plot.))
	# 6-open-east coded as Plot# 7,8,9 in three separate years
	# Plots 7 and 8 only appear once each
	# Plot 3 has east and west for open (3 years) but not for shrub
		# (but mismathes are all north-south)

subset(pa,PlotID=="6-open-east")
	# -> 6-open-east coded as Plot# 7,8,9 instead of 6 (typo)

subset(pa,Plot.==180 & Habitat=="open" & Replicate=="north")
with(subset(cd,plot=="180-open-north"),table(year))
	# 1990-1991 genuinely missing from 180-open-north even though 180-open-south present
	# (both datasets)

subset(pa,PlotID=="0/3-shrub-north")
	# Plot# = NA for 2001-2014; accidentally ommited?

with(subset(pa,grepl("3-", PlotID)),table(as.factor(as.character(PlotID)),year))
with(subset(cd,grepl("3-", plot)),table(as.factor(as.character(plot)),year))
subset(cd,plot=="3-open-east" & year==2012)

cdplots <- levels(cd$plot)
paplots <- levels(pa$PlotID)

subset(pa, PlotTypeAbbr=="gfp")
	# What is the difference between germination fraction plot and long term plot?

setdiff(cdplots,paplots) # length 8
setdiff(paplots,cdplots) # length 15 (probably not same plots)

miscd <- cd[cd$plot %in% paplots==F,]
miscd$plot <- as.factor(as.character(miscd$plot))
summary(miscd)
with(miscd, table(plot,year))
	# 30 only present in 1991-1992
	# 9 present for most years 1990-2008, except shrub absent until 1994
	# neither present past 2008
with(pa,tapply(area_m2,list(PlotID,year),mean))

### PLOT SIZES
with(pa,tapply(area_m2,list(Plot.,year,Habitat),mean))
subset(pa,Plot.==60 & year %in% 2005:2008)
	# In 2007, 60-open-north has area of 0

prescd <-  cd[cd$plot %in% paplots==T,]
prescd$plot <- as.factor(as.character(prescd$plot))

cdpytab <- with(prescd,table(plot,year)>0)
papytab <- with(pa,table(PlotID,year)>0)

pytabdiff <- cdpytab - papytab
pytabdiff>0
subset(pa,PlotID %in% c("3-open-east","3-shrub-north"))
	# 3-open-east and 3-shrub-north in 2012 are present in cd but not in pa

######################################################
### CORRECT VENABLE DATA AND SAVE NEW FILE VERSION ###
######################################################

# NB: Some changes made to CD in "dataprocess" file

# 2006 already removed from pa

### OLD ERRORS
# some no longer present; others must still be edited

# pa$area_m2[pa$PlotID=="60-open-north" & pa$year==2007] <- 0.1
# pa$Plot.[pa$PlotID=="0/3-shrub-north"] <- "0.3"
# pa$Plot.[pa$PlotID=="6-open-east" & pa$Plot. %in% 7:9] <- 6
# pa$current[pa$PlotID %in% c("3-open-east", "3-open-west")] <- F

### NEW ERRORS
pa$PlotID[pa$PlotID=="30-shurb-south"] <- "30-shrub-south"

names(pa)[names(pa)=="PlotID"] <- "plot"
names(pa)[names(pa)=="plot.type"] <- "type"
names(pa)[names(pa)=="PlotTypeAbbr"] <- "typeabbr"
names(pa)[names(pa)=="Plot."] <- "plotno"
names(pa)[names(pa)=="Habitat"] <- "habitat"
names(pa)[names(pa)=="Replicate"] <- "replicate"
names(pa)[names(pa)=="area_m2"] <- "area"

write.csv(pa,file=paste0("venablePlots_processed_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)




