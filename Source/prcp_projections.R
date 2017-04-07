###################################
### Running mean of prcp values ###
###################################

library(reshape2)
library(tidyr)
library(plyr)
library(zoo)

# Load Data ---------------------------------------------------------------

allfiles <- list.files(paste0(getwd(),"/Data/Tucson_prcp_monthly/bcsd5"))
pfiles <- allfiles[grep(".csv",allfiles)]
pnames <- gsub(".csv","",pfiles)
pnames <- gsub("pr.","",pnames)
pmeta <- colsplit(pnames,pattern=".rcp",names=c("model","scenario"))

plist <- list()
np <- length(pnames)
for(i in 1:np){
  plist[[i]] <- read.csv(paste0(getwd(),"/Data/Tucson_prcp_monthly/bcsd5/",pfiles[i]),
                         header = FALSE,
                         col.names = c("year","mon","prcp","prcp_exl") 
                         )
  plist[[i]] <- subset(plist[[i]],select=c("year","mon","prcp"))
    # 2nd pr excluded because from further-away location
  plist[[i]]$model <- pmeta$model[i]
  plist[[i]]$scenario <- pmeta$scenario[i]
  }

pp <- rbind.fill(plist)
  # prcp = mean daily precipitation (mm/day)

# Calculate period --------------------------------------------------------

# startday <- 364-59 
#   yday has 364 when not leap year
#    = 3rd Nov
#   (gprcp starts 30 days before, so = 4th Oct)
# endday <- 112
#   = 23rd Apr
# earlyendday <- 32 
#   90% quantile for germination dates
#   = 1st Feb
# conclusion: 
#   growing season: Nov - Apr
#   germination season: Oct - Jan

startmon <- 11
endmon <- 4
earlystartmon <- 11 - 1 # 30 days = 1 mon
earlyendmon <- 1
  # these months in usual notation, rather than R notation (which starts at 0)
pdur_d <- 59 + 112
gdur_d <- (59+30) + 32
pdur_m <- 6
gdur_m <- 4

ppa <- ddply(pp, .(year,model,scenario), summarise,
  springprcp = sum(prcp[mon<=endmon]),
  wintprcp = sum(prcp[mon>=startmon]),
  espringprcp = sum(prcp[mon<=earlyendmon]),
  ewintprcp = sum(prcp[mon>=earlystartmon])
  )

ppa$prevyear <- ppa$year - 1
ppy <- ddply(ppa, .(model,scenario), summarise,
  year = year,
  model = model,
  scenario = scenario,
  seasprcp = (springprcp + wintprcp[match(prevyear,year)]) / pdur_m * pdur_d,
  germprcp = (espringprcp + ewintprcp[match(prevyear,year)]) / gdur_m * gdur_d
  )
ppy <- na.omit(ppy)
  # removes 1950 (no winter from 1949)

ppy$year <- as.numeric(ppy$year)
ppy$model <- as.factor(ppy$model)
ppy$scenario <- as.factor(ppy$scenario)

ppy$yearcat <- with(ppy,factor(
  ifelse(year<2015,"0",ifelse(year<2050,"50","100")),
  levels=c("0","50","100"))
  )

pps <- ddply(ppy, .(yearcat,model,scenario), summarise,
  pam = mean(log(seasprcp)),
  psd = sd(log(seasprcp)),
  gam = mean(log(germprcp)),
  gsd = sd(log(germprcp))
  )
  # p -> seasprcp, g -> germprcp

ppm = ddply(pps, .(model,scenario), summarise,
  yearcat=yearcat,
  mpam = pam/pam[yearcat=="0"],
  mpsd = psd/psd[yearcat=="0"],
  mgam = gam/gam[yearcat=="0"],
  mgsd = gsd/gsd[yearcat=="0"]
  )

myhist <- function(x){
  hist(x,breaks=100,main="")
  abline(v=1,col="red")
  }

ppm2 <- gather(ppm,measure,value,-c(model,scenario,yearcat))
ppma <- ddply(subset(ppm2,yearcat!="0"),.(measure,yearcat,scenario),summarise,median=median(value))

par(mfcol=c(2,4))
myhist(ppm$mpam[ppm$yearcat=="50"])
myhist(ppm$mpam[ppm$yearcat=="100"])
myhist(ppm$mpsd[ppm$yearcat=="50"])
myhist(ppm$mpsd[ppm$yearcat=="100"])
myhist(ppm$mgam[ppm$yearcat=="50"])
myhist(ppm$mgam[ppm$yearcat=="100"])
myhist(ppm$mgsd[ppm$yearcat=="50"])
myhist(ppm$mgsd[ppm$yearcat=="100"])

write.csv(ppma,file=paste0("Output/prcp_projection_summaries_",format(Sys.Date(),"%d%b%Y"),".csv"),
  row.names=F)

# Observed climate --------------------------------------------------------

ncy <- read.csv("Output/ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
	# removes first value (missing because no previous winter)
  # prcp values given in total mm (NOT 10ths of a mm) over whole period

pmroll <- rollmean(ncy$seasprcp,20)
pcroll <- rollapply(ncy$seasprcp,20,sd)
yroll <- round(rollmean(ncy$year,20),0)
plot(pmroll~yroll)
plot(pcroll~yroll)
