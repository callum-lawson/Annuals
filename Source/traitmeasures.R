# Setup -------------------------------------------------------------------

library(plyr)
library(reshape2)

### Load data

source("Source/figure_functions.R")
source("Source/prediction_functions.R")
source("Source/trait_functions.R")

msy <- read.csv("Output/msy_15Jan2016.csv",header=T)
csyp <- read.csv("Output/csyp_15Jan2016.csv",header=T)
cd <- read.csv("Output/cd_15Jan2016.csv",header=T)
cdpos <- read.csv("Output/cdpos_15Jan2016.csv",header=T)
Tvalues <- read.csv("Output/Tvalues_31Jul2015.csv",header=T)

# Define params -----------------------------------------------------------

prdat <- subset(csyp, ngerm>0 & is.na(nmeas)==F & is.na(germd)==F)
rsdat <- cdpos

### Data params
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
yrvals <- unique(msy$year)
nyear <- length(yrvals)

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])
# T in years

tau_d <- 100 # adjustment for density
tau_p <- 100 # adjustment for rainfall
tau_s <- 100 # adjustment for seed density

### Prediction params

nseq <- 10^3
niter <- 500 # nrow(rs$sig_y_r)

### Data used for model-fitting

msy$spyear <- with(msy, as.numeric(factor(species:as.factor(year))))
msy_merge <- subset(msy,select=c("species","year","spyear"))
prdat <- merge(prdat,msy_merge,by=c("species","year"),all.x=T)
rsdat <- merge(rsdat,msy_merge,by=c("species","year"),all.x=T)

prdat$species <- match(prdat$species,spvals)
prdat$year <- match(prdat$year,yrvals)
rsdat$species <- match(rsdat$species,spvals)
rsdat$year <- match(rsdat$year,yrvals)

prdat <- subset(prdat,select=c(species,year,prcp,germd))
rsdat <- subset(rsdat,select=c(species,year,prcp,germd))

expandscale <- 50 # for increasing range of z limits
minprcp <- min(prdat$prcp)*expandscale
maxprcp <- max(prdat$prcp)/expandscale

prdat$prcp[match(1:nspecies,prdat$species)] <- minprcp
prdat$prcp[match(1:nspecies,prdat$species)+1] <- maxprcp
prdat$prcp[(prdat$prcp %in% range(prdat$prcp))==F] <- NA
# rsdat$prcp[match(1:nspecies,rsdat$species)] <- minprcp
# rsdat$prcp[match(1:nspecies,rsdat$species)+1] <- maxprcp
# rsdat$prcp[(rsdat$prcp %in% range(rsdat$prcp))==F] <- NA
  # hack to make sure that each species has min and max prcp (rest of rows NA)

# using ngerm/area for each csyp plot - so taking mean of this

# PR

pr_prcp_lp <- r_lp(pr$beta_p,"pr","prcp",avmeth=median)
pr_dens_lp <- r_lp(pr$beta_p,"pr","dens")

# RS

rs_prcp_lp <- r_lp(rs$beta_r,"rs","prcp",avmeth=median)
rs_dens_lp <- r_lp(rs$beta_r,"rs","dens")

# PCR

pcr_prcp_lp <- list(
  xarr = pr_prcp_lp$xarr, 
  # pr = arbitrary
  pred_lp = log( plogis(pr_prcp_lp$pred_lp) * exp(rs_prcp_lp$pred_lp) )
  )
pcr_dens_lp <- list(
  xarr = pr_dens_lp$xarr,
  # pr = arbitrary
  pred_lp = log( plogis(pr_dens_lp$pred_lp) * exp(rs_dens_lp$pred_lp) )
  )

# Bet-hedging trait plots -------------------------------------------------

zseq <- makeseq(log(seq(minprcp,maxprcp,length.out=2)/tau_p)) 
  # was used in calculations but had not yet been saved as an object

rz <- pcr_prcp_lp$pred_lp
lz <- exp(rz)

zopt <- apply(rz,c(1,3),zopt_f,z=zseq)
zwid <- apply(rz,c(1,3),zwid_f,z=zseq)
rmax <- apply(rz,c(1,3),max,na.rm=T)
lsum <- apply(lz,c(1,3),sum)

rzrel <- aaply(exp(rz),c(1,3),function(x) x/sum(x))
rzrel <- log(aperm(rzrel,c(1,3,2)))

zopt_med <- apply(zopt,2,median,na.rm=T) # but check why getting NAs
zwid_med <- apply(zwid,2,median,na.rm=T)
rmax_med <- apply(rmax,2,median,na.rm=T)
lsum_med <- apply(lsum,2,median,na.rm=T)

rmed <- apply(rz,c(2,3),median)
rrelmed <- apply(rzrel,c(2,3),median)

### Plot traits

boxplot(zopt,type="l",horizontal=T,range=0)
boxplot(zwid,type="l",horizontal=T,range=0)

plot(zwid_med~zopt_med)
summary(lm(zwid_med~zopt_med))

plot(zwid_med~rmax_med)
summary(lm(zwid_med~rmax_med))

plot(lsum_med~rmax_med)

matplot(zseq[zseq>-1 & zseq<2],rmed[zseq>-1 & zseq<2,10:20],type="l",xlim=,lty=1)
matplot(zseq[zseq>-1 & zseq<2],exp(rmed[zseq>-1 & zseq<2,]),type="l")

matplot(zseq,rrelmed,type="l")
matplot(zseq,exp(rrelmed),type="l",xlim=c(-1,2),lty=1)

# not surprising that differ in seed numbers (e.g. if make smaller seeds)
# calculate relative reproduction instead?

# Mortality functions -----------------------------------------------------

# divide by K to re-scale curves?
# distinction between per-capita reproduction and total density
# DD in germination forces us to use simulation-derived traits? 

m0 <- exp(go$alpha_m[1:niter,])
m1 <- exp(go$beta_m[1:niter,])
lKn <- log(Kncalc(m0,m1,T3=T3))
lhn <- log(hncalc(m0,m1,T3=T3))

m0med <- apply(m0,2,median)
m1med <- apply(m1,2,median)
lKnmed <- apply(lKn,2,median)
lhnmed <- apply(lhn,2,median)

godmean <- with(go,godmean_f(alpha_G,beta_Gz)) # might be affected by density
godvar <- with(go,godvar_f(beta_Gz))
godmeanmed <-  apply(godmean,2,median)
godvarmed <-  apply(godvar,2,median)

exp(lhnmed[1])*Sncalc(exp(lhnmed[1]),m0=m0med[1],m1=m1med[1],T3=T3)
exp(lKnmed[1])/2
   # checking half-saturation density - looks correct

nd <- 1000
dseq <- exp(seq(-5,10,length.out=nd))
Snseq <- matrix(nr=nd,nc=nspecies)
for(i in 1:nspecies){
  Snseq[,i] <- qlogis(Sncalc(dseq,m0=m0med[i],m1=m1med[i],T3=T3))
  }
d2Sdn2 <- apply(Snseq,2,function(x) diff(diff(x)))
d2Sdn2_maxpos <- apply(Sncurve,2,function(x) which(x==min(x)))
ld2Sdn2max <- log(c(rep(NA,2),dseq)[maxcurve_dpos]) 
  # highly correlated with K and half-saturation

boxplot(lKn,type="l",horizontal=T,range=0)
boxplot(lhn,type="l",horizontal=T,range=0)
plot(lKnmed~lhnmed)
  # K and half-saturation basically equivalent
  # so large maximum seed density -> larger density increase needed to reach half of this

par(mfrow=c(1,2))
library(fields)
loglims <- 10^c(0,4)
i <- 1
curve(x*Sncalc(x,m0=m0med[i],m1=m1med[i],T3=T3),
  xlim=loglims,ylim=loglims,log="xy"
  )
# abline(v=exp(lhnmed[i]),lty=3)
for(i in 2:nspecies){
  curve(x*Sncalc(x,m0=m0med[i],m1=m1med[i],T3=T3),
      add=T,col=tim.colors(nspecies)[i]
      )
  }
# abline(h=Kn[i],col="red",lty=2)

library(fields)
loglims <- 10^c(0,4)
i <- 1
curve(qlogis(Sncalc(x,m0=m0med[i],m1=m1med[i],T3=T3)),
  xlim=loglims,log="x"
  )
# abline(v=exp(lhnmed[i]),lty=3)
for(i in 2:nspecies){
  curve(qlogis(Sncalc(x,m0=m0med[i],m1=m1med[i],T3=T3)),
    add=T,col=tim.colors(nspecies)[i]
    )
  }

# Save traits to RDS file -------------------------------------------------

medtraits <- data.frame(
  zopt=zopt_med,
  zwid=zwid_med,
  rmax=rmax_med,
  lsum=lsum_med,
  lKnmed=lKnmed,
  lhnmed=lhnmed,
  ld2Sdn2max=ld2Sdn2max,
  godmean=godmeanmed,
  godvar=godvarmed
  )
# because zwid is quantile, no difference if take absolute or 
# relative growth rates

saveRDS(medtraits,
  paste0("Output/medtraits_",format(Sys.Date(),"%d%b%Y"),".rds")
  )
