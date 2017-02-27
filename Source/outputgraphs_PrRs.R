################################################
### Output graphs for Stan models, including ###
### go, reproduction, and population growth  ###
################################################

library(plyr)
library(reshape2)
library(truncdist)

### LOAD DATA 

source("Source/figure_functions.R")
source("Source/prediction_functions.R")

msy <- read.csv("Output/msy_15Jan2016.csv",header=T)
csyp <- read.csv("Output/csyp_15Jan2016.csv",header=T)
cd <- read.csv("Output/cd_15Jan2016.csv",header=T)
cdpos <- read.csv("Output/cdpos_15Jan2016.csv",header=T)
Tvalues <- read.csv("Output/Tvalues_31Jul2015.csv",header=T)

#####################
### DEFINE PARAMS ###
#####################

prdat <- subset(csyp, ngerm>0 & is.na(nmeas)==F & is.na(germd)==F)
rsdat <- cdpos

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
yrvals <- unique(msy$year)
nyear <- length(yrvals)

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])
	# T in years

	# CHANGE THESE

tau_d <- 100 # adjustment for density
tau_p <- 100 # adjustment for rainfall

### MCMC AND PREDICTION PARAMS

nseq <- 100
# niter <- max(length(goparams$lp__),length(prparams$lp__),length(rsparams$lp__))
	# go is longer than reproduction
niter <- nrow(rs$sig_y_r)
myprobs <- c(0.025,0.5,0.975)
nquant <- length(myprobs)

### DATA USED FOR MODEL-FITTING

msy$spyear <- with(msy, as.numeric(factor(species:as.factor(year))))
msy_merge <- subset(msy,select=c("species","year","spyear"))
prdat <- merge(prdat,msy_merge,by=c("species","year"),all.x=T)
rsdat <- merge(rsdat,msy_merge,by=c("species","year"),all.x=T)

prdat$spsite <- with(prdat, as.numeric(factor(species:plot)))
prdat$species <- match(prdat$species,spvals)
prdat$year <- match(prdat$year,yrvals)

rsdat$spsite <- with(rsdat, as.numeric(factor(species:plot)))
rsdat$species <- match(rsdat$species,spvals)
rsdat$year <- match(rsdat$year,yrvals)

# PR

# smalliter <- 200 # smaller number of samples for quicker processing during code development
pr_prcp_lp <- r_lp(pr$beta_p,"pr","prcp") # maxiter=smalliter
pr_dens_lp <- r_lp(pr$beta_p,"pr","dens") # maxiter=smalliter

pr_prcp_quant <- plogis(r_quant(pr_prcp_lp))
pr_dens_quant <- plogis(r_quant(pr_dens_lp))

# RS

rs_prcp_lp <- r_lp(rs$beta_r,"rs","prcp") # maxiter=smalliter
rs_dens_lp <- r_lp(rs$beta_r,"rs","dens") # maxiter=smalliter

rs_prcp_quant <- exp(r_quant(rs_prcp_lp))
rs_dens_quant <- exp(r_quant(rs_dens_lp))

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

### non-linear averaging thing going on?

pcr_prcp_quant <- exp(r_quant(pcr_prcp_lp))
pcr_dens_quant <- exp(r_quant(pcr_dens_lp))

# remove(pr_prcp_lp,pr_dens_lp,pr_prcp_quant,pr_dens_quant,
#	rs_prcp_lp,rs_dens_lp,rs_prcp_quant,rs_dens_quant,
#	pcr_prcp_lp,pcr_dens_lp,pcr_prcp_quant,pcr_dens_quant); gc()

# pcr_cal_quant

###########################
### CALIBRATION FIGURES ###
###########################

prdat$prpred <- plogis(pr$pi)
rsdat$seedspred <- rs$mu_r

# mu_r calculated in Stan model-fitting R script
# (already back-transformed, truncation already accounted for)

prdat$prpos <- 1:nrow(prdat)
rsdat$rspos <- 1:nrow(rsdat)
prrsdat <- merge(prdat,rsdat,by=c("species","year","plot"))
  # multiple other columns with same names but different entries

prrsdat$pcrpred <- with(prrsdat,prpred*seedspred)

rsdat0 <- ddply(rsdat,.(species),summarise,
	year = year,
	seedspred = seedspred,
	seedsint0 = ifelse(seedsint==0,min(seedsint[seedsint!=0]),seedsint)
	)
prrsdat0 <- ddply(prrsdat,.(species),summarise,
	year = year,
	pcrpred = pcrpred,
	pcr0 = ifelse(pcr==0,min(seedsint[pcr!=0]),pcr)
	)

# Aggregate, then plot...?

pdf(paste0("Plots/Pr_Rs_calibration_nozeroes_",format(Sys.Date(),"%d%b%Y"),".pdf"),
		width=9,height=3)
par(mfrow=c(1,3),mar=c(5,5,2,2),las=1,bty="l")

calplotall(
	xdat=prdat$prpred,
	ydat=prdat$pr,
	xlab="Predicted Pr(Y(t)>1)",
	ylab="Observed Pr(Y(t)>1)"
	)
calplotall(
	xdat=log(rsdat0$seedspred),
	ydat=log(rsdat0$seedsint0),
	xlab="Predicted ln(Y(t)|Y(t)>1)",
	ylab="Observed ln(Y(t)|Y(t)>1)"
	)
calplotall(
	xdat=log(prrsdat0$pcrpred),
	ydat=log(prrsdat0$pcr0),
	xlab="Predicted ln(Y(t))",
	ylab="Observed ln(Y(t))"
	)

	# Does truncated model do a better job of this?

dev.off()

########################
### SIMULATED COUNTS ###
########################

# csyp_num <- csyp # numeric so that can match to prrsdat
# csyp_num$species <- match(csyp_num$species,spvals)
# csyp_num$year <- match(csyp_num$year,yrvals)
# csyp_num$int <- rep(1,nrow(csyp_num))
# csyp_num$tprcp <- log(csyp_num$prcp/tau_p)
# csyp_num$tprcp2 <- csyp_num$tprcp^2
# csyp_num$tgermd <- log(csyp_num$germd/tau_p)
# 
# idvars <- c("species","year","plot")
# measvars <- c("int","tprcp","tprcp2","germd")
# 
# xprep <- subset(csyp_num,
#   select=c(idvars,measvars)
#   )
# xmelt <- melt(xprep,measure=measvars)
# xarr <- acast(xmelt,... ~ variable ~ species)
# 
# beta_p_med <- array(dim=c(1,dim(pr$beta_p)[-1]))
# beta_p_med[] <- apply(pr$beta_p,c(2,3),median)
# 
# beta_r_med <- array(dim=c(1,dim(rs$beta_r)[-1]))
# beta_r_med[] <- apply(rs$beta_r,c(2,3),median)
# 
# pr_d_lp <- r_lp(beta_p_med,"pr","data",xarr=xarr,maxiter=1)$pred_lp
# rs_d_lp <- r_lp(beta_r_med,"rs","data",xarr=xarr,maxiter=1)$pred_lp
# pcr_d_lp <- plogis(pr_d_lp) * exp(rs_d_lp)
# 
# csyp_num$pcr_exp <- as.vector(aperm(pcr_d_lp,c(1,3,2)))
# 
# head(subset(csyp_num,
#   germd>0 & !is.na(germd),
#   select=c(measvars[-c(1,3)],"pcr_exp")
#   ),
#   100)
# 
# nprrs <- nrow(prrsdat)
# 
# sig_o_p_med <- median(pr$sig_o_p)
# phi_r_med <- median(rs$phi_r)
# 
# prnoise <- rnorm(nprrs,mean=0,sd=sig_o_p_med)
# 
# prnoise <- rnorm(nprrs,mean=0,sd=sig_o_p_med)
# prrsdat$prsimprob <- plogis(qlogis(prrsdat$prpred) + prnoise)
# prrsdat$prsim <- rbinom(nprrs, size=1, prob=prrsdat$prsimprob)
# prrsdat$seedssim <- rnbinom(nprrs,mu=prrsdat$seedspred,size=phi_r_med)
# # truncation already accounted for (see above)
# prrsdat$pcrsim <- with(prrsdat,prsim*seedssim)
# 
# prrssum <- ddply(prrsdat,.(species,year,plot),summarise,
#   area_r = sum(area.x), # equivalent to area.y
#   newseed_sim=sum(pcrsim)
# )
# 
# csyp_sim <- merge(csyp_num,prrssum,all.x=T)
# 
# 
# write.csv(csyp_sim,
#           file=paste0("csyp_sim",format(Sys.Date(),"%d%b%Y"),".csv"),
#           row.names=F
#           )

############################
### RELATIONSHIP FIGURES ###
############################

### PREDICTIONS AND DATA

csyp_nz <- ddply(subset(csyp,year!="2006"),.(species),summarise,
	year=year,
	plot=plot,
	prcp=prcp,
	germd=germd,
	pr=pr,
	rsplus=ifelse(rseed==0,min(rs[rseed>0],na.rm=T),rseed),
	pcrplus=ifelse(pcr==0,min(pcr[pcr>0],na.rm=T),pcr)
	)
	# csyp no zeroes
	# no 2006 because no germination in that year (so no reproduction measures)

### PRCP 

# PR

preddatplot(
	plotname="pr_prcp",
	xdat=myjitter(log(csyp_nz$prcp)),
	ydat=myjitter(csyp_nz$pr),	
	sind=csyp_nz$species,
	xpred=pr_prcp_lp$xarr[,2,] + log(tau_p), # 2 = prcp
	ypred=pr_prcp_quant,
	xptype="matrix",
	pcol=myblue,
	xname=expression(z(t)),
	yname=expression(Pr(Y(t)>1)),
	ylim=c(0,1),
	xlim=c(3.3,5.8)
	)

# RS

preddatplot(
	plotname="rs_prcp",
	xdat=myjitter(log(csyp_nz$prcp)),
	ydat=myjitter(log(csyp_nz$rsplus)),	
	sind=csyp_nz$species,
	xpred=rs_prcp_lp$xarr[,2,] + log(tau_p), # 2 = prcp
	ypred=log(rs_prcp_quant),
	xptype="matrix",
	pcol=myblue,
	xname=expression(z(t)),
	yname="ln(Y(t)|Y(t)>1)",
	xlim=c(3.3,5.8)
	)

# PCR 

preddatplot(
	plotname="pcr_prcp",
	xdat=myjitter(log(csyp_nz$prcp)),
	ydat=myjitter(log(csyp_nz$pcrplus)),	
	sind=csyp_nz$species,
	xpred=pcr_prcp_lp$xarr[,2,] + log(tau_p), # 2 = prcp
	ypred=log(pcr_prcp_quant),
	xptype="matrix",
	pcol=myblue,
	xname=expression(z(t)),
	yname=expression(ln(Y(t))),
	xlim=c(3.3,5.8)
	)

### DENS 

# PR

preddatplot(
	plotname="pr_dens",
	xdat=myjitter(log(csyp_nz$germd)),
	ydat=myjitter(csyp_nz$pr),	
	sind=csyp_nz$species,
	xpred=pr_dens_lp$xarr[,4,] + log(tau_d), # 2 = prcp # 4 because 4d with sq term
	ypred=pr_dens_quant,
	xptype="matrix",
	pcol=myred,
	xname=expression(ln(N[g](t))),
	yname=expression(Pr(Y(t)>1)),
	ylim=c(0,1)
	)

# RS

preddatplot(
	plotname="rs_dens",
	xdat=myjitter(log(csyp_nz$germd)),
	ydat=myjitter(log(csyp_nz$rsplus)),	
	sind=csyp_nz$species,
	xpred=rs_dens_lp$xarr[,4,] + log(tau_d), # 2 = prcp
	ypred=log(rs_dens_quant),
	xptype="matrix",
	pcol=myred,
	xname=expression(ln(N[g](t))),
	yname="ln(Y(t)|Y(t)>1)"
	)

# PCR 

preddatplot(
	plotname="pcr_dens",
	xdat=myjitter(log(csyp_nz$germd)),
	ydat=myjitter(log(csyp_nz$pcrplus)),	
	sind=csyp_nz$species,
	xpred=pcr_dens_lp$xarr[,4,] + log(tau_d), # 3 = dens
	ypred=log(pcr_dens_quant),
	xptype="matrix",
	pcol=myred,
	xname=expression(ln(N[g](t))),
	yname=expression(ln(Y(t)))
	)

#######################
### PROPORTIONAL DD ###
#######################

xarr_largedens <- array(NA,dim=c(nseq,3,22))
xarr_largedens[] <- c(
	rep(1,nseq),
	rep(mean(msy$prcp/tau_p,na.rm=T),length.out=nseq), 
	seq(
		from=log(min(0.01*msy$germdhat[msy$germdhat>0]/tau_p,na.rm=T)),	
		to=log(max(100*msy$germdhat[msy$germdhat>0]/tau_p,na.rm=T)),
		length.out=nseq
		)
	)

pr_largedens_lp <- r_lp(pr$beta_p,"pr",xtype="data",xarr=xarr_largedens)
pr_largedens_quant <- plogis(r_quant(pr_largedens_lp))
rs_largedens_lp <- r_lp(rs$beta_r,"rs",xtype="data",xarr=xarr_largedens)
rs_largedens_quant <- exp(r_quant(rs_largedens_lp))
pcr_largedens_lp <- list(
	xarr = pr_largedens_lp$xarr,
		# pr = arbitrary
	pred_lp = log( plogis(pr_largedens_lp$pred_lp) * exp(rs_largedens_lp$pred_lp) )
	)
pcr_largedens_quant <- exp(r_quant(pcr_largedens_lp))

nn_largedens_quant <- array(NA,dim=dim(pcr_largedens_quant))

for(i in 1:3){
	nn_largedens_quant[,i,] <- log(pcr_largedens_quant[,i,]) + xarr_largedens[,3,] + log(tau_d)
	}

preddatplot(
	plotname="pcr_prop_largedens",
	xdat=myjitter(log(csyp_nz$germd)),
	ydat=myjitter(log(csyp_nz$pcrplus*csyp_nz$germd)),	
	sind=csyp_nz$species,
	xpred=xarr_largedens[,3,] + log(tau_d), # 3 = dens
	ypred=nn_largedens_quant,
	xptype="matrix",
	pcol=myred,
	xname=expression(ln(N[g](t))),
	yname=expression(ln(N[n](t))),
	oneoneline=T,
	xlim=c(-5,15),
	ylim=c(-10,20)
	)


