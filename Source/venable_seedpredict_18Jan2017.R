###################################################
### Exploring observation-level variance models ###
###################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"
# setwd("/mnt/data/home/NIOO/calluml/Source/")

library(plyr)
library(reshape2)
library(matrixStats)
library(truncdist)

# setwd(paste0(maindir,"Analyses/Venable"))
# source("venable_figure_functions_23Apr2016.R")

csy <- read.csv("csy_15Jan2016.csv",header=T)
msy <- read.csv("msy_15Jan2016.csv",header=T)

csyp <- read.csv("csyp_15Jan2016.csv",header=T)
ssyp <- read.csv("ssyp_15Jan2016.csv",header=T)

rsdat <- read.csv("cdpos_15Jan2016.csv",header=T)

#####################
### NEW VARIABLES ###
#####################

tau_d <- 100 	# adjustment for density
tau_p <- 100 	# adjustment for rainfall
tau_s <- 100 # Density per 10 cm square

nsim <- 100

spvals <- levels(msy$species)
keepvars <- c("species","year","plot","area")

gdata <- subset(csyp,is.na(ngerm)==F,
  select=c("species","year","plot","ngerm","pr","nfert","area","prcp","newseed")
  )

exlsp <- with(subset(msy,year<2001),
  names(which(tapply(totlive,species,sum,na.rm=T)<5))
  )
  # species to be omitted

odata <- subset(ssyp,is.na(totlive)==F
    & (
      ((species %in% exlsp)==T & year>=2001)
      |
        ((species %in% exlsp)==F & year>=1990)
    ),
    select=c("species","year","plot","totlive","area_o")
  )

prdat <- subset(csyp, ngerm>0 & is.na(nmeas)==F & is.na(germd/tau_d)==F)

spvals <- levels(msy$species)
# yrvals_g <- levels(as.factor(gdata$year))
# yrvals_o <- levels(as.factor(odata$year))
# 
# gdata$year <- as.factor(gdata$year)
# odata$year <- as.factor(odata$year)
# gdata$spyear <- with(gdata, as.numeric(factor(species:year)))
# odata$spyear <- with(odata, as.numeric(factor(species:year)))
# 	# for G and O models spyear was fitted separately
# gdata$spsite_g <- with(gdata, as.numeric(factor(species:plot)))
# odata$spsite_o <- with(odata, as.numeric(factor(species:plot)))
# 
# go <- readRDS("gnzhh_onhh_pars_medians_26Oct2015.rds")
# go2 <- readRDS("onhh_pars_medians_naspecies_28Mar2016.rds")
#   # separate omod run, omitting na species
pr <- readRDS("pr_pars_yearhet_squared_pc_02Mar2016.rds")
rs <- readRDS("rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")

beta_p <- apply(pr$beta_p,c(2,3),median)
beta_r <- apply(rs$beta_r,c(2,3),median)
sig_s_p <- median(pr$sig_s_p) # only one variance param for p
sig_s_r <- median(rs$sig_s_r) # also only one site variance param for r
sig_o_p <- median(pr$sig_o)
phi_r <- median(rs$phi)	
eps_y_p <- pr$eps_y_p
eps_y_r <- rs$eps_y_r
eps_s_r <- rs$eps_s_r
	# pr predictions not needed because pr obs only 
  # missing from 1% of ngerm>0 plots

# sig_s_g <- colMedians(go$sig_s_g)
# sig_s_o <- colMedians(go2$sig_s_o)
# 
# theta_g <- colMedians(go$theta_g)
# 
# phi_g <- median(go$phi_g)
# phi_o <- median(go2$phi_o)

# L <- nlevels(msy$species)
# I_g <- nrow(gdata)
# I_o <- nrow(odata)
# J_g <- length(unique(gdata$year))
# J_o <- length(unique(odata$year))
# K_g <- length(unique(gdata$spsite_g))
# K_o <- length(unique(odata$spsite_o))

# gdata$species_num <- match(gdata$species,levels(msy$species))
# odata$species_num <- match(odata$species,levels(msy$species))
prdat$species_num <- match(prdat$species,levels(msy$species))

prdat$spsite_p <- with(prdat, as.numeric(factor(species:plot)))
rsdat$spsite_r <- with(rsdat, as.numeric(factor(species:plot)))

msy$spyear_m <- with(msy, as.numeric(factor(species:as.factor(year))))

msy_merge <- subset(msy,select=c("species","year","spyear_m"))
# gdata <- merge(gdata,msy_merge,by=c("species","year"),all.x=T)
prdat <- merge(prdat,msy_merge,by=c("species","year"),all.x=T)
rsdat_merge <- unique(subset(rsdat,select=c("species","year","plot","spsite_r")))
prdat <- merge(prdat,rsdat_merge,by=c("species","year","plot"),all.x=T)
	# spsite p not needed for g, because either
	# encapsulated in lam or sites sampled anew

# ### CALCULATIONS
# 
# loglam_g_marg <- site_g <- matrix(nr=I_g,nc=nsim)
# loglam_o_marg <- site_o <- matrix(nr=I_o,nc=nsim)
# 
# eps_s_g <- matrix(nr=K_g,nc=nsim)
# eps_s_o <- matrix(nr=K_o,nc=nsim)
# 
# loglam_g_marg[] <- go$loglam_g_marg
# loglam_o_marg[] <- go2$loglam_o_marg
# 
# eps_species_g <- gdata$species_num[match(unique(gdata$spsite_g),gdata$spsite_g)]
# eps_species_o <- odata$species_num[match(unique(odata$spsite_o),odata$spsite_o)]
# 
# eps_s_g[] <- rnorm(K_g*nsim,0,sig_s_g[rep(eps_species_g,nsim)])
# eps_s_o[] <- rnorm(K_o*nsim,0,sig_s_o[rep(eps_species_o,nsim)])
# 
# site_g[] <- eps_s_g[gdata$spsite_g,] 
# site_o[] <- eps_s_o[odata$spsite_o,] 
# 
# lam_g <- exp(loglam_g_marg + site_g)
# lam_o <- exp(loglam_o_marg + site_o)
# 
# ### PREDICT G AND O DENS
# 
# # G
# 
# lam_g_molt <- melt(lam_g)
# names(lam_g_molt) <- c("obsno","simno","lam_g")
# gdat_sim <- data.frame(gdata[rep(1:I_g,times=nsim),],lam_g_molt)
# gdat_sim$theta_g <- theta_g[gdat_sim$species_num]
#   # matrix of median theta_g samples given species matching each row of gdata
# 
# gdat_sim$ngerm_hat <- with(gdat_sim,
# 	rbinom(I_g*nsim,prob=1-theta_g,size=1)
# 	* rnbinom(I_g*nsim,mu=lam_g,size=phi_g)
# 	)
# 	# zinfl, so no zero-truncation
# 
# sum(gdata$ngerm)
# hi <- with(gdat_sim,tapply(ngerm_hat,simno,sum))
# 	# something to do with minor skew in dist?
# 
# # O 
# 
# lam_o_molt <- melt(lam_o)
# names(lam_o_molt) <- c("obsno","simno","lam_o")
# odat_sim <- data.frame(odata[rep(1:I_o,times=nsim),],lam_o_molt)
# 
# odat_sim$totlive_hat <- with(odat_sim, 
# 	rnbinom(I_o*nsim,mu=lam_o,size=phi_o)
# 	)
# # no theta_o here because no zero-inflation for odist
# 
# sum(odata$totlive)
# hi2 <- with(odat_sim,tapply(totlive_hat,simno,sum))

# FILL IN MISSING VALUES WITH ZEROES

prdat$site_r <- eps_s_r[prdat$spsite_r] 
	# site index irrelevant later because sites simulated anyway
	# rs$ needed here to prevent naming conflicts later
prdat$site_r[is.na(prdat$spsite_r)] <- 0
	# just using marginal prediction for missing sites
prdat$year_r <- eps_y_r[prdat$spyear_m]
	# m because using spyear from msy 
  # (checked, still the case)

## RS - calculate true values

xdat_pred <- data.frame(
	intfake=rep(1,nrow(prdat)),
	ltprcp=log(prdat$prcp/tau_p),
  ltprcp2=(log(prdat$prcp/tau_p))^2,
	ltgermd=log(prdat$ngerm/(prdat$area*tau_d))
	)

beta_r_mat_pred <- as.matrix(xdat_pred) * beta_r[prdat$species_num,]
beta_r_marg_pred <- rowSums(beta_r_mat_pred)
prdat$lam_r <- exp(beta_r_marg_pred + prdat$year_r + prdat$site_r)

gposition_pred <- rep(1:nrow(prdat),prdat$nfert)
	# one g position for every fertile germinant in that position
# seed_hat <- prdat$lam_r[gposition_pred]
n_gposition <- length(gposition_pred)
n_seedsim <- 25
seed_sim <- matrix(nr=n_gposition,nc=n_seedsim)
seed_sim[] <- rtrunc(n_gposition*n_seedsim, spec="nbinom",a=0,
	mu=rep(prdat$lam_r[gposition_pred],n_seedsim),size=phi_r
	)
seed_sim_mean <- rowMeans(seed_sim)
	# simulated and took averages as easy way to get expected values
seed_hat <- tapply(seed_sim_mean,gposition_pred,sum)
prdat$seed_hat <- vector(length=nrow(prdat))
prdat$seed_hat[unique(gposition_pred)] <- seed_hat
prdat$isseed_hat <- prdat$seed_hat>0

prdat_agg <- ddply(prdat,.(species,year),summarise,
	totseed_hat = round(sum(seed_hat),0),
	totisseed = sum(seed_hat>0),
	area_r = sum(area[seed_hat>0])
	)
# area_n = mysum(area[(is.na(nfert)==F & is.na(nbarren)==F) | ngerm==0]), 
# (checked, correct and *includes germ counts of 0*)

msy_new <- merge(msy,prdat_agg,by=c("species","year"),all.x=T)
msy_new$totseed_hat[is.na(msy_new$totseed_hat)] <- 0
	# replacing NAs for years in which all pr is NA 
	# (usually zero germinants, but three obs of 1,7,4 germinants)
msy_new$isseed_hat <- msy_new$totseed_hat > 0
msy_new$area_r0 <- with(msy_new,area_n - area_r)
msy_new$totnoseed <- with(msy_new,
  ifelse(is.na(totisseed), nplot_n, nplot_n - totisseed)
  )
  # some msy rows have no associated pr - assume zero seed
msy_new$nsdhat <- with(msy_new,totseed_hat/area_n)

plot(log(nsdhat)~log(nsdbar),data=msy_new)
abline(0,1,col="red")

# ## PR / RS
# 
# sp_bysite_g <- with(gdata, species_num[match(unique(spsite_g),spsite_g)])
# 
# eps_s_p <- rnorm(K_g*nsim,0,sig_s_p) # sig_s all same for pr
# eps_s_r <- rnorm(K_g*nsim,0,sig_s_r)# sig_s all same for rs now too
# 
# gdat_sim$site_p <- eps_s_p[gdat_sim$spsite_g] 
# gdat_sim$site_r <- eps_s_r[gdat_sim$spsite_g] 
# 	# just need to indicate different sites here 
# 	# - don't need same as from pr
# 
# gdat_sim$year_p <- eps_y_p[gdat_sim$spyear_m] 
# gdat_sim$year_r <- eps_y_r[gdat_sim$spyear_m] 	
# 	# year effects were fitted for all combos in msy
# 
# xdat <- data.frame(
# 	int=rep(1,nrow(gdat_sim)),
# 	prcp=log(gdat_sim$prcp/tau_p),
#   prcp2=(log(gdat_sim$prcp/tau_p))^2,
# 	ngerm=log(gdat_sim$ngerm_hat/(gdat_sim$area*tau_d))
# 	)
# 
# beta_p_mat <- as.matrix(xdat) * beta_p[gdat_sim$species_num,]
# beta_p_marg <- rowSums(beta_p_mat)
# gdat_sim$pi <- plogis(beta_p_marg + gdat_sim$year_p + gdat_sim$site_p + rnorm(I_g*nsim,0,sig_o_p))
# 
# beta_r_mat <- as.matrix(xdat) * beta_r[gdat_sim$species_num,]
# beta_r_marg <- rowSums(beta_r_mat)
# gdat_sim$lam_r <- exp(beta_r_marg + gdat_sim$year_r + gdat_sim$site_r)
# 
# gdat_sim$nfert_hat <- with(gdat_sim,rbinom(I_g*nsim,prob=pi,size=ngerm_hat))
# 
# gposition <- rep(1:(I_g*nsim),gdat_sim$nfert_hat)
# library(truncdist)
# seed_hat <- rtrunc(length(gposition),spec="nbinom",a=0,b=Inf,mu=gdat_sim$lam_r[gposition],size=phi_r)
# totseed_hat <- tapply(seed_hat,gposition,sum)
# gdat_sim$totseed_hat <- vector(length=I_g*nsim)
# gdat_sim$totseed_hat[unique(gposition)] <- totseed_hat
# 
# with(gdat_sim,
# 	cbind(tapply(ngerm_hat,species,sum),tapply(ngerm,species,sum))
# 	)
# 
# gdat_sim$isseed_hat <- gdat_sim$totseed_hat > 0
# 
# # assuming pr and rs site effect uncorrelated
# #sum(gdata$ngerm)
# #with(gdat_sim,tapply(ngerm_hat,simno,sum))
# #plot(log(gdat_sim$ngerm_hat)~log(gdat_sim$ngerm))
# #abline(0,1,col="red")
# #lines(supsmu(log(gdat_sim$ngerm),log(gdat_sim$ngerm_hat)),col="blue")
# 
# ### AGGREGATE BY SPECIES AND YEAR
# 
# gdat_agg <- ddply(gdat_sim,.(species,year,simno),summarise,
# 	totgerm_obs = sum(ngerm),
# 	totgerm_hat = sum(ngerm_hat),
# 	totarea = sum(area)
# 	)
#   # "area" = area of plot from csyp
# gdat_agg$germd_hat <- with(gdat_agg,totgerm_hat/totarea)
# 
# gdat_sim_r <- subset(gdat_sim,is.na(pr)==F | ngerm==0)
#   # for reproduction stuff
# 
# sdat_agg <- ddply(gdat_sim_r,.(species,year,simno),summarise,
# 	# totpres_obs = sum() # can't do this unless join predicts from msy_new
# 	totpres_hat = sum(isseed_hat==T),
# 	totabs_hat = sum(isseed_hat==F),
# 	totseed_hat = sum(totseed_hat),
# 	totarea = sum(area),
#   nsd_obs = sum(newseed)/sum(area[is.na(newseed)==F])
# 	)
#   # "area" = area of plot from csyp
# sdat_agg$isseed_hat <- sdat_agg$totseed_hat > 0
#   # reproduction
# sdat_agg$nsd_hat <- with(sdat_agg,totseed_hat/totarea)
# 
# odat_agg <- ddply(odat_sim,.(species,year,simno),summarise,
# 	totlive_hat = sum(totlive_hat),
# 	totarea = sum(area_o)
# 	)
# odat_agg$olsd_hat <- with(odat_agg,totlive_hat/totarea)
# 	# no need to adjust because predictions (lam) already incorporate area
# 	# trying to predict variation in actual count
# 
# 	# fit distribution to these errors and then use that as error?
# 	# process error handled by random G etc...?
# 
# plot(log(nsd_obs)~log(nsd_hat),data=sdat_agg)
# abline(0,1,col="red")
# with(sdat_agg,lines(supsmu(log(nsd_hat),log(nsd_obs)),col="blue"))
# 
# ### FIGURES
# 
# msy$olsdhat <- with(msy,totlive/area_o)
#   # this missing, for some reason (hat -> total count / total seeds)
# 
# pdf(paste0("germd_bootstrap_hists_log_",format(Sys.Date(),"%d%b%Y"),".pdf"),
# 	width=10,height=10)
# 	for(i in 1:L){
# 	par(mfrow=c(6,6),mar=c(3,3,1,1),oma=c(0,0,4,0))
# 		for(j in 1:J_g){
# 			d <- subset(gdat_agg,species==spvals[i] & year==yrvals_g[j])
# 			germd_obs <- msy$germdhat[msy$species==spvals[i] & msy$year==yrvals_g[j]]
# 			hist(log(d$germd_hat),breaks=100,
# 					xlim=c(0,max(log(c(d$germd_hat,germd_obs))))
# 					)
# 			abline(v=log(germd_obs),col="red")
# 			mtext(d$species[1],side=3,outer=T,line=1)
# 			}
# 		}
# dev.off()
# 
# pdf(paste0("olsd_bootstrap_hists_log_",format(Sys.Date(),"%d%b%Y"),".pdf"),
# 	width=10,height=10)
# 	for(i in 1:L){
# 	par(mfrow=c(6,6),mar=c(3,3,1,1),oma=c(0,0,4,0))
# 		for(j in 8:J_g){
# 			d <- subset(odat_agg,species==spvals[i] & year==yrvals_o[j])
# 			if(nrow(d)>0 & sum(d$olsd_hat)>0) {
# 			  olsd_obs <- msy$olsdhat[msy$species==spvals[i] & msy$year==yrvals_o[j]]
# 				hist(log(d$olsd_hat),breaks=100,
# 					xlim=c(0,max(log(c(d$olsd_hat,olsd_obs))))
# 					)
# 				abline(v=log(olsd_obs),col="red")
# 				mtext(d$species[1],side=3,outer=T,line=1)
# 				}
# 			}
# 		}
# dev.off()
# 
# pdf(paste0("nsd_bootstrap_hists_",format(Sys.Date(),"%d%b%Y"),".pdf"),
# 	width=10,height=10)
# 	for(i in 1:L){
# 	par(mfrow=c(6,6),mar=c(3,3,1,1),oma=c(0,0,4,0))
# 		for(j in 1:J_g){
# 			d <- subset(sdat_agg,species==spvals[i] & year==yrvals_o[j])
# 			nsd_obs <- msy_new$newseedhat[msy_new$species==spvals[i] & msy_new$year==yrvals_g[j]]
# 			if(nrow(d)>0 & sum(d$nsd_hat)>0) {
# 				hist(log(d$nsd_hat),breaks=100,
# 					xlim=c(0,log(max(c(d$nsd_hat,nsd_obs))))
# 					)
# 				abline(v=log(nsd_obs),col="red")
# 				mtext(d$species[1],side=3,outer=T,line=1)
# 				}
# 			}
# 		}
# dev.off()
#   # red lines (observed values) may not be correct here because 
#   # don't account for unmeasured plants...?
# 
# ### TOTAL MODELS
# 
# # SIM SDS PER SPECIES
# 
# stddevs <- function(obj){
#     unlist(lapply(VarCorr(obj), function(m) sqrt(diag(m))))
# 	}
# 
# gsds <- sapply(split(gdat_agg,gdat_agg$species),function(d){
# 	tgm <- glmer(totgerm_hat ~ offset(log(totarea*tau_s)) 
# 		+ (1|year/simno),
# 		family="poisson",
# 		data=d
# 		)
# 	obssd <- stddevs(tgm)["simno:year.(Intercept)"]
# 	return(obssd)
# 	})
# 
# osds <- sapply(split(odat_agg,odat_agg$species),function(d){
# 	tlm <- glmer(totlive_hat ~ offset(log(totarea*tau_s)) 
# 		+ (1|year/simno),
# 		family="poisson",
# 		data=d
# 		)
# 	obssd <- stddevs(tlm)["simno:year.(Intercept)"]
# 	return(obssd)
# 	})
# 
# ssds <- sapply(split(subset(sdat_agg,totseed_hat > 0),sdat_agg$species[sdat_agg$totseed_hat>0]),function(d){
# 	tsm_r <- glmer(totseed_hat ~ offset(log(totarea*tau_s))
# 		+ (1|year/simno),
# 		family="poisson",
# 		data=d
# 		)
# 	obssd <- stddevs(tsm_r)["simno:year.(Intercept)"]
# 	return(obssd)
# 	})
# 
# sds <- data.frame(species=levels(msy$species),gsd=gsds,osd=osds,ssd=ssds,row.names=NULL)
# write.csv(sds,
# 	paste0("observation_error_byspecies_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
# 
# # SIM SDS OVER ALL SPECIES
# 
# tgm <- glmer(totgerm_hat ~ offset(log(totarea*tau_s)) + (1|species/year/simno),family="poisson",data=gdat_agg)
# tlm <- glmer(totlive_hat ~ offset(log(totarea*tau_s)) + (1|species/year/simno),family="poisson",data=odat_agg)
# tsm_p <- glmer(cbind(totpres_hat,totabs_hat) ~ (1|species/year/simno),family="binomial",data=sdat_agg)
# 	# doesn't account for different plot areas
# 
# tsm_p2 <- glmer(isseed_hat ~ log(totarea*tau_s) + (1|species/year/simno),
#                 family=binomial(link=cloglog),data=sdat_agg)
#   # offset(area) doesn't work
#   
# tsm_r <- glmer(totseed_hat ~ offset(log(totarea*tau_s)) + (1|species/year/simno),
#                family="poisson",data=subset(sdat_agg,totseed_hat > 0))
# 
# with(subset(sdat_agg,totseed_hat > 0),tapply(log(totseed_hat/(totarea*10^4)),list(species,year),sd))
# 
# ### ZERO-INFLATION
# 
# totseed_ispos <- with(sdat_agg,tapply(totseed_hat,list(species,year),function(x) sum(x>0)))
# 	# not always zero or non-zero
# totseed_medpos <- with(subset(sdat_agg,totseed_hat>0),tapply(totseed_hat,list(species,year),function(x) mean(log(x))))
# 	# median | >0 - not calculable if all are zeroes
# totseed_sdpos <- with(subset(sdat_agg,totseed_hat>0),tapply(totseed_hat,list(species,year),function(x) sd(log(x))))
# 
# plot(totseed_ispos ~ log(totseed_medpos))
# plot(totseed_ispos ~ log(totseed_sdpos))
# 
# cor.test(totseed_ispos,log(totseed_medpos))
# cor.test(totseed_ispos[is.na(log(totseed_sdpos))==F],log(totseed_sdpos[is.na(log(totseed_sdpos))==F]))
# 
# t(apply(totseed_sdpos,1,mean,na.rm=T))
# 
# hist(totseed_ispos,breaks=100)
# sum(totseed_ispos %in% c(0:5,95:100))/length(totseed_ispos)
# 
# # OTHER
# 
# 
# totgerm_ispos <- with(gdat_agg,tapply(totgerm_hat,list(species,year),function(x) sum(x>0)))
# totlive_ispos <- with(odat_agg,tapply(totlive_hat,list(species,year),function(x) sum(x>0)))
# 
# hist(log(gdat_agg$totgerm_hat),breaks=100)
# hist(log(odat_agg$totlive_hat),breaks=100)
# hist(log(sdat_agg$totseed_hat),breaks=100)
# 
# with(csyp,hist(germd,breaks=10^4,xlim=c(0,500)))
# with(ssyp,hist(olsdbar,breaks=10^4,xlim=c(0,1000)))
# 
# # FRACTION OF UNSUITABLE SITES (WITH NO DENSITIES >0)
# 
# big0sum <- function(x) sum(x>0,na.rm=T)
# big0prop <- function(x) sum(x>0)/length(x)
# 
# g0 <- with(csyp,tapply(germd,list(plot,species),big0sum))
# apply(g0,2,big0prop)
# # na.rm because some plots missing in some years, so have NA germd values
# # NB: germdhat is for ALL plots for that species in that year
# # so need germd instead
# 
# s0 <- with(ssyp,tapply(olsdbar,list(plot,species),big0sum))
# apply(s0,2,big0prop)
# 
# ### CALC SDS
# 
# loglam_g_sd <- ddply(gdat_agg,.(species,year),summarise,
# 	sd_totgerm = sd(log(totgerm_hat))
# 	)
# loglam_s_sd <- ddply(gdat_agg,.(species,year),summarise,
# 	sd_totseed = sd(log(totseed_hat))
# 	)
# loglam_o_sd <- ddply(odat_agg,.(species,year),summarise,
# 	sd_totlive = sd(log(totlive_hat))
# 	)

### SAVE DATA

# OBSERVED DATA

write.csv(msy_new,
	paste0("msy_seedests_",format(Sys.Date(),"%d%b%Y"),".csv"),
	row.names=F
	)

# write.csv(prdat,
#   paste0("pr_seedests_",format(Sys.Date(),"%d%b%Y"),".csv"),
#   row.names=F
#   )
# 
# # SIMULATED DATA
# 
# write.csv(gdat_agg,
# 	paste0("gdat_agg_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
# write.csv(sdat_agg,
# 	paste0("sdat_agg_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
# write.csv(odat_agg,
# 	paste0("odat_agg_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
# 
# # SDS
# 
# write.csv(loglam_g_sd,
# 	paste0("loglam_g_sd_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
# write.csv(loglam_o_sd,
# 	paste0("loglam_o_sd_",format(Sys.Date(),"%d%b%Y"),".csv"),
# 	row.names=F
# 	)
