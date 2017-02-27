######################################################################
### Simulations of population dynamics for infinite area           ###
### Includes purely analytical calculations of DI long-term growth ###
######################################################################

# *Update file directories*

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

library(plyr)
library(reshape2)
library(mvtnorm) # for MVN

### LOAD DATA 

setwd(paste0(maindir,"Analyses/Venable"))

source("venable_figure_functions_23Apr2016.R")
source("venable_prediction_functions_06Dec2016.R")

msy <- read.csv("msy_seedests_28Mar2016.csv",header=T)
Tvalues <- read.csv("Tvalues_31Jul2015.csv",header=T)

#####################
### DEFINE PARAMS ###
#####################

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
nyear <- length(unique(msy$year))

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])
	# T in years

tau_s <- 100 # adujstment for seeds
tau_d <- 100	# adjustment for density
tau_p <- 100	# adjustment for rainfall

ncy <- read.csv("ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
# removes first value (missing because no previous winter)

msy$year_num <- as.numeric(as.character(msy$year))
exlsp <- with(subset(msy,year<2001),
  names(which(tapply(totlive,species,sum,na.rm=T)<5))
  )
msy90 <- subset(msy,
  ((species %in% exlsp)==T & year_num>=2001)
  |
  ((species %in% exlsp)==F & year_num>=1990)
  )
msy90$year <- as.factor(as.character(msy90$year))

# Assumptions:
# - plant and seed plant densities stay constant at median values
# - current rainfall conditions (prcp,gprcp,cor)
# - no spatial variation in reproduction
# - (spatial variation in germination irrelevant)
# - no random interannual variation in any variable
# - constant mean rainfall; no variance (add later)

###################
### QUICK PLOTS ###
###################

### SIMULATE CLIMATE

z_sig <- sd(log(ncy$seasprcp))
z_mu <- mean(log(ncy$seasprcp))
z_min <- min(log(ncy$seasprcp))
z_max <- max(log(ncy$seasprcp))

w_lm <- glm(log(germprcp) ~ log(seasprcp),data=ncy)

# w_sig <- sd(log(ncy$germprcp))
# w_mu <- mean(log(ncy$germprcp))
  # should predict using lm instead

# IGNORING: add when have variability
# rho <- 0.82
# 
# z <- w <- array(NA,c(ni,nt))
# 
# zw_mu <- c(z_mu,w_mu)
# zw_sig <- matrix(c(z_sig^2,rep(rho*z_sig*w_sig,2),w_sig^2),nr=2,nc=2)
# 
# zw <- dmvnorm(x=, mean=zw_mu, sigma=zw_sig)
# 
# z[] <- zw[,1] - log(tau_p)
# # transformed log winter rainfall
# w[] <- zw[,2] - log(tau_p)
# # transformed log germination rainfall

### DO THE REST

nseq <- 10^3
ni <- 10^4

z <- seq(z_min,z_max,length.out=nseq) - log(tau_s)
z_rel <- z + log(tau_s) - z_mu
  # = log(prcp/mean prcp)
# w <- w_mu + z_rel - log(tau_s)
  # w values scaled in same way as g values 
  # i.e. assuming 1:1 change in means of z and w
# w <- rep(w_mu - log(tau_s),nseq) # ADDED TAU HERE
zp <- exp(z+log(tau_s))
w <- predict(w_lm,newdata=list(seasprcp=zp),type="response") - log(tau_s)
int <- rep(1,nseq)

loglam_q <- array(dim=c(nseq,3,nspecies))

for(j in 1:nspecies){
  curmsy <- subset(msy90,species==spvals[j])
  lsd <- with(curmsy, 
    rep(log(median(germdhat+totlive/area_o,na.rm=T)),nseq) - log(tau_s)
    )
    # repeat this for different density values
  x <- cbind(int,z,z^2,lgd)
  
  loglam <- matrix(NA,nr=nseq,nc=ni)
  for (i in 1:ni){
    beta_p <- pr$beta_p[i,j,]
    beta_r <- rs$beta_r[i,j,]
    phi_r <- rs$phi_r[i] # needed to calculate truncation
    alpha_G <- go$alpha_G[i,j]
    beta_Gz <- go$beta_Gz[i,j]
    alpha_m <- go$alpha_m[i,j]
    beta_m <- go$beta_m[i,j]
    m0 <- exp(alpha_m)
    m1 <- exp(beta_m)
    
    G <- Gcalc(w,alpha_G,beta_Gz)
    p <- t(plogis(beta_p %*% t(x)))
    r0 <- t(exp(beta_r %*% t(x)))
    r <- nbtmean(r0,phi_r)
    Y <- p*r
    So <- exp(-m0)
    
    lgd <- lsd + log(G)
      
    Sn <- Sncalc(nd,m0,m1,T3) # over different time interval to So
      # need to change this too
    
    loglam[,i] <- log( G*Y*Sn + (1-G)*So )
    }
  
  loglam_q[,,j] <- t(apply(loglam,1,quantile,probs=c(0.025,0.5,0.975),na.rm=T))
  }

zrange <- quantile(log(ncy$seasprcp),probs=c(0.25,0.75))
zcv_new <- with(ncy,1.2*sd(seasprcp)/mean(seasprcp))
zam_new <- with(ncy,0.8*mean(seasprcp))
z_sig_new <- sqrt(log(zcv_new^2+1)) # has to be done first
z_mu_new <- log(zam_new) - z_sig_new^2/2
  
zrange <- quantile(log(ncy$seasprcp),probs=c(0.25,0.75))
zrange_new <- qnorm(c(0.25,0.75),mean=z_mu_new, sd=z_sig_new)

rpol <- function(zrange,yrange,polcol){
  polygon(x=c(zrange,rev(zrange)),
    y=rep(yrange,each=2),
    col=polcol,
    border=NA
    )
  }

rgrowplots <- function(y,xname="ln rainfall",yname="Population growth rate", npol=0,...){
  par(mfrow=c(4,6),mar=c(2.5,2.5,4,2.5),oma=c(4,4,1,1),las=1,bty="l")
  for(i in 1:nspecies){
    g <- G_ord[i]
    matplot(I(z+log(tau_s)),y[,,g],
      lty=c(2,1,2),
      col="black",
      ...
      )
    yrange <- c(min(y),max(y))
    if(npol>0) rpol(zrange,yrange,tranteal)
    if(npol>1) rpol(zrange_new,yrange,tranbrown)
    
    mtext(bquote(bold((.(letters[i])))~.(spvals[g])),side=3,line=0.5,xpd=T,adj=0.05)
    
    if(i %in% 17:23) addxlab(xname) 
    if(i %in% seq(1,23,6)) addylab(yname) 
    }
  }

msy$olsdhat <- with(msy,totlive/area_o)
msy$godhat <- with(msy,germdhat+olsdhat)
G_med <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,median,na.rm=T)
  )
G_ord <- order(G_med)

tranteal <- rgb(226,240,214,maxColorValue=255,alpha=150)
tranbrown <- rgb(139,69,19,maxColorValue=255,alpha=25)

pdf(paste0("allspecies_growthcurves_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=14,height=10)
  rgrowplots(loglam_q, npol=0, type="n")
  rgrowplots(loglam_q, npol=0, type="l")
  rgrowplots(loglam_q, npol=1, type="l")
  rgrowplots(loglam_q, npol=2, type="l")
dev.off()

rgrowplot_scba <- function(npol=0,...){
  par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")
  y <- loglam_q[,,spvals=="scba"]
  matplot(I(z+log(tau_s)),y,
    lty=c(2,1,2),
    col="black",
    xlab="ln rainfall",
    ylab="Population growth rate",
    ...
  )
  yrange <- c(min(y),max(y))
  if(npol>0) rpol(zrange,yrange,tranteal)
  if(npol>1) rpol(zrange_new,yrange,tranbrown)
  }

pdf(paste0("scba_growthcurve_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=4.5,height=4.5)
  rgrowplot_scba(npol=0,type="n")
  rgrowplot_scba(npol=0,type="l")
  rgrowplot_scba(npol=1,type="l")
  rgrowplot_scba(npol=2,type="l")
dev.off()

########################################################################

# ln_cons <- function(mean){
# 	mzd <- abs(z - log(mean))
# 	mz <- which(mzd==min(mzd))
# 	pz <- rep(0,nseq) 
#  	pz[mz] <- 1
# 	return(pz)
# 	}
# 
# ln_var <- function(mean,cv){
# 	sig <- sqrt(log(cv^2+1))
# 	mu <- log(mean) - sig^2/2
# 	pz <- dnorm(z,mean=mu,sd=sig)
# 	return(pz)
# 	}

m <- 2
	# multiplicative factor

pz <- cbind(
	pz_m1c0 = ln_cons(z_mean),
	pz_m2c0 = ln_cons(z_mean*m),
	pz_m1c1 = ln_var(z_mean,z_cv),
	pz_m2c1 = ln_var(z_mean*m,z_cv),
	pz_m1c2 = ln_var(z_mean,cv*m),
	pz_m2c2 = ln_var(z_mean*m,z_cv*m),
	pz_m05c0 = ln_cons(z_mean/m),
	pz_m05c1 = ln_var(z_mean/m,z_cv),
	pz_m05c2 = ln_var(z_mean/m,z_cv*m)
	)

rbar <- matrix(nr=ncol(pz),nc=22,dimnames=list(colnames(pz),levels(msy$species)))
for(i in 1:22){
	rbar[,i] <- ( t(pz) %*% rt[,i] ) / colSums(pz)
	}

m1hc0 <- rbar["pz_m1c1",] - rbar["pz_m1c0",]
m1hc <- rbar["pz_m1c2",] - rbar["pz_m1c1",]
hmc1 <- rbar["pz_m2c1",] - rbar["pz_m1c1",]
lmc1 <- rbar["pz_m05c1",] - rbar["pz_m1c1",]
hmhc <- rbar["pz_m2c2",] - rbar["pz_m1c1",]
lmhc <- rbar["pz_m05c2",] - rbar["pz_m1c1",]

par(mfrow=c(3,2),mar=c(4,4,1,1))
plot(m1hc0~G_med)
plot(m1hc~G_med)
plot(lmc1~G_med)
plot(lmhc~G_med)
plot(hmc1~G_med)
plot(hmhc~G_med)

par(mfrow=c(1,2))
plot(m1hc~G_med,pch=16)
plot(lmc1~G_med,pch=16)

mycols <- colorRampPalette(c("saddlebrown", "palegreen"))

pdf(paste0("mean_variance_comparison_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=5.3,height=5.3)

par(mfrow=c(1,1),mar=c(8,8,2,2),las=1,bty="l")
plot(m1hc~lmc1,
	bg=mycols(22)[rank(G_med)],
	col="black", 
	pch=21,
	cex=1.5,
	ann=F
	)
mtext("Effect of decreased mean rainfall \n on population growth rate",side=1,line=4,cex=1.1)
mtext("Effect of increased rainfall variability \n on population growth rate",side=2,line=4,las=3,cex=1.1)

dev.off()

### SCBA PLOTS FOR OSLO TALK
# run go and pr/rs model fits first
# venable_outputgraphs_Stan_GO_naspecies_noerr_GDD_07Apr2016
# venable_outputgraphs_Stan_PrRs_07Mar2016

tranteal <- rgb(190,235,159,maxColorValue=255,alpha=100)
j <- which(spvals=="scba")
stanplot <- function(xd,yd,xp,yp,...){
  plot(yd~xd,
    bg=tranteal,
    col="black",
    pch=21,
    cex=1.5,
    cex.lab=1.4,
    cex.axis=1.2,
    ...
    )
  matplot(xp,yp,add=T,type="l",lty=c(2,1,2),col="black")
  }

pdf(paste0("scba_relationships_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=4.5,height=4.5)
par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")

stanplot(
  xd=log(msy90$prevnsdhat[msy90$species=="scba"]),
  yd=log(msy90$Shat1[msy90$species=="scba"]),
  xp=lnsdseqs[,j],
  yp=log(Squant_new[,,j]),
  xlab="ln seed density",
  ylab="ln survival"
  )
stanplot(
  xd=msy90$tprcp[msy90$species=="scba"] + log(tau_p),
  yd=msy90$Ghat[msy90$species=="scba"],
  xp=zseq + log(tau_p),
  yp=Gzquantadj[,,j],
  xlab="ln rainfall",
  ylab="Germination probability"
  )
stanplot(
  xd=log(msy90$godhatadj[msy90$species=="scba"]),
  yd=msy90$Ghat[msy90$species=="scba"],
  xp=dmatseq[,j]+log(tau_d),
  yp=Gdquant[,,j],
  xlab="ln seed density",
  ylab="Germination probability"
  )
stanplot(
  xd=myjitter(log(csyp_nz$prcp[csyp_nz$species=="scba"])),
  yd=myjitter(log(csyp_nz$pcrplus[csyp_nz$species=="scba"])),
  xp=pcr_prcp_lp$xarr[,2,j] + log(tau_p),
  yp=log(pcr_prcp_quant[,,j]),
  xlab="ln rainfall",
  ylab="ln seed production"
  )
stanplot(
  xd=myjitter(log(csyp_nz$germd[csyp_nz$species=="scba"])),
  yd=myjitter(log(csyp_nz$pcrplus[csyp_nz$species=="scba"])),
  xp=pcr_dens_lp$xarr[,4,j] + log(tau_d),
  yp=log(pcr_dens_quant[,,j]),
  xlab="ln plant density",
  ylab="ln seed production"
  )

dev.off()
