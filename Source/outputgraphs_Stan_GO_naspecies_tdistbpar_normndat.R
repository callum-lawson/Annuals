####################################################
### Population-level models for Pr, Rs, G, and O ###
####################################################

library(plyr)
library(matrixStats)
library(reshape2)

source("Source/figure_functions.R")
source("Source/prediction_functions.R")

msy <- read.csv("Output/msy_seedests_28Mar2016.csv",header=T)
# no sds attached yet
Tvalues <- read.csv("Output/Tvalues_31Jul2015.csv",header=T)
obserr <- read.csv("Output/observation_error_byspecies_28Mar2016.csv",header=T)

nspecies <- nlevels(msy$species)
spvals <- levels(msy$species)

################
### ADD DATA ###
################

# Taken from:
# venable_Stan_poplevel_binomialG_tdistbpar_normndat_allm1prior_16Feb2016

tau_s <- 100		# adjustment for seed density
tau_p <- 100		# adjustment for rainfall
tau_d <- 100

msy$year_num <- as.numeric(as.character(msy$year))
msyerr <- merge(msy,obserr,by="species")
exlsp <- with(subset(msy,year<2001),
  names(which(tapply(totlive,species,sum,na.rm=T)<5))
  )
msy90 <- subset(msyerr,
  ((species %in% exlsp)==T & year_num>=2001)
  |
  ((species %in% exlsp)==F & year_num>=1990)
  )
msy90$year <- as.factor(as.character(msy90$year))
msy90pos <- subset(msy90,isseed_hat)

dl <- list()

dl$I <- nrow(msy90)
dl$N <- nrow(msy90pos)
# all years/species have germinant surveys
# surveys start in 1990, so first year with prevolsd is 1991

dl$year <- as.numeric(msy90$year)
dl$species <- as.numeric(msy90$species)
dl$tprcp <- log(msy90$gprcp/tau_p)

dl$J <- length(unique(msy90$year))
# effectively obslevel ranef
dl$J0 <- with(msy90,tapply(as.numeric(year),species,min))
dl$L <- max(dl$species)

dl$totgerm <- msy90$totgerm
dl$totlive <- msy90$totlive
dl$isseed <- as.numeric(msy90$isseed_hat)
dl$totseed <- msy90pos$totseed_hat

dl$larea_g <- log(msy90$area_g * tau_s)
dl$larea_o <- log(msy90$area_o * tau_s)
dl$larea_n <- log(msy90pos$area_n * tau_s)

dl$larea_go <- with(dl, larea_g - larea_o)

dl$totlivemult <- with(dl, round(totlive*exp(larea_go),0))
dl$totboth <- with(dl,totgerm+totlivemult)
# area for both = area_g

dl$totboth_p <- with(dl,totboth[totboth>0])
dl$totgerm_p <- with(dl,totgerm[totboth>0])
dl$M <- length(dl$totboth_p)
dl$mpos <- which(dl$totboth>0)

dl$sig_b <- msy90$sig_b 
dl$nu_b <- msy90$nu_b 
dl$sig_n <- msy90pos$sig_n 

dl$T1 <- with(Tvalues,duration[period=="T1"])
dl$T2 <- with(Tvalues,duration[period=="T2"])
dl$T3 <- with(Tvalues,duration[period=="T3"])
# T in years

msy90pos$npos <- 1:nrow(msy90pos)
msy90B <- merge(msy90,msy90pos,by=c("species","year"),all.x=T,all.y=T)
dl$npos <- msy90B$npos
dl$npos[is.na(dl$npos)] <- 0 # placeholder

dl$lseeddens <- with(msy90pos, log(totseed_hat / (area_n * tau_s)))

#################################
### CALCULATE EXTRA VARIABLES ###
#################################

msy90$tprcp <- with(msy90,log(gprcp/tau_p))
msy90$olsdhat <- with(msy90,totlive/area_o)
msy90$Ghat <- with(msy90,germdhat/(olsdhat+germdhat))
msy90$nsdhat <- with(msy90,totseed_hat/area_n)
msy90$csdhat <- with(msy90,olsdhat+nsdhat)

msy90$prevyear <- as.numeric(as.character(msy90$year)) - 1
ms90 <- ddply(msy90, .(species), summarize,
	species = species,
	year = year,
	prevolsdhat = olsdhat[match(prevyear,year)],
	prevnsdhat = nsdhat[match(prevyear,year)],
	germdhat0 = ifelse(germdhat==0,min(germdhat[germdhat>0]),germdhat),
	olsdhat0 = ifelse(olsdhat==0,min(olsdhat[olsdhat>0]),olsdhat)
	)
	# 0 -> replace zeroes with min non-zero obs

msy90 <- merge(msy90,ms90,by=c("species","year"))

msy90$prevcsdhat <- with(msy90,prevolsdhat+prevnsdhat)

msy90$Shat <- with(msy90,(germdhat+olsdhat)/(prevolsdhat+prevnsdhat))
msy90$Shat1 <- with(msy90,ifelse(Shat>1,1,Shat))

#####################
### PLOT RAW DATA ###
#####################

pdf(paste0("Plots/germhat_olsdhat_equalscale_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=7,height=7)
with(msy90,plot(log(germdhat+1)~log(olsdhat+1),
  xlim=c(-2,10),ylim=c(-2,10))
  )
abline(0,1,col="red")
abline(h=0,col="blue",lty=3)
abline(v=0,col="blue",lty=3)
dev.off()

pdf(paste0("Plots/whist_ngerm0_",format(Sys.Date(),"%d%b%Y"),".pdf"),width=5,height=5)
par(mfrow=c(1,1),bty="l")
mybreaks <- with(msy90,hist(log(gprcp),breaks=50,main=""))
with(subset(msy90,germdhat==0),
  hist(log(gprcp),breaks=mybreaks$breaks,add=T,col="blue")
  )
with(subset(msy90,germdhat==0 & olsdhat==0),
  hist(log(gprcp),breaks=mybreaks$breaks,add=T,col="red")
  )
legend("topleft",
  fill=c("white","blue","red"),
  legend=c("All","Ng>0","Ng>0 & No=0"),
  bty="n"
  )
dev.off()

###############################
### CALCULATE RELATIONSHIPS ###
###############################

### DATA

niter <- length(go$alpha_G_mu) # arbitrary
nsp <- dl$L
nseq <- 100
probs <- c(0.025,0.50,0.975)

### MORTALITY

dseqcalc <- function(d,type){
	# rat <- with(subset(d,nsdhat>0 | olsdhat>0), mean(nsdhat/(nsdhat+olsdhat)))
	rat <- 1
	subd <- subset(d, nsdhat > 0 & olsdhat > 0 )
	lcsd <- log(subd$csdhat)
	lcsdseq <- seq(min(lcsd),max(lcsd),length.out=nseq)
	if(type=="old"){
		return( (1-rat)*exp(lcsdseq) / tau_s )
		}
	if(type=="new"){
		return( rat*exp(lcsdseq) / tau_s )
		}
	}

dseqcalc_comb <- function(d){
	with(d, seq(from=log(min(csdhat[csdhat>0])),to=log(max(csdhat)),length.out=nseq))
	}
dseqcalc_new <- function(d){
	with(d, seq(from=log(min(nsdhat[nsdhat>0])),to=log(max(nsdhat)),length.out=nseq))
	}
dseqcalc_old <- function(d){
	with(d, seq(from=log(min(olsdhat[olsdhat>0])),to=log(max(olsdhat)),length.out=nseq))
	}

lolsdseqs <- sapply(split(msy90,msy90$species),dseqcalc_old)
lnsdseqs <- sapply(split(msy90,msy90$species),dseqcalc_new)

m0 <- exp(go$alpha_m)
m1 <- exp(go$beta_m)

Spred_old <- Spred_new <-array(dim=c(niter,nseq,nsp))

for(i in 1:niter){
	for(j in 1:nsp){
		Spred_old[i,,j] <- Scalc(exp(lolsdseqs[,j]),0,m0[i,j],m1[i,j])
		Spred_new[i,,j] <- Scalc(0,exp(lnsdseqs[,j]),m0[i,j],m1[i,j])
		}
	}

Squant_old <- apply(Spred_old,c(2,3),quantile,probs=probs)
Squant_new <- apply(Spred_new,c(2,3),quantile,probs=probs)

Squant_old <- aperm(Squant_old,c(2,1,3))
Squant_new <- aperm(Squant_new,c(2,1,3))
# transpose to fit graph function

### GERMINATION

# RAINFALL

m0_med <- exp(colMedians(go$alpha_m))
msy90$godhat <- with(msy90,germdhat+olsdhat)
msy90$godhatadj <- NA
for(i in 1:nsp){
  select <- msy90$species==spvals[i]
  msy90$godhatadj[select] <- with(msy90[select,],
    germdhat+olsdhat*exp(m0_med[i]*dl$T1)
    )
  }

zseq <- with(dl,seq(min(tprcp),max(tprcp),length.out=nseq))

dmatcons <- matrix(NA,nr=nseq,nc=nspecies)
for(i in 1:nsp){
  dmatcons[,i] <- with(subset(msy90,species==spvals[i]),
    rep(median(log(godhat)-log(tau_s)),nseq)
    )
  } # still not adjusted

Gzpred <- Gzpredadj <- array(dim=c(niter,nseq,nsp))
	# Galt = adjusted for old seed mortality

for(i in 1:niter){
	for(j in 1:nsp){
		Gzpred[i,,j] <- Gcalc(zseq,dmatcons[,j],
		  go$alpha_G[i,j],go$beta_Gz[i,j] # ,go$beta_Gd[i,j]
		  )
		Gzpredadj[i,,j] <- Gzpred[i,,j]/( Gzpred[i,,j] + (1-Gzpred[i,,j])*exp(-m0[i,j]*dl$T1) )
		}
	}

Gzquant <- apply(Gzpred,c(2,3),quantile,probs=probs)
Gzquantadj <- apply(Gzpredadj,c(2,3),quantile,probs=probs)
Gzquant <- aperm(Gzquant,c(2,1,3)) 
Gzquantadj <- aperm(Gzquantadj,c(2,1,3)) 
# transpose to fit graph function

# DENSITY
# 
# Gdpred <- array(dim=c(niter,nseq,nsp))
#   # adjusted for old seed mortality
# 
# dmatseq <- matrix(NA,nr=nseq,nc=nspecies)
# 
# for(j in 1:nsp){
#   curdat <- subset(msy90,species==spvals[j])
#   zcons <- with(curdat,rep(mean(tprcp),length.out=nseq))
#   dvals <- log(na.omit(curdat$godhatadj)/tau_d)
#   dmatseq[,j] <- seq(min(dvals[dvals>-Inf]),max(dvals[dvals>-Inf]),length.out=nseq)
#   for(i in 1:niter){
#     Gdpred[i,,j] <- Gcalc(zcons,dmatseq[,j],
#       go$alpha_G[i,j],go$beta_Gz[i,j],go$beta_Gd[i,j]
#       )
#     }
#   }
# 
# Gdquant <- apply(Gdpred,c(2,3),quantile,probs=probs)
# Gdquant <- aperm(Gdquant,c(2,1,3)) 
# # transpose to fit graph function

#####################
### RELATIONSHIPS ###
#####################

### GDD relationships

gdp <- list()
for(i in 1:nspecies){
  gdp[[i]] <- with(subset(msy90,species==spvals[i]&is.na(Ghat)==F&is.na(germdhat)==F&is.na(olsdhat)==F
    &(germdhat+olsdhat)>0),
    cor.test(Ghat,log(germdhat+olsdhat))
    )
  }
t(sapply(gdp,function(x) c(x$estimate,p=x$p.value,e=x$estimate<0&x$p.value<0.1)))

### PLOTS

Gpow <- 0.125 
mpow <- 0.175 # larger = more uneven sizes
Gq <- 0.9 
mq <- 0.4 # larger = smaller sizes

Gsize <- with(msy90,(totgerm+olsdhat*area_o)^Gpow/quantile(sqrt(totgerm+olsdhat*area_o)^Gpow,Gq))
msize <- with(msy90,
  ifelse(prevolsdhat==0,1,
    (prevnsdhat/(prevnsdhat+prevolsdhat))^mpow
    / quantile((prevnsdhat/(prevnsdhat+prevolsdhat))^mpow,mq,na.rm=T)
  )
)

with(msy90,
  preddatplot("S_relationships_max1_newonly_naspecies_tdistbpar_normndat_",
    xdat=log(prevnsdhat),
    ydat=log(Shat1),
    sind=species,
    xpred=lnsdseqs,
    ypred=log(Squant_new),
    xptype="matrix",
    pcol="purple",
    xname=expression(ln(N[n])),
    yname=expression(ln(S[3])),
    mycex=msize,
    ylim=c(-8,0)
  )
)

with(msy90,
  preddatplot("S_relationships_max1_oldonly_naspecies_tdistbpar_normndat_",
    xdat=log(prevolsdhat),
    ydat=log(Shat1),
    sind=species,
    xpred=lolsdseqs,
    ypred=log(Squant_old),
    xptype="matrix",
    pcol="purple",
    xname=expression(ln(N[n])),
    yname=expression(ln(S[T3])),
    mycex=msize,
    ylim=c(-8,0)
  )
)

with(msy90,
  preddatplot("G_z_relationships_olsdadj_naspecies_tdistbpar_normndat_",
    xdat=tprcp + log(tau_p),
    ydat=Ghat,
    sind=species,
    xpred=zseq + log(tau_p),
    ypred=Gzquantadj,
    xptype="vector",
    pcol="green",
    xname="w",
    yname="G",
    ylim=c(0,1),
    mycex=Gsize
  )
)

# with(msy90,
#   preddatplot("G_d_relationships_olsdadj_naspecies_tdistbpar_normndat_",
#     xdat=log(godhatadj),
#     ydat=Ghat,
#     sind=species,
#     xpred=dmatseq+log(tau_d),
#     ypred=Gdquant,
#     xptype="matrix",
#     pcol="green",
#     xname=expression(N[s]),
#     yname="G",
#     ylim=c(0,1),
#     mycex=Gsize
#   )
# )

#################
### RESIDUALS ###
#################

msy90$Ghat_fake <- msy90$Ghat
msy90$Ghat_fake[msy90$Ghat==1] <- max(msy90$Ghat[msy90$Ghat<1],na.rm=T)
msy90$Ghat_fake[msy90$Ghat==0] <- min(msy90$Ghat[msy90$Ghat>0],na.rm=T)

nyear <- nlevels(msy90$year)
coefd <- matrix(NA,nr=22,nc=2)
resid <- matrix(NA,nr=22,nc=nyear)
for(i in 1:nspecies){
	m <- glm(qlogis(Ghat_fake)~log(gprcp/tau_p),	
		data=subset(msy90,species==spvals[i] & is.na(Ghat_fake)==F)
		)
	coefd[i,] <- coef(m)
	resid[i,as.numeric(m$data$year)] <- resid(m)
	}
plot(coefd)
cor.test(coefd[,1],coefd[,2])

alpha_G_med <- apply(go$alpha_G,2,median)
beta_Gz_med <- apply(go$beta_Gz,2,median)
with(go,points(alpha_G_med,beta_Gz_med,col="red"))

par(mfrow=c(2,2),mar=c(4,4,1,1))
plot(coefd[,1]~alpha_G_med)
abline(0,1)
plot(coefd[,2]~beta_Gz_med)
abline(0,1)
plot(coefd[,1]~alpha_G_med)
abline(0,1)
text(alpha_G_med,coefd[,1],labels=spvals,pos=4)
plot(coefd[,2]~beta_Gz_med)
abline(0,1)
text(beta_Gz_med,coefd[,2],labels=spvals,pos=4)
  # separate glms overestimate, for some reason?

par(mfrow=c(1,1),las=2)
boxplot(resid,names=levels(msy90$year))
abline(h=0,col="red",lty=3)

################################
### CALIBRATION - PARAMETERS ###
################################

eps_G <- with(go,colMedians(eps_G))
eps_m <- with(go,colMedians(eps_m))
sig_G_eps <-  with(go,median(sig_G_eps))
sig_m_eps <-  with(go,median(sig_m_eps))

lgd <- go$lgd
lod <- go$lod

sig_g <- with(go,colMedians(sig_g))
sig_o <- with(go,colMedians(sig_o))
# table(sig_o/sig_g > 1)

nu_g <- with(go,colMedians(nu_g))
nu_o <- with(go,colMedians(nu_o))

par(mfrow=c(1,1))
plot(qlogis(msy90$Ghat)~qlogis(go$G))
abline(0,1,col="red")
table(msy90$Ghat>go$G)

with(msy90,plot(qlogis(tapply(Ghat,species,median,na.rm=T))~qlogis(tapply(go$G,species,median))))
abline(0,1,col="red")

par(mfrow=c(1,1))
plot(msy90$Ghat~go$G)

alphaG_med <- with(go,colMedians(alpha_G))
beta_Gz_med <- with(go,colMedians(beta_Gz))
# beta_Gd_med <- with(go,colMedians(beta_Gd))

alpha_m_med <- with(go,colMedians(alpha_m))
beta_m_med <- with(go,colMedians(beta_m))

G_med <- plogis(
  alphaG_med[dl$species] 
  + beta_Gz_med[dl$species]*dl$tprcp 
  # + beta_Gd_med[dl$species]*(log(msy90$germdhat+msy90$olsdhat)-log(tau_s))
  )

with(msy90,plot(qlogis(tapply(Ghat,species,median,na.rm=T))~qlogis(tapply(G_med,species,median))))
abline(0,1,col="red")

plot(qlogis(msy90$Ghat)~qlogis(G_med))
abline(0,1,col="red")
lines(supsmu(qlogis(G_med),qlogis(msy90$Ghat)),col="blue")
table(msy90$Ghat>G_med)/nrow(msy90)

pairs(cbind(alphaG_med,beta_Gz_med,beta_Gd_med,alpha_m_med,beta_m_med),
  lower.panel=panel.smooth, upper.panel=panel.cor)

# also - calibration of median G values

###################################
### CALIBRATION - DISTRIBUTIONS ###
###################################

tgreen <- rgb(red=0,green=1,blue=0,alpha=0.3,maxColorValue=1)
tblue <- rgb(red=0,green=0,blue=1,alpha=0.3,maxColorValue=1)

# make distribution of residuals too
# edit into histogram functions

# par(mfrow=c(1,2))
# hist(eps_G,breaks=50,prob=T)
# curve(dnorm(x,0,sig_G_eps),col="red",add=T)
# hist(eps_m,breaks=50,prob=T)
# curve(dnorm(x,0,sig_m_eps),col="red",add=T)

# pdf(paste0("obserror_hists_",format(Sys.Date(),"%d%b%Y"),".pdf"),
#   width=plotwidth,height=plotheight)
# par(mfrow=c(6,4))
# for(i in 1:nspecies){
#   cursp <- msy90$species==spvals[i]
#   mybreaks <- hist(eps_g[cursp],breaks=50,prob=T,main=spvals[i])
#   hist(rt(10^5,df=nu_g)*msy90$gsd[cursp][1],add=T,col=tgreen,prob=T,breaks=c(-Inf,mybreaks$breaks,Inf))
# }
# par(mfrow=c(6,4))
# for(i in 1:nspecies){
#   cursp <- msy90$species==spvals[i]
#   mybreaks <- hist(eps_o[cursp],breaks=50,prob=T,main=spvals[i])
#   hist(rt(10^5,df=nu_o)*msy90$osd[cursp][1],add=T,col=tblue,prob=T,breaks=c(-Inf,mybreaks$breaks,Inf))
# }
# dev.off()

npar <- nrow(msy90)
nsim <- 10^3
ntot <- npar*nsim
gsimmat <- osimmat <- area_g_mat <- area_o_mat <- matrix(NA,nr=npar,nc=nsim)
area_g_mat[] <- rep(msy90$area_g,nsim)
area_o_mat[] <- rep(msy90$area_o,nsim)
epsg <- rt(ntot,df=rep(median(go$nu_g),nsim))*rep(go$sig_g_vec,nsim)
epso <- rt(ntot,df=rep(median(go$nu_o),nsim))*rep(go$sig_o_vec,nsim)
  # is this wrong?

gsimmat[] <- rpois(ntot,
  lambda=exp(rep(go$lgd+log(area_g_mat[,1]*tau_s),times=nsim)+epsg)
) / area_g_mat
osimmat[] <- rpois(ntot,
  lambda=exp(rep(go$lod+log(area_o_mat[,1]*tau_s),times=nsim)+epso)
) / area_o_mat

pdf(paste0("hists_lnpois_BH_naspecies_tdistbpar_normndat_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=12,height=8)
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5))
mybreaks <- hist(msy90$germdhat,breaks=50,prob=T,density=20,xlab=expression(N[g]),main="")
hist(exp(go$lgd+log(tau_s)),add=T,col=tgreen,prob=T,breaks=c(0,mybreaks$breaks,Inf))
mybreaks <- hist(msy90$germdhat,breaks=50,prob=T,density=20,xlab=expression(N[g]),main="")
hist(gsimmat,add=T,col=tgreen,prob=T,breaks=c(0,mybreaks$breaks,Inf))
plot(log(msy90$germdhat)~I(go$lgd+log(tau_s)))
abline(0,1,col="red")

mybreaks <- hist(msy90$olsdhat,breaks=50,prob=T,density=20,xlab=expression(N[o]),main="")
hist(exp(go$lod+log(tau_s)),add=T,col=tblue,prob=T,breaks=c(0,mybreaks$breaks,Inf))
mybreaks <- hist(msy90$olsdhat,breaks=50,prob=T,density=20,xlab=expression(N[o]),main="")
hist(osimmat,add=T,col=tblue,prob=T,breaks=c(0,mybreaks$breaks,Inf))
plot(log(msy90$olsdhat)~I(go$lod+log(tau_s)))
abline(0,1,col="red")
dev.off()

comphist("germd_mean_hists_",x=as.numeric(msy90$species),y1=msy90$germdhat,y2=exp(go$lgd+log(tau_s)),xname=expression(N[o]))
comphist("olsd_mean_hists_",x=as.numeric(msy90$species),y1=msy90$olsdhat,y2=exp(go$lod+log(tau_s)),xname=expression(N[o]))

obsG <- gsimmat/(osimmat+gsimmat)
ng <- exp(rep(go$lgd+log(tau_s),times=nsim))
no <- exp(rep(go$lod+log(tau_s),times=nsim))
trueG <- ng/(no+ng)

yeah <- log(obsG/trueG)
hist(yeah,breaks=100)
table(yeah>0)
  # differences in distributions don't explain bias

library(moments)
with(msy90,cbind(
  tapply(germdhat,species,skewness)/
  tapply(olsdhat,species,skewness)
  ))

#################################
### CALIBRATION - PREDICTIONS ###
#################################

plot(log((msy90$germdhat+msy90$olsdhat+1)) ~ I(go$lbd+log(tau_s)))
abline(0,1,col="red")

par(mfrow=c(1,2))
plot(log(msy90$germdhat)~I(lgd+log(tau_s)))
abline(0,1,col="red")
plot(log(msy90$olsdhat)~I(lod+log(tau_s)))
abline(0,1,col="red")

x <- seq(-10,10,0.01)
ymed <- array(rep(x,22),dim=c(length(x),2,22))
nomat <- array(NA,dim=c(length(x),2,22))
library(abind)
ymat <- abind(nomat,ymed,nomat,along=2)

# Cal plots for each species
with(msy90,
  preddatplot("G_calibration_nozeroes_byspecies_naspecies_tdistbpar_normndat_",
    xdat=(go$lgd + log(tau_s)),
    ydat=log(msy90$germdhat0),
    sind=species,
    xpred=x,
    ypred=ymat,
    xptype="vector",
    pcol="black",
    xname=expression(N[g]~pred),
    yname=expression(N[g]~obs)
  )
)

with(msy90,
  preddatplot("G_Vs_O_byspecies_naspecies_tdistbpar_normndat_",
    xdat=log(msy90$olsdhat0),
    ydat=log(msy90$germdhat0),
    sind=species,
    xpred=x,
    ypred=ymat,
    xptype="vector",
    pcol="black",
    xname=expression(N[o]),
    yname=expression(N[g])
  )
)

with(msy90,
  preddatplot("O_calibration_nozeroes_byspecies_naspecies_tdistbpar_normndat_",
    xdat=(go$lod + log(tau_s)),
    ydat=log(msy90$olsdhat0),
    sind=species,
    xpred=x,
    ypred=ymat,
    xptype="vector",
    pcol="black",
    xname=expression(N[o]~pred),
    yname=expression(N[o]~obs)
  )
)

pdf(paste0("GO_calibration_nozeroes_naspecies_tdistbpar_normndat_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=8,height=4)
par(mfrow=c(1,2),mar=c(5,5,2,2),las=1,bty="l")

calplotall(
  xdat=(go$lgd + log(tau_s)),
  ydat=log(msy90$germdhat0),
  xlim=c(-6,8),
  ylim=c(-6,8),
  xlab=expression(Predicted~ln(N[g])),
  ylab=expression(Observed~ln(N[g]))
)
calplotall(
  xdat=(go$lod + log(tau_s)),
  ydat=log(msy90$olsdhat0),
  xlim=c(-3,11),
  ylim=c(-3,11),
  xlab=expression(Predicted~ln(N[o])),
  ylab=expression(Observed~ln(N[o]))
)

dev.off()

#####################################
### SPECIES-BY-YEAR OUTPUT MATRIX ###
#####################################

msy_molt <- melt(msy90,measure=c("totgerm","totlive"))
totgerm_mat <- acast(subset(msy_molt,variable=="totgerm"),species~year)
totlive_mat <- acast(subset(msy_molt,variable=="totlive"),species~year)

pgerm_mat <- totgerm_mat > 0
plive_mat <- totlive_mat > 0

ngerm <- colSums(pgerm_mat)
nlive <- colSums(plive_mat)

agerm <- with(subset(msy,species=="scba"), tapply(area_g,year,sum,na.rm=T))
alive <- with(subset(msy,species=="scba"), tapply(area_o,year,sum,na.rm=T))

pdf(paste0("Species_richness_timeseries_",format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=7,height=8)
par(mfrow=c(2,1),mar=c(6,6,2,2),las=1,bty="l")
plot(ngerm~as.numeric(names(ngerm)),type="b",pch=16,cex=sqrt(agerm)/max(sqrt(agerm)),
	xlab="Year",
	ylab="N species in germination surveys \n(max = 22)"
	)
plot(nlive~as.numeric(names(nlive)),type="b",pch=16,
	xlab="Year",
	ylab="N species in seed surveys \n(max = 22)"
	)
dev.off()

write.csv(totlive_mat[,(colnames(totlive_mat) %in% 1983:1989)==F],
	paste0("totlive_species_year_mat_",format(Sys.Date(),"%d%b%Y"),".csv")
	)

write.csv(totgerm_mat[,(colnames(totgerm_mat) %in% 1983:1989)==F],
	paste0("totgerm_species_year_mat_",format(Sys.Date(),"%d%b%Y"),".csv")
	)

#########################
### DI SURVIVAL TABLE ###
#########################

Stab <- adply(exp(-exp(go$alpha_m)),2,quantile,probs=probs)
names(Stab) <- c("Species",probs)
Stab$Species <- levels(msy$species)
Stab[,-1] <- round(Stab[,-1],2)
write.csv(Stab,file=paste0("DI_survival_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)

