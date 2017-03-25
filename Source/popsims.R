##########################################################
### Simulations of population dynamics for finite area ###
##########################################################

library(plyr)
library(reshape2)
library(parallel)
library(abind)
library(RColorBrewer)

### LOAD DATA 

source("Source/simulation_functions.R")
source("Source/figure_functions.R")
source("Source/prediction_functions.R")
source("Source/trait_functions.R")

msy <- read.csv("Output/msy_15Jan2016.csv",header=T)
Tvalues <- read.csv("Output/Tvalues_31Jul2015.csv",header=T)

# Define params -----------------------------------------------------------

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
nyear <- length(unique(msy$year))

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])
	# T in years

tau_s <- 100  # adujstment for seeds
tau_d <- 100	# adjustment for density
tau_p <- 100	# adjustment for rainfall

### CLIMATE

# DATA
pp <- read.csv("Output/prcp_projection_summaries_08Apr2016.csv",header=T)
mpam <- with(pp, median[measure=="mpam" & scenario==60 & yearcat==100])
mpcv <- with(pp, median[measure=="mpcv" & scenario==60 & yearcat==100])
  # using projected season precipitation for germination season precipitation change
  # (both very similar)
  # year = 2100
  # Representative Concentration Pathway 6.0

ncy <- read.csv("Output/ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
	# removes first value (missing because no previous winter)

zam <- mean(ncy$seasprcp)
zcv <- sd(ncy$seasprcp)/zam

wam <- mean(ncy$germprcp)
wcv <- sd(ncy$germprcp)/wam

### MODEL PARAMS (already permuted)

pl <- list(
	go = readRDS("Models/go_pars_tdistpois_naspecies_noerr_noGDD_loglik_BH_01Mar2017.rds"),
	gs = readRDS("Models/gnzhh_onhh_pars_medians_26Oct2015.rds"),
		# gs = g site level
		# source script: venable_Stan_GO_descriptive_gnzhh_onhh_26Oct2015
		# uses tau_s = 100
		# but tau actually irrelevant because all multiplicative?
	pr = readRDS("Models/pr_pars_yearhet_squared_pc_02Mar2016.rds"),
	rs = readRDS("Models/rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")
	)

# Sims --------------------------------------------------------------------

# 1000 per 0.1m^2 - *what does this mean?*

maml <- as.list(c(1,1,mpam,1,mpam))
mcvl <- as.list(c(0,1,1,mpcv,mpcv))

nclim <- length(maml)
cpc <- 2 # CORES per CLIMATE
ncores <- nclim*cpc
mpos <- rep(1:nclim,each=cpc)

nstart <- rep(10000,nspecies)
ni <- 500 # iterations PER CORE
nt <- 50
nj <- 22
nk <- 10000

# ni and nk must be >1
nit <- ni*cpc

set.seed(1)
maxiter <- 10000 # max number of itertions in PARAMETERISATION
cpos <- rep(1:cpc,times=nclim)
cipos <- rep(1:cpc,each=ni)
itersetl <- split(1:(ni*cpc),cipos)
  # requires that ni < maxiter

simp <- function(l){
  lapply(l,function(x){
    signif(x,2)
  })
}
  
cnames_unique <- gsub("\\.","",paste0("mu",simp(maml),"_cv",simp(mcvl)))
cnames_bycore <- paste0(rep(cnames_unique,each=cpc),"_s",rep(1:cpc,times=nclim))
cnames_merged <- paste(cnames_unique,collapse="_")

# !!! Change to PAIR climate replicates !!!

# Each core focuses on one climate
# Iterations for one climate can be split up over multiple cores
# (controlled by cpc)
system.time({
CL = makeCluster(ncores)
clusterExport(cl=CL, c("popsim","pl",
  "zam","zcv","wam","wcv",
  "mpos","maml","mcvl","cpos",
  "nstart","ni","nt","nj","nk",
  "Tvalues","itersetl","cnames_bycore"
  )) 
parLapply(CL, 1:ncores, function(n){
	mam <- maml[[mpos[n]]]
	mcv <- mcvl[[mpos[n]]]
	popsim(pl=pl,ni=ni,nt=nt,nj=nj,nk=nk,
		nstart=nstart,zam=zam*mam,zcv=zcv*mcv,
		wam=wam*mam,wcv=wcv*mcv,rho=0.82,
		Tvalues=Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
		iterset=itersetl[[cpos[n]]],
		savefile=paste0(cnames_bycore[n])
		)
	})
stopCluster(CL)
})

# Read back in ------------------------------------------------------------

### Small RAM read-in

# cnames_bycore_small <- cnames_bycore[grep("mu1_cv1",cnames_bycore)]
# ncores_small <- length(cnames_bycore_small)
# psl <- as.list(rep(NA,ncores_small))
# for(n in 1:ncores_small){
#   psl[[n]] <- readRDS(paste0(cnames_bycore_small[n],"_28Apr2016.rds"))
#   }
# names(psl) <- cnames_bycore_small

psl <- as.list(rep(NA,ncores))
for(n in 1:ncores){
  psl[[n]] <- readRDS(paste0("Sims/",cnames_bycore[n],"_17Mar2017.rds"))
  }
names(psl) <- cnames_bycore

varl <- psl[[1]]
nvar <- length(varl)
dimvar <- sapply(varl,dim)
ndimvar <- sapply(dimvar,length)
firstmpos <- match(1:nclim,mpos)

# Combine sim matrices by sim type
psla <- vector("list", nvar)
names(psla) <- names(varl)
for(i in 1:nvar){
  if(ndimvar[i]>0){
    for(j in 1:nclim){
      for(k in 1:cpc){
        if(k==1) ast <- psl[[firstmpos[j]]][[i]]
        else ast <- abind(ast,psl[[j+(k-1)]][[i]],along=1)
      }
      if(j==1){
        psla[[i]] <- array(ast,dim=c(dim(ast),1))
      }
      else{
        psla[[i]] <- abind(psla[[i]],ast,along=ndimvar[i]+1)
      }
    }
    dimnames(psla[[i]])[[ndimvar[i]+1]] <- cnames_unique
  }
}
psla <- psla[ndimvar>0]

# Derived parameters ------------------------------------------------------

### From input parameters

pl$pr$alpha_p <- pl$pr$beta_p[,,1]
pl$pr$beta_d_p <- pl$pr$beta_p[,,4]

pl$rs$beta_d_r <- pl$rs$beta_r[,,4]

pl$go$iota_mu <- with(pl$go, godmean_f(alpha_G,beta_Gz) )
pl$go$iota_sig <- with(pl$go, godvar_f(beta_Gz) )

pl$go$rho <- with(pl$go, alpha_G + beta_Gz*log(zam/tau_p))

pl$go$m0 <- exp(pl$go$alpha_m)
pl$go$m1 <- exp(pl$go$beta_m)
pl$go$Kn <- with(pl$go, Kncalc(m0,m1,T3))
pl$go$hn <- with(pl$go, hncalc(m0,m1,T3))

### From simulations

psla$Y <- with(psla,nn/ng)
psla$Ye <- with(psla,nnb/ng)
aNA <- array(dim=c(nit,1,nj,nclim))
psla$r <- log(psla$ns) - log(abind(aNA,psla$ns[,-nt,,],along=2))
# popualtion growth rates
psla$pY <- apply(psla$nn,2:4,function(x) sum(x>0)/length(x))
# probability of at least one new seed
# (could also use to calculate extinction risk)

# not surprising that differ in seed numbers (e.g. if make smaller seeds)
# calculate relative reproduction instead?

### From combination of input parameters and simulations

Knarr <- array(dim=c(nit,nj,nt)) # flip nj and nt later
Knarr[] <- pl$go$Kn[as.vector(unlist(itersetl)),]*(nk/10*tau_s) 
  # total K, adjusting for number of sites (nk)
  # fill in same for all nt
Knarr <- aperm(Knarr, c(1,3,2)) # flip nj and nt
psla$nsK <- psla$ns/rep(Knarr,nclim)

# to calculate: median G and So for timestep 0

rcatm <- cbind(rcaa[,,"ns"],mta)
rcata <- array(dim=c(dim(rcatm),1),
  dimnames=c(dimnames(rcatm),list("ns"))
)
# putting on extra dimension so that works with pairplot
rcata[,,1] <- unlist(rcatm)
# filling-in only works because rcatm has dim3 of 1

# Aggregate data ----------------------------------------------------------

seriesquant <- function(a,probs=c(0.05,0.50,0.95),keepdims=2:4){
	qarr <- apply(a,keepdims,quantile,prob=probs,na.rm=T)
	return(qarr)
}

q_ns <- seriesquant(log(psla$ns))
q_nn <- seriesquant(log(psla$nn))
q_G <- seriesquant(qlogis(psla$G))
q_Sn <- seriesquant(qlogis(psla$Sn))
  # Includes cases where no seeds were produced (DI survival)
q_Y <- seriesquant(log(psla$nn/psla$ng))
  # Ignoring cases where Ng=0  
q_Ye <- seriesquant(log(psla$Ye))
  # Ignoring cases where Ng=0  
q_nnb <- seriesquant(log(psla$nnb))
q_no <- seriesquant(log(psla$no))
q_ng <- seriesquant(log(psla$ng))
q_So <- seriesquant(qlogis(psla$So))
	# applies quantile function to list of densities for each climate scenario
	# dims: (quantile,time,species,clim)

# Plot results ------------------------------------------------------------

seriesplot(q_ns,"ns",yname=expression(ln(N[s])))
seriesplot(q_nn,"nn",yname=expression(ln(N[n])))
seriesplot(q_G,"G",yname="logit(G)")
seriesplot(q_Sn,"Sn",yname="logit(Sn)")
seriesplot(q_Y,"Y",yname="ln Y")
seriesplot(q_Ye,"Ye",yname="ln Yeff")
seriesplot(q_nnb,"nnb",yname=expression(ln~N[nb]))
seriesplot(q_no,"no",yname=expression(ln~N[o]))
seriesplot(apYf,"pY",yname="Pr(Y>0)")

# Relative change between scenarios ---------------------------------------

scenbase <- "mu1_cv1"
scennew <- c("mu081_cv1", "mu1_cv12", "mu081_cv12")
tpos <- 15 # averaging over range of z due to different iterations
keepsp <- (spvals %in% c("plpa","vuoc"))==F

qal <- list(ns=q_nsf,G=q_Gf,ng=q_ngf,Y=q_Yf,Ye=q_Yef,nn=q_nnf,Sn=q_Snf)

rcl <- lapply(qal,relchange,scenbase=scenbase,scennew=scennew,keepsp=keepsp)
rca <- array(dim=c(dim(rcl[[1]]),length(qal)),
  dimnames=c(dimnames(rcl[[1]]),list(names(qal)))
  )
rca[] <- unlist(rcl)

ccl <- lapply(qal,relchange,scenbase="mu1_cv0",scennew="mu1_cv1",keepsp=keepsp)
cca <- array(dim=c(dim(ccl[[1]]),length(qal)),
  dimnames=c(dimnames(ccl[[1]]),list(names(qal)))
  )
cca[] <- unlist(ccl)

### Quantifying mean-var interactions

mvint <- cca
mvint[] <- NA
dimnames(mvint)[[2]] <- "mvint"
mvint[] <- rca[,"mu081_cv12",] - (rca[,"mu081_cv1",] + rca[,"mu1_cv12",])
  # filling-in only works because mvint has 1 column
rcaa <- abind(rca,cca,mvint,along=2)

### P(nd>1)

pna <- abind(pYf,pYf,along=4)
pna <- aperm(pna,c(1,4,2,3))
  # hack to stack two copies of pYf along 2nd dimension so that relchange works

rpna <- relchange(qlogis(pna),scenbase="mu1_cv0",scennew="mu1_cv1",keepsp=keepsp)

### Pairs plots

pairplot("popchange_pairs_vitalratechange_newenv",rca,2)
pairplot("popchange_pairs_vitalratechange_consenv",cca,2)
pairplot("popchange_pairs_allclimscenarios_",rcaa,3)
pairplot("popchange_pairs_allclimscenarios_alltraits",rcata,3,w=21,h=21)

hist(mvint[,"ns"],breaks=30)
abline(v=0,col="red")

plot(cca[,1,"Sn"] ~ mta$lKnmed)
plot(cca[,1,"Sn"] ~ rpna)

### Parameters from different iterations

alpha_G <- go$alpha_G[iter_go,]
# calculate derived params, e.g. K, for each species

### Climate graphs

pdf(paste0("Plots/zdists",format(Sys.Date(),"%d%b%Y"),".pdf"),width=4.5,height=4.5)
plot(density(psls$mu1_cv1$z),xlim=c(-2,2),main="")
lines(density(psls$mu08_cv1$z),col="red")
lines(density(psls$mu1_cv12$z),col="blue")
lines(density(psls$mu08_cv12$z),col="purple")
abline(v=psls$mu1_cv0$z[1],col="green",lty=2)

plot(density(tau_p*exp(psls$mu1_cv1$z)),xlim=c(0,250),ylim=c(0,0.011),main="")
lines(density(tau_p*exp(psls$mu08_cv1$z)),col="red")
lines(density(tau_p*exp(psls$mu1_cv12$z)),col="blue")
lines(density(tau_p*exp(psls$mu08_cv12$z)),col="purple")
abline(v=tau_p*exp(psls$mu1_cv0$z[1]),col="green",lty=2)
dev.off()

# Optimal parameters within species ---------------------------------------

withinplot(pl$go,psls,"alpha_G","ns",simtrans_fun=log)
withinplot(pl$go,psls,"beta_Gz","ns",simtrans_fun=log)
withinplot(pl$go,psls,"alpha_m","ns",simtrans_fun=log)
withinplot(pl$go,psls,"beta_m","ns",simtrans_fun=log)

withinplot(pl$pr,psls,"alpha_p","ns",simtrans_fun=log)
withinplot(pl$pr,psls,"beta_d_p","ns",simtrans_fun=log)

withinplot(pl$go,psls,"iota_mu","ns",simtrans_fun=log)
withinplot(pl$go,psls,"iota_sig","ns",partrans_fun=log,simtrans_fun=log)
withinplot(pl$go,psls,"rho","ns",simtrans_fun=log)

withinplot(pl$go,psls,"alpha_G","nsK",simtrans_fun=log)
withinplot(pl$go,psls,"beta_Gz","nsK",simtrans_fun=log)

# Oslo talk plots ---------------------------------------------------------

tpos <- 15
myylim <- c(-5,0)

# MEANS AND VARIANCES

cclimpos <- which(cnames_unique=="mu1_cv1")
reld <- function(m){
  m[,2,tpos,] - rep(m[cclimpos,2,tpos,],each=dim(m)[1])
  }
md_nsf <- reld(q_nsf)
  # cclimpos = constant climate position (here mu1_cv1, later mu1_cv0)
  # 2 = median

keepnames <- c("mu08_cv1","mu1_cv12","mu08_cv12")
# keepnames <- c("mu09_cv1","mu1_cv11","mu09_cv11")
keepclim <- cnames_unique %in% keepnames
keepsp <- (spvals %in% c("plpa","vuoc"))==F
md_nsf2 <- md_nsf[keepclim,keepsp]
md_nsd <- melt(md_nsf2)
names(md_nsd) <- c("scenario","species","dlN")
  # dlN = difference relative to mu1_cv1

msy$olsdhat <- with(msy,totlive/area_o)
msy$godhat <- with(msy,germdhat+olsdhat)
G_med <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,median,na.rm=T)
  )
md_nsd$G_med <- G_med[md_nsd$species]

medvarbet <- function(G,...){
  median(G*(1-G),...)
  }
G_var <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,medvarbet,na.rm=T)
  )
md_nsd$G_var <- G_var[md_nsd$species]

myteal <- rgb(190,235,159,maxColorValue=255)

pdf(paste0("Output/mean_variance_comparison_",
  paste(keepnames,collapse="_"),
  "_t",tpos,"_",
  format(Sys.Date(),"%d%b%Y"),".pdf"
  ),
  width=4.5,height=4.5)

par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")
plot(dlN ~ scenario, 
  data=md_nsd,
  ylim=myylim,
  range=0,
  medlwd=2,
  boxwex=0.35,
  lty=1,
  ylab=expression(Delta~ln~population~size),
  xlab="Rainfall change scenario",
  names=c("mean","CV","both"),
  border="white" # makes invisible
  )

abline(h=0,lty=2)

plot(dlN ~ scenario, 
  data=md_nsd,
  ylim=myylim,
  range=0,
  col=myteal,
  medlwd=2,
  boxwex=0.35,
  lty=1,
  ylab=expression(Delta~ln~population~size),
  xlab="Rainfall change scenario",
  names=c("mean","CV","both")
  )

abline(h=0,lty=2)

# plot(dlN ~ qlogis(G_med),
#   data=subset(md_nsd,scenario==keepnames[3]),
#   ylim=myylim,
#   bg=myteal,
#   col="black",
#   pch=21,
#   cex=1.5,
#   ylab=expression(Delta~ln~population~size),
#   xlab="Species germination probability\n(logit scale)"
#   )
# 
# abline(h=0,lty=2)

dev.off()

# Other -------------------------------------------------------------------

cclimpos <- which(cnames_unique=="mu1_cv0")
  # changing reference climate to mu1_cv0
md_nsf3 <- reld(q_nsf)
md_nsd3 <- melt(md_nsf3[,keepsp])
names(md_nsd3) <- c("scenario","species","dlN")
G_mod <- plogis(q_Gf["mu1_cv1",2,tpos,])[keepsp]
consdiff <- md_nsd3$dlN[md_nsd3$scenario=="mu1_cv1"]

plot(dlN ~ qlogis(G_med[keepsp]), data=subset(md_nsd3,scenario=="mu1_cv1"))

iota_mu <- -go$alpha_G/go$beta_Gz
iota_mu_med <- apply(iota_mu,2,median)
iota_sig <- pi^2/(3*go$beta_Gz^2)
iota_sig_med <- apply(beta_G_var,2,median)
  # from Godfray & Rees 2002

pairs(cbind(G_mod,mu=iota_mu_med[keepsp],lsig=log(iota_sig_med[keepsp])))
  # lower G mostly reflects higher threshold
  # (i.e. not necessarily bet-hedging - take longer to germinate but do so at
  # same time)

plot(G_mod,consdiff)
plot(iota_mu_med[keepsp],consdiff)
plot(log(iota_sig_med[keepsp]),consdiff)
  # nothing predicts whether a bet-hedger (cur env definition) 

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
  # less var -> seem to do slightly better under mean change
  # note change from md_nsd3 to md_nsd

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
  # less var -> do better under variance change

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
  # less var -> do slightly better under overall change

matplot(qlogis(G_med),t(q_Gf[,2,tpos,]),pch=16)

diffplot <- function(x,y,xname,...){
  plot(y ~ x,bg=myteal,col="black",pch=21,ann=F,cex=1.5,...)
  addxlab(xname)
  }
# custlettlab <- function(i,label){
#   mtext(bquote(bold((.(letters[i])))~.(label)),
#     side=3,line=0.5,xpd=T,adj=0.05)
#   }
# addtoplab <- function(tname){
#   mtext(bquote(bold(tname)),side=3,line=2,las=1,xpd=T,cex=1.1)
#   }

diffplot_all <- function(...){
  xlabcur <- expression(Delta~ln~N~current~variability)
  par(mfrow=c(2,3),mar=c(5,2.5,5,2.5),oma=c(4,4,2,2),las=1,bty="l")
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu08_cv1"],xlabcur,...)
  mtext(bquote(bold(mean)),side=3,line=3,las=1,xpd=T,cex=1.3)
  addylab(expression(Delta~ln~population~size~(N)))
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu1_cv12"],xlabcur,...)
  mtext(bquote(bold(cv)),side=3,line=3,las=1,xpd=T,cex=1.3)
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],xlabcur,...)
  mtext(bquote(bold(both)),side=3,line=3,las=1,xpd=T,cex=1.3)
  
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu08_cv1"],"germination probability",...)
  addylab(expression(Delta~ln~population~size~(N)))
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu1_cv12"],"germination probability",...)
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],"germination probability",...)
  }

pdf(paste0("Plots/bet_hedging_prediction_",
    paste(keepnames,collapse="_"),
    "_t",tpos,"_",
    format(Sys.Date(),"%d%b%Y"),".pdf"
    ),
  width=10,height=8
  )
diffplot_all(type="n")
diffplot_all()
dev.off()

pdf(paste0("Plots/bet_hedging_prediction_K_So",
  "_t",tpos,"_",
  format(Sys.Date(),"%d%b%Y"),".pdf"
  ),width=9,height=4.5)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(4,4,1,1),las=1,bty="l")

diffplot(mta$lKnmed,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="ln max seed density",type="n")
addylab(expression(Delta~ln~population~size))
diffplot(mta$So,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="logit dormant seed survival",type="n")

diffplot(mta$lKnmed,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="ln max seed density")
addylab(expression(Delta~ln~population~size))
diffplot(mta$So,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="logit dormant seed survival")

dev.off()

summary(glm(md_nsd$dlN[md_nsd$scenario=="mu08_cv12"]~consdiff))
summary(glm(md_nsd$dlN[md_nsd$scenario=="mu08_cv12"]~G_mod))

plot(md_nsd$dlN[md_nsd$scenario=="mu1_cv12"] ~ md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])

# Bet-hedging definition plots --------------------------------------------

tpos <- 25
myylim <- c(-0.8,-0.2)

# MEANS AND VARIANCES

reld <- function(m){
  m[,2,tpos,] - rep(m[1,2,tpos,],each=dim(m)[1])
  }
md_nsf <- reld(q_nsf)

keepnames <- c("mu1_cv0","mu1_cv1")
# keepnames <- c("mu09_cv1","mu1_cv11","mu09_cv11")
keepclim <- cnames_unique %in% keepnames
md_nsf2 <- md_nsf[keepclim,]
md_nsd <- melt(md_nsf2)
names(md_nsd) <- c("scenario","species","dlN")

msy$olsdhat <- with(msy,totlive/area_o)
msy$godhat <- with(msy,germdhat+olsdhat)
G_med <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,median,na.rm=T)
  )
md_nsd$G_med <- G_med[md_nsd$species]

myteal <- rgb(190,235,159,maxColorValue=255)

csyp <- read.csv("Output/csyp_15Jan2016.csv",header=T)
cdpos <- read.csv("Output/cdpos_15Jan2016.csv",header=T)
prdat <- subset(csyp, ngerm>0 & is.na(nmeas)==F & is.na(germd)==F)
rsdat <- cdpos

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
yrvals <- unique(msy$year)
nyear <- length(yrvals)

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

xarr <- array(dim=c(nseq,4,nspecies))
for(i in 1:nspecies){
  subdat <- subset(prdat,species==i)
  xarr[,,i] <- c(
    rep(1,nseq),
    makeseq(log(subdat$prcp/(tau_p/2))), # change if change variable
    makeseq(log(subdat$prcp/(tau_p/2)))^2, # change if change variable
    makemean(log(subdat$germd/tau_d)) # change if change variable
    )
  }

# PR

pr_prcp_lp <- r_lp(pr$beta_p,ytype="pr",xtype="data",xarr=xarr)
pr_dens_lp <- r_lp(pr$beta_p,ytype="pr",xtype="data",xarr=xarr)

pr_prcp_quant <- plogis(r_quant(pr_prcp_lp))
pr_dens_quant <- plogis(r_quant(pr_dens_lp))

# RS

rs_prcp_lp <- r_lp(rs$beta_r,ytype="rs",xtype="data",xarr=xarr)
rs_dens_lp <- r_lp(rs$beta_r,ytype="rs",xtype="dens",xarr=xarr)

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

pcr_med <- apply(pcr_prcp_lp$pred_lp,c(2,3),median)
tprcp_med <- pcr_prcp_lp$xarr[,2,]
  # [,2,] -> select tprcp from x matrix

maxprcp <- lowprcp <- highprcp <- vector()
for(i in 1:nspecies){
  mp <- tprcp_med[pcr_med[,i]==max(pcr_med[,i]),i]
  is_lp <- which(tprcp_med[,i] < mp)
  is_hp <- which(tprcp_med[,i] > mp) 
  lp_0pos <- which( abs(pcr_med[is_lp,i]) == min(abs(pcr_med[is_lp,i])) )
  hp_0pos <- which( abs(pcr_med[is_hp,i]) == min(abs(pcr_med[is_hp,i])) )
  
  maxprcp[i] <- mp + log(tau_p)
  lowprcp[i] <- tprcp_med[is_lp[lp_0pos],i] + log(tau_p)
  highprcp[i] <- tprcp_med[is_hp[hp_0pos],i] + log(tau_p)
  }
 
rangeprcp <- highprcp - lowprcp

plottest <- function(x,...){
  plot(dlN ~ x,
    data=subset(md_nsd,scenario==keepnames[2]),
    bg=myteal,
    col="black",
    pch=21,
    cex=1.5,
    ylab=expression(Delta~ln~population~size),
    ...
    )
  mod <- glm(dlN ~ x,data=subset(md_nsd,scenario==keepnames[2]))
  summary(mod)
  }

library(matrixStats)
alpha_m <- colMedians(pl$go$alpha_m)
beta_m <- colMedians(pl$go$beta_m)
beta_Gz <- colMedians(pl$go$beta_Gz)

par(mfrow=c(2,3))
plottest(G_med,main="G_med")
plottest(beta_Gz,main="beta_Gz")
plottest(maxprcp,main="maxprcp")
plottest(rangeprcp,main="rangeprcp")
plottest(alpha_m,main="alpha_m")
plottest(beta_m,main="beta_m")

pairs(cbind(G_med,beta_Gz,maxprcp,rangeprcp,alpha_m,beta_m),
  lower.panel=panel.smooth, upper.panel=panel.cor)

summary( glm(dlN ~ maxprcp + rangeprcp,data=subset(md_nsd,scenario==keepnames[2])) )
resid <- resid(glm(dlN ~ rangeprcp,data=subset(md_nsd,scenario==keepnames[2])))
par(mfrow=c(1,1))
plot(maxprcp,resid)

# germination prob only makes a difference for 3 species; others opposite
# optimum environment makes no difference
# DI survival makes no difference
# strength of DD makes no difference

md_nsds <- subset(md_nsd,scenario=="mu1_cv1")
md_nsds$G_med <- G_med
md_nsds$beta_Gz <- beta_Gz
md_nsds$maxprcp <- maxprcp
md_nsds$rangeprcp <- rangeprcp
md_nsds$alpha_m <- alpha_m
md_nsds$beta_m <- beta_m

saveRDS(md_nsds,
  paste0(cnames_merged,"_keyparams_",format(Sys.Date(),"%d%b%Y"),".rds")
  )

### VENABLE 2007 ANALYSES

gmdat <- ddply(subset(msy,germdhat>0),.(species),summarise,
  gm = exp(sd(log((newseedhat+2)/(germdhat*area_n))))
  )
par(mfrow=c(1,1))
plot(gmdat$gm,G_med)
text(gmdat$gm,G_med,labels=gmdat$species,pos=4)

gmdat2 <- ddply(subset(msy,germdhat>0 & newseedhat>0),.(species),summarise,
  gm = exp(sd(log((newseedhat)/(germdhat*area_n))))
  )
par(mfrow=c(1,1))
plot(gmdat2$gm,G_med)
text(gmdat2$gm,G_med,labels=gmdat$species,pos=4)

# artefacts, e.g. low germ affects Y values that can be measured?

### NEW SEED SURVIVAL CHECKS

# looking at variance and skew in nn effects on Sn

library(fields)
library(sn)

BHSeg <- function(n,m0=exp(alpha_m[15]),m1=exp(beta_m[15]),T=T3){
  exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T))*n )	
  }

### VARIANCE

nmseq <- 100 
nnseq <- 1000
meanseq <- seq(9,10,length.out=nmseq)
sdseq <- seq(0,3,length.out=nmseq)
meanvec <- rep(meanseq,times=nmseq)
sdvec <- rep(sdseq,each=nmseq)

nnarr <- array(NA,dim=c(nnseq,nmseq,nmseq))
  # [entries,means,sds]
nnarr[] <- exp(rnorm(nnseq*nmseq*nmseq,rep(meanvec,each=nnseq),rep(sdvec,each=nnseq)))

smat <- qlogis(apply(nnarr,c(2,3),function(x) median(BHSeg(x/(500/10*tau_s)))))
image.plot(x=meanseq,y=sdseq,z=smat)

# variance alone can't explain why higher median seeds -> higher median survival
# (especially since variance seems to be bigger for red [variable] scenario)

scaleseq <- seq(1,2,length.out=nmseq)
shapeseq <- seq(1,2,length.out=nmseq)
scalevec <- rep(scaleseq,times=nmseq)
shapevec <- rep(shapeseq,each=nmseq)
nnarr <- array(NA,dim=c(nnseq,nmseq,nmseq))
  # [entries,means,sds]
nnarr[] <- exp(rsn(
  n=nnseq*nmseq*nmseq,
  xi=9,
  omega=rep(scalevec,each=nnseq),
  alpha=rep(shapevec,each=nnseq)
  ))

smat <- qlogis(apply(nnarr,c(2,3),function(x) median(BHSeg(x/(500/10*tau_s)))))
image.plot(x=scaleseq,y=shapeseq,z=smat)

plot(density(log(nnarr[,1,1])))
lines(density(log(nnarr[,1,50])),col="blue")
lines(density(log(nnarr[,1,100])),col="red")

# the more RIGHT(upward)-skewed the distribution, the lower the survival
# the more LEFT(downward)-skewed the distribution, the higher the survival
# (the latter inferred, not tested)
# (effect stronger at low skew than at high skew values)
# but in variable scenario, distribution is DOWNWARD-skewed, so should
# have even higher survival than constant scenario

plot(psls[[1]]$Sn[1:100,tpos,]~log(psls[[1]]$nn[1:100,tpos,]))
plot(qlogis(psls[[1]]$Sn[1:100,tpos,])~log(psls[[1]]$nn[1:100,tpos,]))

# Post-hoc analyses -------------------------------------------------------

pscur <- psls$mu1_cv1
ipos <- 1:1000
tpos <- 15

### r ~ z

x_in <- psls$mu1_cv1$z
y_in <- psls$mu1_cv1$r
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/popgrowth_rainfall_",format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=plotwidth,height=plotheight)
  
  plotsetup()
  
  for(j in 1:nspecies){

    y <- y_in[ipos,tpos,j]
    x <- x_in[ipos,tpos-1] # clim from previous "year"
    
    plot(y~x,pch=16,col=myblack)
    abline(h=0,col="red",lty=2)
    lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
    myrho <- round(cor.test(x,y)$estimate,2)
    legend("bottomright",bty="n",
      legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
      )
    
    lettlab(j)
    
    if(j %in% 19:22) addxlab("z") 
    if(j %in% seq(1,23,4)) addylab("r") 
    
    }
  
dev.off()
  
### r ~ G

x_in <- qlogis(psls$mu1_cv1$G)

pdf(paste0("Plots/popgrowth_G_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos-1,j] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab("logit(G)") 
  if(j %in% seq(1,23,4)) addylab("r") 
  
}

dev.off()


### r ~ ns

x_in <- log(psls$mu1_cv1$ns)
y_in <- psls$mu1_cv1$r
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/popgrowth_density_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos-1,j] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab(expression(ln(N[s])))
  if(j %in% seq(1,23,4)) addylab("r") 
  
}

dev.off()


pdf(paste0("Plots/popgrowth_hists_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  hist(y,breaks=1000,main="")
  abline(v=log(plogis(So_full[j])),col="red",lty=3)
  lettlab(j)

  if(j %in% 19:22) addxlab(expression(ln(N[s])))
  if(j %in% seq(1,23,4)) addylab("r") 
  
  }

dev.off()

### Y ~ z

x_in <- psls$mu1_cv1$z
y_in <- log(psls$mu1_cv1$nn/psls$mu1_cv1$ng)
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/Y_rainfall_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  abline(v=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab("z") 
  if(j %in% seq(1,23,4)) addylab("ln(Y)") 
  
}

dev.off()

y_in <- log(psls$mu1_cv1$nnb/psls$mu1_cv1$ng)
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/Ye_rainfall_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  abline(v=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab("z") 
  if(j %in% seq(1,23,4)) addylab(expression(ln(Y[e]))) 
  
}

dev.off()

### Y ~ gdens

x_in <-  log(psls$mu1_cv1$ng)
y_in <- log(psls$mu1_cv1$nn/psls$mu1_cv1$ng)
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/Y_gd_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos,j] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  abline(v=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab(expression(ln(N[g])))
  if(j %in% seq(1,23,4)) addylab("ln(Y)") 
  
  }

dev.off()


### Pr(Y>0) ~ gdens

### UNFINISHED!

x_in <-  log(psls$mu1_cv1$ng)
y_in <- apply(psls$mu1_cv1$nn[,nt,],2,function(x) sum(x>0)/length(x))
y_in[!is.finite(y_in)] <- NA

pdf(paste0("Plots/Y_gd_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  
  y <- y_in[ipos,tpos,j]
  x <- x_in[ipos,tpos,j] # clim from previous "year"
  
  plot(y~x,pch=16,col=myblack)
  abline(h=0,col="red",lty=2)
  abline(v=0,col="red",lty=2)
  lines(supsmu(x[is.na(y)==F],y[is.na(y)==F]),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab(expression(ln(N[g])))
  if(j %in% seq(1,23,4)) addylab("ln(Y)") 
  
  }

dev.off()

### Y~G

pdf(paste0("Plots/Y_vs_G_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  x <- with(pscur, qlogis(G[ipos,tpos,j]))
  y <- with(pscur, log(Y[ipos,tpos,j]))
  x <- x[y>-Inf]
  y <- y[y>-Inf]
  
  plot(y~x,pch=16,col=myblack)
  abline(lm(y~x),col="red",lty=2)
  lines(supsmu(x,y),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
    )
  
  lettlab(j)

  if(j %in% 19:22) addxlab("logit(G)") 
  if(j %in% seq(1,23,4)) addylab("ln Y") 
  
  }

plotsetup()

for(j in 1:nspecies){
  x <- with(pscur, qlogis(G[ipos,tpos,j]))
  y <- with(pscur, log(Y[ipos,tpos,j]))
  yb <- y > -Inf # y binary
  
  boxplot(x~yb,pch=16,col=myblack,horizontal=T,range=0,lty=1,boxwex=0.5)
  #m <- glm(yb~x,family="binomial")

  lettlab(j)
  
  if(j %in% 19:22) addxlab("logit(G)") 
  if(j %in% seq(1,23,4)) addylab("Y>0") 
  
  }

plotsetup()

for(j in 1:nspecies){
  x <- with(pscur, qlogis(G[ipos,tpos,j]))
  y <- with(pscur, log(Ye[ipos,tpos,j]))
  x <- x[y>-Inf]
  y <- y[y>-Inf]
  
  plot(y~x,pch=16,col=myblack)
  abline(lm(y~x),col="red",lty=2)
  lines(supsmu(x,y),col="blue")
  myrho <- round(cor.test(x,y)$estimate,2)
  legend("bottomright",bty="n",
    legend=substitute(paste(rho," = ",myrho,sep=""),list(myrho=myrho))
  )
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab("logit(G)") 
  if(j %in% seq(1,23,4)) addylab(expression(ln~Y[eff])) 
  
  }

plotsetup()

for(j in 1:nspecies){
  x <- with(pscur, qlogis(G[ipos,tpos,j]))
  y <- with(pscur, log(Ye[ipos,tpos,j]))
  yb <- y > -Inf # y binary
  
  boxplot(x~yb,pch=16,col=myblack,horizontal=T,range=0,lty=1,boxwex=0.5)
  #m <- glm(yb~x,family="binomial")
  
  lettlab(j)
  
  if(j %in% 19:22) addxlab("logit(G)") 
  if(j %in% seq(1,23,4)) addylab(expression(Y[eff]>0)) 
  
  }

dev.off()
