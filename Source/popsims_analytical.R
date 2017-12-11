### Simulations of population dynamics for infinite area ###
# Large edits from older version - grid ESS calculations no longer function

library(plyr)
library(reshape2)
library(parallel)
library(abind)
library(RColorBrewer)

### LOAD DATA 

source("Source/simulation_functions_analytical.R")
# source("Source/simulation_functions_stochastic.R")
source("Source/figure_functions.R")
source("Source/prediction_functions.R")
source("Source/trait_functions.R")
source("Source/invasion_functions.R")

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
pp <- read.csv("Output/prcp_projection_summaries_03Sep2017.csv",header=T)
mpam <- with(pp, median[measure=="mpam" & scenario==60 & yearcat==100])
mpsd <- with(pp, median[measure=="mpsd" & scenario==60 & yearcat==100])
  # using projected season precipitation for germination season precipitation change
  # (both very similar)
  # year = 2100
  # Representative Concentration Pathway 6.0

ncy <- read.csv("Output/ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
	# removes first value (missing because no previous winter)

zamo <- mean(log(ncy$seasprcp))
zsdo <- sd(log(ncy$seasprcp))

wamo <- mean(log(ncy$germprcp))
wsdo <- sd(log(ncy$germprcp))

### MODEL PARAMS (already permuted)

pl <- list(
  go = readRDS("Models/go_pars_lnGtnt_BH_25Nov2017.rds"),
	# go = readRDS("Models/go_pars_tdistpois_naspecies_noerr_noGDD_loglik_BH_01Mar2017.rds"),
	gs = readRDS("Models/gnzhh_onhh_pars_medians_26Oct2015.rds"),
		# gs = g site level
		# source script: venable_Stan_GO_descriptive_gnzhh_onhh_26Oct2015
		# uses tau_s = 100
		# but tau actually irrelevant because all multiplicative?
	pr = readRDS("Models/pr_pars_yearhet_squared_pc_02Mar2016.rds"),
	rs = readRDS("Models/rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")
	)

### MEDIAN PARAMETERS

# pls <- pl
# for(i in 1:length(pl)){
#   for(j in 1:length(pl[[i]])){
#     ndim <- length(dim(pl[[i]][[j]])) 
#     if(ndim==1){
#       pls[[i]][[j]] <- rep(median(pl[[i]][[j]]),length(pl[[i]][[j]]))
#     }
#     if(ndim>1){
#       keepdim <- 2:ndim # remove first dim, which is always iter
#       eachrep <- dim(pl[[i]][[j]])[1] # select iterdim
#       pls[[i]][[j]][] <- rep(apply(pl[[i]][[j]],keepdim,median),each=eachrep)
#     }
#     # if is.null(ndim), do nothing
#   }
# }

### PARAMS FOR SENSITIVITY ANALYSES

# transformed climate range approx. -1 -> 1
# tmrange <- c(-1.5,0.5)
# tsrange <- c(-3,1)
# plasticity <- T

# if(plasticity==F){
#   nsens <- 100
#   Gsens <- data.frame(
#     # alpha_G=qlogis(seq(0.001,0.999,length.out=nsens)),
#     alpha_G=seq(-5,5,length.out=nsens),
#     beta_Gz=rep(0,nsens)
#   )
# }
# 
# if(plasticity==T){
#   neach <- 10
#   tau_mu <- seq(tmrange[1],tmrange[2],length.out=neach)
#   tau_sd <- 2^seq(tsrange[1],tsrange[2],length.out=neach)
#   Gsens <- expand.grid(tau_mu=tau_mu,tau_sd=tau_sd)
#   nsens <- nrow(Gsens) # = neach^2
#   Gsens$alpha_G <- with(Gsens,godalpha_f(tau_mu,tau_sd))
#   Gsens$beta_Gz <- with(Gsens,godbeta_f(tau_sd))
# }

# if(plasticity==T){
#   neach <- 20
#   tau_scaled <- seq(tmrange[1],tmrange[2],length.out=neach)
#   tau_sd <- exp(seq(tsrange[1],tsrange[2],length.out=neach))
#   Gsens <- expand.grid(tau_scaled=tau_scaled,tau_sd=tau_sd)
#   Gsens$tau_mu <- with(Gsens, tau_scaled * tau_sd)
#   nsens <- nrow(Gsens) # = neach^2
#   Gsens$alpha_G <- with(Gsens,godalpha_f(tau_mu,tau_sd))
#   Gsens$beta_Gz <- with(Gsens,godbeta_f(tau_sd))
# }

# if(plasticity==T){
#   neach <- 10
#   tau_mu <- seq(-2,2,length.out=neach)
#   beta_Gz <- seq(0.1,3,length.out=neach)
#   Gsens <- expand.grid(tau_mu=tau_mu,beta_Gz=beta_Gz)
#   Gsens$tau_sd <- with(Gsens,sqrt(pi^2/(3*beta_Gz^2)))
#   Gsens$alpha_G <- with(Gsens,godalpha_f(tau_mu,tau_sd))
#   nsens <- nrow(Gsens) # = neach^2
# }

# if(plasticity==T){
#   neach <- 10
#   alpha_G <- seq(-2,2,length.out=neach)
#   tau_sd <- seq(0.1,3,length.out=neach)
#   Gsens <- expand.grid(alpha_G=alpha_G,tau_sd=tau_sd)
#   Gsens$beta_Gz <- with(Gsens,sqrt(pi^2/(3*tau_sd^2)))
#   Gsens$tau_mu <- with(Gsens,-alpha_G/beta_Gz)
#   nsens <- nrow(Gsens) # = neach^2
# }

ESG <- readRDS("Sims/ESSall_pTRUE_nk0_10Dec2017.rds")

niE <- ESG$ni
nclimE <-  3 # number of ESS climates
rpi <- 10 # number of replicated climate runs for each parameter set

ESG$alpha_G <- with(ESG$psla, abind(am[,,1],am[,,2],am[,,3],along=1))
ESG$beta_Gz <- with(ESG$psla, abind(bm[,,1],bm[,,2],bm[,,3],along=1))
nit <- rpi * niE * nclimE
iseq <- rep(1:niE, each=rpi*nclimE)
iseqG <- rep(1:(niE*nclimE), each=rpi)
  # repetitions vary fastest, ESS climate varies slowest
iseqz <- rep(1:(rpi*niE),times=nclimE) 

pls <- pl
pls$go$alpha_G <- ESG$alpha_G[iseqG,]
pls$go$beta_Gz <- ESG$beta_Gz[iseqG,]

# nplot <- 4
# Gshow <- round(
#   rep(seq(1,neach,length.out=nplot),times=nplot)
#   + rep(seq(0,neach^2-neach,length.out=nplot),each=nplot),
#   0)
# Gplot <- Gsens[Gshow,]
# 
# pdf(paste0("Plots/Germination_functions_", format(Sys.Date(),"%d%b%Y"),".pdf"),
# 		width=7,height=7)
# par(mfrow=c(nplot,nplot),mar=c(3,3,1,1),las=1,ann=F,bty="l")
# for(i in 1:nrow(Gplot)){
#   with(Gplot, curve(plogis(alpha_G[i]+beta_Gz[i]*x),
#     xlim=c(-2,2),ylim=c(0,1)) )
#   if(plasticity==T){
#     legend("topleft",bty="n",
#       legend=c(paste0("mu=",signif(Gplot$tau_mu[i],2)),
#         paste0("sd=",signif(Gplot$tau_sd[i],2))
#         )
#       )
#   }
# }
# dev.off()

# Sims --------------------------------------------------------------------

# Each core focuses on one climate
# Iterations for one climate can be split up over multiple cores
# (controlled by cpc)

# maml <- as.list(c(1,1,mpam,1,mpam,mpam))
# msdl <- as.list(c(0,1,1,mpsd,mpsd,0))
maml <- as.list(c(0,mpam)) 
msdl <- as.list(c(1,mpsd)) 
  # scaling mean log rainfall (zamo) only works because sign stays the same

nclim <- length(maml)
cpc <- 30 # CORES per CLIMATE (assumed equal for resident and invader)
ncores <- nclim*cpc
mpos <- rep(1:nclim,each=cpc)

nstart <- rep(1,nspecies)
nt <- 1025
nj <- 22
  # min invader iterations per core = rpi * nsens

cpos <- rep(1:cpc,times=nclim)
cipos <- rep(1:cpc,each=nit/cpc)

itersetl <- split(iseq,cipos)
itersetlG <- split(iseqG,cipos)
itersetlz <- split(iseqz,cipos)

ni <- length(itersetl[[1]]) # iterations PER CORE for RESIDENT simulations
# nii <- length(itersetli[[1]]) # iterations per core for INVADER simulations

simp <- function(l){
  lapply(l,function(x){
    signif(x,2)
  })
}
 
cnames_unique <- gsub("\\.","",paste0("mu",simp(maml),"_sd",simp(msdl)))
cnames_bycore <- paste0(rep(cnames_unique,each=cpc),"_s",rep(1:cpc,times=nclim))
cnames_merged <- paste(cnames_unique,collapse="_")
  # same for both invaders and residents

# maxiter <- 10000 # max number of iterations in PARAMETERISATION

# Simulate normalised climate and year effects ----------------------------

require(MASS)
set.seed(1)
zwn_mu <- rep(0,2)
zwn_sig <- matrix(c(1,rep(0.82,2),1),nr=2,nc=2) # correlation = 0.82

zwyols <- list()
for(i in 1:(rpi*niE)){
  zwn <- mvrnorm(n=nt, mu=zwn_mu, Sigma=zwn_sig)
  zwyols[[i]] <- data.frame(
    zn = zwn[,1],
    wn = zwn[,2],
    eps_y_pn = rnorm(nt,0,1),
    eps_y_rn = rnorm(nt,0,1)
  )
}
zwyol <- rep(zwyols,nclimE)
  # different climate / reproduction sequence for each rep and ESS starting pars
  # same sequences for ESSs from different climates
  
zwyoli <- list()
for(i in 1:cpc){
  zwyoli[[i]] <- do.call("rbind", zwyol[itersetlz[[i]]])
}

# CALCULATE GERMINATION BEFOREHAND
# REDUCE S0 TO NI*NJ VALUES
# CYCLE ITERATIONS WHEN >1000

# Resident simulations ----------------------------------------------------

set.seed(1)
system.time({
CL = makeCluster(ncores)
clusterExport(cl=CL, c("popana","pls", 
  "ni","nt","nj","nstart",
  "zwyoli","zamo","zsdo","wamo","wsdo",
  "mpos","maml","msdl","cpos",
  "Tvalues","cnames_bycore",
  "itersetl","itersetlG"
  )) 
parLapply(CL, 1:ncores, function(n){
	mam <- maml[[mpos[n]]]
	msd <- msdl[[mpos[n]]]
	popana(pl=pls,ni=ni,nt=nt,nj=nj,
		nstart=nstart,
	  zam=zamo+mam*zsdo,zsd=zsdo*msd,
		zwy=zwyoli[[cpos[n]]],
	  wam=wamo+mam*wsdo,wsd=wsdo*msd,
		Tvalues=Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
		iterset=itersetl[[cpos[n]]],itersetG=itersetlG[[cpos[n]]],
		savefile=paste0("res_",cnames_bycore[n]) # res -> residents
		)
	})
stopCluster(CL)
})
  # 10 mins: 30 cores, nit=3000, nt=1025

# Read resident simulations back in ---------------------------------------

### Small RAM read-in

# cnames_bycore_small <- cnames_bycore[grep("mu1_sd1",cnames_bycore)]
# ncores_small <- length(cnames_bycore_small)
# psl <- as.list(rep(NA,ncores_small))
# for(n in 1:ncores_small){
#   psl[[n]] <- readRDS(paste0(cnames_bycore_small[n],"_28Apr2016.rds"))
#   }
# names(psl) <- cnames_bycore_small

psl <- as.list(rep(NA,ncores))

dir <- paste0(getwd(),"/Sims/")
files <- paste0(dir,list.files(dir))

for(n in 1:ncores){
  curname <- paste0("Sims/res_",cnames_bycore[n],"_10Dec2017.rds")
  finished <- grep(curname,files)
  if(length(finished)!=0){
    psl[[n]] <- readRDS(curname)
  }
}
names(psl) <- cnames_bycore

# psl[is.na(psl)] <- psl[1]
psla <- simcombine(psl,nclim=nclim,cpc=cpc)

# Invader simulations -----------------------------------------------------

# Clunky input because saves on total RAM used

nchunk <- 10
lchunk <- ncores/nchunk
corechunk <- split(1:ncores,rep(1:nchunk,each=lchunk))
for(c in 1:nchunk){ # 1:length(corechunk)
  curchunk <- corechunk[[c]]
  system.time({
    CL = makeCluster(ncores)
    clusterExport(cl=CL, c(
      "popinv","mpos","cpos",
      "pls", "psla",
      "itersetli","itersetlr",
      "nii","nt","nj","nstart",
      "cnames_bycore",
      "tmin",
      "curchunk"
    )) 
    parLapply(CL, curchunk, function(n){ # doing piece-by-piece to save RAM
      climpos <- mpos[n]
      iterseti <- itersetli[[cpos[n]]]
      itersetr <- itersetlr[[cpos[n]]]
      popinv(  
        z=psla$z[itersetr,,climpos],
        w=psla$w[itersetr,,climpos],
        alpha_G=pls$go$alpha_G[iterseti,], # change list names if needed
        beta_Gz=pls$go$beta_Gz[iterseti,], 
        Y=psla$Y[itersetr,,,climpos],
        Sn=psla$Sn[itersetr,,,climpos],
        So=psla$So[itersetr,,,climpos],
        ni=nii,nt=nt,nj=nj,      
        tmin=tmin,
        full=F,
        savefile=paste0("inv_",cnames_bycore[n]) # inv -> invaders
      )
    })
    stopCluster(CL)
  })
  # 18 mins
}
  
# Read invader simulations back in ----------------------------------------

psl2 <- as.list(rep(NA,ncores))
for(n in 1:ncores){
  psl2[[n]] <- readRDS(paste0("Sims/inv_",cnames_bycore[n],"_23Aug2017.rds"))
}
names(psl2) <- cnames_bycore

psla2 <- simcombine(psl2)

# Time series graphs ------------------------------------------------------

j <- 19
require(fields)
ex <- psla$ns[,nt,j,1] < 0.001  
  # remove automatically extinct resident / invader
matplot(t(log(psla$ns[!ex,,j,1])),type="l",lty=1,col="black")
matplot(t(log(psla2$ns[rep(!ex,each=nsens),,j,1])),type="l",
  lty=rep(1:nsens,each=nsens),
  col=rep(tim.colors(nsens)[!ex],each=nsens),
  add=T
  )
 # cols = RESIDENT population
 # lines = different strategies


# PIPs --------------------------------------------------------------------

psla2$rbari <- with(psla2, apply(ri[,tmin:nt,,],c(1,3,4),mean))

pinvf <- function(x){
  mean(x>0)
}
  # fraction to have grown by end of interval
  
pip <- array(dim=c(nsens,nsens,nj,nclim))
for(m in 1:nclim){
  for(j in 1:nj){
    pip[,,j,m] <- tapply(
      psla2$rbari[,j,m],
      list(psla$G[iseqres,nt,j,m],psla2$G[,nt,j,m]),
      pinvf
      )
  }
}

# pip[ex,] <- NA
# pip[,ex] <- NA
  # remove automatically extinct resident / invader

Gmed <- with(msy,tapply(qlogis(germdhat/(olsdbar+germdhat)),species,median,na.rm=T))
Gest <- apply(pl$go$alpha_G,2,median)
pipplot(z=pip,xname=expression(G[r]),yname=expression(G[i]),
  pointvals=list(Gmed,Gest),
  x=alphaGseq,y=alphaGseq
  )

alphaGseq <- Gsens$alpha_G
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,19,1])
abline(0,1)
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,19,2])
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,19,3])
  # resident on x, invader on y
  # patterns at high G for constant environment unclear because 
  # take longer to displace?
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,15,1])
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,15,2])
image.plot(x=alphaGseq,y=alphaGseq,z=pip[,,15,3])

# Explaining PIP patterns -------------------------------------------------

psla$Y <- with(psla,nn/ng)
psla$Ye <- with(psla,nnb/ng)
hist(log(psla$Y[,nt,19,1]),breaks=100)
hist(log(psla$Ye[psla$G[,nt,19,1]>0.5,nt,19,1]),breaks=1000)

# PIPs - multiple variables -----------------------------------------------

pinvf <- function(x){
  mean(x>0)
}

tau_mu_res <- rep(Gsens$tau_mu,each=rpi)[iseqres]
tau_mu_inv <- rep(Gsens$tau_mu,each=rpi)[iseqinv]
tau_sd_res <- rep(Gsens$tau_sd,each=rpi)[iseqres]
tau_sd_inv <- rep(Gsens$tau_sd,each=rpi)[iseqinv]

intr <- rep(1:nsens,each=rpi)[iseqres]
inti <- rep(1:nsens,each=rpi)[iseqinv]

pip <- array(dim=c(nsens,nsens,nj,nclim))
for(m in 1:nclim){
  for(j in 1:nj){
    pip[,,j,m] <- tapply(
      psla2$rbari[,j,m],
      list(intr,inti),
      pinvf
    )
  }
}

pipplot(z=pip,xname=expression(G[r]),yname=expression(G[i]))

pipmin <- apply(pip,2:4,min,na.rm=T)
pipopt <- apply(pipmin,2:3,function(x) which(x==max(x)))  
  # taking only first optimum!

# Optimal species parameters ----------------------------------------------

require(fields)
tmin <- 15

popmed <- apply(log(psla$ns[,tmin:nt,,]),c(1,3,4),median)
pop <- array(dim=c(neach,neach,nj,nclim))
for(m in 1:nclim){
  for(j in 1:nj){
    pop[,,j,m] <- tapply(
      popmed[,j,m],
      as.list(Gsens[rep(1:nsens,each=rpi),c("tau_mu","tau_sd")]),
      median
    )
  }
}

goodpop <- pop
goodpop[goodpop < -1] <- NA
image.plot(x=tau_mu,y=log(tau_sd),z=goodpop[,,j,1])
image.plot(x=tau_mu,y=log(tau_sd),z=goodpop[,,j,2])

popopt <- apply(pop,3:4,function(x) which(x==max(x)))  

# Compare individual and species optima -----------------------------------

j <- 13

for(m in 1:nclim){
  
  curpipopt <- unlist(pipopt[j,m])
  curpopopt <- unlist(popopt[j,m])
  
  for(i in 1:length(curpipopt)){
    with(Gsens[unlist(curpipopt[i]),],{
      if(m==1) curve(plogis(alpha_G + beta_Gz*x),col=m,xlim=c(-2,0),ylim=c(0,1))
      if(m!=1) curve(plogis(alpha_G + beta_Gz*x),col=m,add=T)
    })
  }
  
  for(i in 1:length(curpopopt)){
    with(Gsens, 
      curve(plogis(alpha_G[curpopopt[i]] + beta_Gz[curpopopt[i]]*x),n=10^4,add=T,lty=2,col=m)
      )
  }
  
}
  # small differences in tau_mu only matter when lots of years to distinguish them



# OLD STUFF #


# Separate-out different ES G iterations ----------------------------------

barquant <- function(a,probs=c(0.25,0.50,0.75),keepdims=3:4){
  qarr <- apply(a,keepdims,quantile,prob=probs,na.rm=T)
  return(qarr)
}
# 50% quantiles, not 90% quantiles!

Epos <- rep(1:nclimE,each=rpi*niE)
lns <- array(dim=c(nit,nt,nj,nclim,nclimE))
for(i in 1:nclimE){
  lns[,,,,i] <- log(psla$ns[Epos==i,,,])
}

dlns <- lns[,,,2,] - lns[,,,1,]

q_lns <- barquant(lns[,,,1,])
q_dlns <- barquant(dlns)
qlist <- list(lns=q_lns,dlns=q_dlns)

saveRDS(qlist,paste0("Sims/ESS_medpops_",format(Sys.Date(),"%d%b%Y"),".rds"))
qlist <- readRDS(paste0("Sims/ESS_medpops_",format(Sys.Date(),"%d%b%Y"),".rds"))

Escen <- c("low","medium","high")
myteal <- rgb(190,235,159,maxColorValue=255)

purples <- brewer.pal(9,"Purples")[5] 
blues <- brewer.pal(9,"Blues")[5] 
greens <- brewer.pal(9,"Greens")[5] 
oranges <- brewer.pal(9,"Oranges")[5]
reds <- brewer.pal(9,"Reds")[5] 
greys <- brewer.pal(9,"Greys")[5] 

cols <- c(purples,blues,greens,oranges,reds,greys)

colpos <- 1:nclimE

tran <- 25
cols_rgb <- col2rgb(cols)
trancols <- rgb(
  red=cols_rgb[1,],
  green=cols_rgb[2,],
  blue=cols_rgb[3,],
  alpha=tran,
  maxColorValue = 255
)

# pdf(paste0("mean_variance_comparison_",
#            paste(keepnames,collapse="_"),
#            "_t",tpos,"_",
#            format(Sys.Date(),"%d%b%Y"),".pdf"
# ),
# width=4.5,height=4.5)

par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")
boxplot(qlist$lns[2,,], 
     #ylim=myylim,
     range=0,
     medlwd=2,
     boxwex=0.35,
     lty=1,
     ylab=expression(Delta~ln~population~size),
     xlab="Rainfall change scenario",
     names=Escen,
     border="white" # makes invisible
)

abline(h=0,lty=2)

boxplot(qlist$lns[2,,], 
        #ylim=myylim,
        range=0,
        medlwd=2,
        boxwex=0.35,
        lty=1,
        ylab=expression(ln~population~size),
        xlab="Rainfall change scenario",
        names=Escen,
        col=cols[colpos]
        )

par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")
boxplot(qlist$dlns[2,,], 
        #ylim=myylim,
        range=0,
        medlwd=2,
        boxwex=0.35,
        lty=1,
        ylab=expression(ln~population~size),
        xlab="Rainfall change scenario",
        names=Escen,
        border="white" # makes invisible
)

abline(h=0,lty=2)

boxplot(qlist$dlns[2,,], 
        #ylim=myylim,
        range=0,
        medlwd=2,
        boxwex=0.35,
        lty=1,
        ylab=expression(Delta~ln~population~size),
        xlab="Rainfall change scenario",
        names=Escen,
        col=cols[colpos]
)

# abline(h=0,lty=2)


# Climate distributions ---------------------------------------------------
# zam=zamo+mam*zsdo,zsd=zsdo*msd

climcurve <- function(mam,msd){
  curve(dnorm(x,zamo+zsdo*mam,zsdo*msd),xlim=c(zamo-2.5*zsdo,zamo+2.5*zsdo),ylim=c(0,0.85))
}
climcurve(0,1/msdl[[2]])
climcurve(0,1)
climcurve(0,msdl[[2]])
climcurve(maml[[2]],msdl[[2]])


# Extract parameters used in simulations ----------------------------------

iterextract <- function(p){
  pdim <- dim(p)
  if(maxiter %in% pdim){
    if(length(pdim)==1) return(p[unlist(itersetl)])
    if(length(pdim)==2) return(p[unlist(itersetl),])
    if(length(pdim)==3) return(p[unlist(itersetl),,])
  }
}

goi <- lapply(pl$go,iterextract)
pri <- lapply(pl$pr,iterextract)
rsi <- lapply(pl$rs,iterextract)

gois <- lapply(pls$go,iterextract)

# Derived parameters ------------------------------------------------------

### From input parameters

goi$tau_mu <- with(goi, godmean_f(alpha_G,beta_Gz) )
goi$tau_sig <- with(goi, godvar_f(beta_Gz) )
goi$rho <- with(goi, alpha_G + beta_Gz*log(zamo/tau_p))

gois$tau_mu <- with(gois, godmean_f(alpha_G,beta_Gz) )
gois$tau_sig <- with(gois, godvar_f(beta_Gz) )
gois$rho <- with(gois, alpha_G + beta_Gz*log(zamo/tau_p))

goi$m0 <- exp(goi$alpha_m)
goi$m1 <- exp(goi$beta_m)
goi$So3 <- exp(-goi$m0*T3)
goi$Kn <- with(goi, Kncalc(m0,m1,T3))
goi$hn <- with(goi, hncalc(m0,m1,T3))
goi$taun <- with(goi, tauncalc(m0,m1,T3))

beta_b <- pri$beta_p + rsi$beta_r
  # = log(pr/(1-p))

nseq <- 1000
lNgseq <- seq(-100,100,length.out=nseq) # true values: +log(tau_s)
modmat <- matrix(c(rep(1,(nseq*3)),lNgseq),nr=nseq,nc=4)

ylNg <- array(dim=c(nit,nj,nseq))
for(i in 1:nit){
  for(j in 1:nj){
    ylNg[i,j,] <- log(
      plogis(
        pri$beta_p[i,j,] %*% t(modmat)) 
        * nbtmean(exp(rsi$beta_r[i,j,] %*% t(modmat)),rsi$phi_r[i])
        # need to adjust for truncation
      )
  }
}

lKy <- apply(ylNg,c(1,2),function(x){
  xabs <- abs(x-0)
  lNgseq[which(xabs==min(xabs))] 
  }
)

### From simulations

psla$Y <- with(psla,nn/ng)
psla$Ye <- with(psla,nnb/ng)
aNA <- array(dim=c(nit,1,nj,nclim))
psla$r <- log(abind(psla$ns[,-1,,],aNA,along=2)) - log(psla$ns)
  # popualtion growth rates (year t = growth to next year
psla$pY <- apply(psla$nn,2:4,function(x) sum(x>0)/length(x))
  # probability of at least one new seed
  # (could also use to calculate extinction risk)
psla$pYr <- apply(psla$Y,2:4,function(x){
  xnew <- x[!is.nan(x)]
  sum(xnew>=1)/length(xnew)
})
  # self-replacement rate (before seed DD)
psla$lYS <- with(psla,log(Ye)-log(So))
psla$mnmo <- with(psla,log(Sn) - log(So)*T3) 
  # log(S) = -m*T
  # ln(Sn/So) = ln(Sn) - ln(So) = -mn + mo
  # Convert to T3: ln(Sn) - ln(So)*T3
psla$pP <- apply(psla$ns,2:4,function(x) sum(x>0)/length(x))
  # Probability of persistence

psla$ns_p <- psla$ns_e <- psla$ns
for(m in 1:nclim){
  for(j in 1:nspecies){
    psla$ns_p[psla$ns[,50,j,m]==0,,j,m] <- NA # extinct -> NA
    psla$ns_e[psla$ns[,50,j,m]>0,,j,m] <- NA  # alive -> NA
  }
}
  # population density given that survived to t=50

# not surprising that differ in seed numbers (e.g. if make smaller seeds)
# calculate relative reproduction instead?

### From combination of input parameters and simulations

Knarr <- hnarr <- array(dim=c(nit,nj,nt)) # flip nj and nt later
Knarr[] <- goi$Kn*nk/10*tau_s
hnarr[] <- goi$hn*nk/10*tau_s
  # total K: (K per 0.01m^2) * 10 * (number of 0.1m^2 plots)
  # fill in same for all nt
Knarr <- aperm(Knarr, c(1,3,2)) # flip nj and nt
hnarr <- aperm(hnarr, c(1,3,2)) # flip nj and nt

psla$nsK <- psla$ns/rep(Knarr,nclim)
psla$nnK <- psla$nn/rep(Knarr,nclim)

# to calculate: median G and So for timestep 0

rcatm <- cbind(rcaa[,,"ns"],mta)
rcata <- array(dim=c(dim(rcatm),1),
  dimnames=c(dimnames(rcatm),list("ns"))
)
# putting on extra dimension so that works with pairplot
rcata[,,1] <- unlist(rcatm)
# filling-in only works because rcatm has dim3 of 1

# set reference climate
# for each climate, plot difference against parameter for each species

# Aggregate data ----------------------------------------------------------

seriesquant <- function(a,probs=c(0.25,0.50,0.75),keepdims=2:4){
	qarr <- apply(a,keepdims,quantile,prob=probs,na.rm=T)
	return(qarr)
}
  # 50% quantiles, not 90% quantiles!

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
q_lYS <- seriesquant(psla$lYS)
q_nnK <- seriesquant(log(psla$nnK))
q_mnmo <- seriesquant(psla$mnmo)

# Plot results ------------------------------------------------------------

purples <- brewer.pal(9,"Purples")[5] 
blues <- brewer.pal(9,"Blues")[5] 
greens <- brewer.pal(9,"Greens")[5] 
oranges <- brewer.pal(9,"Oranges")[5]
reds <- brewer.pal(9,"Reds")[5] 
greys <- brewer.pal(9,"Greys")[5] 

cols <- c(purples,blues,greens,oranges,reds,greys)

colledgetext <- cnames_unique
detledgetext <- c(
  paste0("nstart=",nstart[1]),
  paste0("ni=",ni),
  paste0("nt=",nt),
  paste0("nj=",nj)
)

seriesplot(q_ns,"ns",yname=expression(ln(N[s])))
seriesplot(q_nn,"nn",yname=expression(ln(N[n])))
seriesplot(q_G,"G",yname="logit(G)")
seriesplot(q_Sn,"Sn",yname="logit(Sn)")
seriesplot(q_So,"So",yname="logit(So)")
seriesplot(q_Y,"Y",yname="ln Y")
seriesplot(q_Ye,"Ye",yname="ln Yeff")
seriesplot(q_nnb,"nnb",yname=expression(ln~N[nb]))
seriesplot(q_no,"no",yname=expression(ln~N[o]))
seriesplot(psla$pY,"pY",yname="Pr(Nn>0)",quantiles=F)
seriesplot(psla$pYr,"pYr",yname="Pr(Y>1)",quantiles=F)
seriesplot(psla$pP,"pP",yname="Pr(N>0)",quantiles=F)
seriesplot(q_lYS,"lYS",yname=expression(ln(Y[e]/S[o])))
seriesplot(q_nnK,"nnK",yname=expression(ln(N[n]/K[n])))
seriesplot(q_mnmo,"mnmo",yname=expression(ln(S[n3]/S[o3])))

# Relative change between scenarios ---------------------------------------

scenbase <- "mu1_sd1"
scennew <- c("mu095_sd1", "mu1_sd12", "mu095_sd12")
tpos <- 15 # averaging over range of z due to different iterations
keepsp <- (spvals %in% c("plpa","vuoc"))==F

qal <- list(ns=q_ns,G=q_G,ng=q_ng,Y=q_Y,Ye=q_Ye,nn=q_nn,Sn=q_Sn)

rcl <- lapply(qal,relchange,scenbase=scenbase,scennew=scennew,keepsp=keepsp)
rca <- array(dim=c(dim(rcl[[1]]),length(qal)),
  dimnames=c(dimnames(rcl[[1]]),list(names(qal)))
  )
rca[] <- unlist(rcl)

ccl <- lapply(qal,relchange,scenbase="mu1_sd0",scennew="mu1_sd1",keepsp=keepsp)
cca <- array(dim=c(dim(ccl[[1]]),length(qal)),
  dimnames=c(dimnames(ccl[[1]]),list(names(qal)))
  )
cca[] <- unlist(ccl)

### Quantifying mean-var interactions

mvint <- cca
mvint[] <- NA
dimnames(mvint)[[2]] <- "mvint"
mvint[] <- rca[,"mu095_sd12",] - (rca[,"mu095_sd1",] + rca[,"mu1_sd12",])
  # filling-in only works because mvint has 1 column
rcaa <- abind(rca,cca,mvint,along=2)

### P(nd>1)

pna <- with(psla,abind(pY,pY,along=4))
pna <- aperm(pna,c(4,1,2,3))
  # hack to stack two copies of pYf along 2nd dimension so that relchange works
rpna <- relchange(qlogis(pna),scenbase="mu1_sd0",scennew="mu1_sd1",keepsp=keepsp)

# Vital rate plots for individual runs ------------------------------------

i <- 1
j <- 19
par(mfrow=c(1,1))
plot(density(qlogis(psla$G[i,,j,2])))
lines(density(qlogis(psla$G[i,,j,5])),col="red")

par(mfrow=c(2,1),mar=c(4,4,1,1))
hist(psla$G[i,,j,2],breaks=1000,main="")
hist(psla$G[i,,j,5],breaks=1000,main="")

hist(qlogis(psla$G[i,,16,2]),breaks=1000)
hist(qlogis(psla$G[i,,19,2]),breaks=1000)

# Optimal parameters within species ---------------------------------------

### Germination
parplot(gois$alpha_G,log(psla$ns),expression(alpha[G]),expression(ln(N[s30])),t=30,type="n")
parplot(gois$beta_Gz,log(psla$ns),expression(beta[G]),expression(ln(N[s15])),t=15,type="n")

parplot(plogis(gois$alpha_G),log(psla$ns_p),expression(G),expression(ln(N[s50])),t=50,type="n",ylim=c(6,12))

quantile(gois$tau_mu[,1],probs=c(0.025,0.975))
hist(log(gois$tau_mu[,1]),breaks=1000,xlim=c(-5,5))
quantile(gois$tau_sig[,1],probs=c(0.95))
hist(log(gois$tau_sig[,1]),breaks=1000,xlim=c(0,7))

pb <- gois$beta_Gz[,1]>=0
parplot(gois$tau_mu[pb,],log(psla$ns[pb,,,]),expression(tau[mu]),expression(ln(N[s15])),t=15,type="n",xlim=c(-2,2))
parplot(log(gois$tau_sig[pb,]),log(psla$ns[pb,,,]),expression(ln(tau[sigma])),expression(ln(N[s15])),t=15,type="n")
  # For G sensitivity simulations (too many points to plot)

parplot(gois$alpha_G,log(psla$ns),expression(alpha[G]),expression(ln(N[s15])),t=15,type="n")
parplot(plogis(gois$alpha_G),log(psla$ns),expression(G),expression(ln(N[s15])),t=15,type="n")
parplot(gois$beta_Gz,log(psla$ns),expression(beta[G]),expression(ln(N[s15])),t=15,type="n")
parplot(gois$rho,log(psla$ns),expression(rho[G]),expression(ln(N[s15])),t=15,type="n",xlim=c(-5,5))
parplot(plogis(gois$rho),log(psla$ns),"G",expression(ln(N[s15])),t=15,type="n")

### Extinction
psla$P <- ifelse(psla$ns>0,1,0)

j <- 19
par(mfrow=c(1,1))
plot(gois$tau_mu[pb,j],psla$ns[pb,50,j,5]>0,pch="+",xlim=c(-5,5))
for(i in 1:nclim){
  lines(mysupsmu(gois$tau_mu[pb,j],ifelse(psla$ns[pb,50,j,i]>0,1,0)),col=cols[i])
}
  # high threshold (low G) better 
plot(gois$tau_mu[!pb,j],psla$ns[!pb,50,j,5]>0,pch="+",xlim=c(-5,5))
lines(mysupsmu(gois$tau_mu[!pb,j],ifelse(psla$ns[!pb,50,j,5]>0,1,0)),col="red")

parplot(gois$tau_mu[pb,],psla$P[pb,,,],expression(tau[mu]),expression(P[50]),t=50,type="n",xlim=c(-2,2))
bs <- gois$tau_sig[,1]>=median(gois$tau_sig[,1])

parplot(gois$tau_mu[pb&!bs,],psla$P[pb&!bs,,,],expression(tau[mu]),"P50smallvar",t=50,type="n",xlim=c(-2,2))
parplot(gois$tau_mu[pb&bs,],psla$P[pb&bs,,,],expression(tau[mu]),"P50bigvar",t=50,type="n",xlim=c(-2,2))

parplot(log(gois$tau_sig[pb,]),psla$P[pb,,,],expression(ln(tau[sigma])),expression(P[50]),t=50,type="n")
bm <- gois$tau_mu[,1]>=median(gois$tau_mu[,1])
parplot(log(gois$tau_sig[pb&bm,]),psla$P[pb&bm,,,],expression(ln(tau[sigma])),expression("P50bigmu"),t=50,type="n")
parplot(log(gois$tau_sig[pb&!bm,]),psla$P[pb&!bm,,,],expression(ln(tau[sigma])),expression("P50smallmu"),t=50,type="n")

ba <- gois$alpha_G[,1]>=median(gois$alpha_G[,1])
parplot(gois$beta_Gz[ba,],psla$P[ba,,,],expression(beta[G]),expression(P[50]*alpha[G]>0),t=50,type="n")
parplot(gois$beta_Gz[!ba,],psla$P[!ba,,,],expression(beta[G]),expression(P[50]*alpha[G]<0),t=50,type="n")

bb <- gois$beta_Gz[,1]>=quantile(gois$beta_Gz[,1],probs=0.75)
parplot(gois$alpha_G[bb,],psla$P[bb,,,],expression(alpha[G]),expression(P[50]*beta[G]>5),t=50,type="n")
parplot(gois$alpha_G[!bb,],psla$P[!bb,,,],expression(alpha[G]),expression(P[50]*beta[G]<5),t=50,type="n")

parplot(gois$beta_Gz[ba,],log(psla$ns[ba,,,]),expression(beta[G]),expression(ln(N[s]*alpha[G]>0)),t=15,type="n")
parplot(gois$beta_Gz[!ba,],log(psla$ns[!ba,,,]),expression(beta[G]),expression(ln(N[s]*alpha[G]<0)),t=15,type="n")

parplot(gois$alpha_G[bb,],log(psla$ns[bb,,,]),expression(alpha[G]),expression(ln(N[s]*beta[G]>5)),t=15,type="n")
parplot(gois$alpha_G[!bb,],log(psla$ns[!bb,,,]),expression(alpha[G]),expression(ln(N[s]*beta[G]<5)),t=15,type="n")

parplot(plogis(gois$alpha_G),psla$P,"G",expression(P[50]),t=50,type="n")

### Seed survival
parplot(goi$alpha_m,log(psla$ns),expression(alpha[m]),expression(ln(N[s15])),t=15)
parplot(goi$beta_m,log(psla$ns),expression(beta[m]),expression(ln(N[s15])),t=15)
parplot(log(goi$Kn),log(psla$ns),expression(log(K[N])),expression(ln(N[s15])),t=15)
parplot(log(goi$hn),log(psla$ns),expression(log(h[N])),expression(ln(N[s15])),t=15)
parplot(goi$taun,log(psla$ns),expression(tau[N]),expression(ln(N[s15])),t=15)

### Reproduction
parplot(pri$beta_p[,,1],log(psla$ns),expression(beta[p1]),expression(ln(N[s15])),t=15)
parplot(rsi$beta_r[,,1],log(psla$ns),expression(beta[r1]),expression(ln(N[s15])),t=15)

parplot(pri$beta_p[,,4],log(psla$ns),expression(beta[p4]),expression(ln(N[s15])),t=15)
parplot(rsi$beta_r[,,4],log(psla$ns),expression(beta[r4]),expression(ln(N[s15])),t=15)

parplot(beta_b[,,4],log(psla$ns),expression(beta[b4]),expression(ln(N[s15])),t=15)
parplot(lKy,log(psla$ns),expression(lKy),expression(ln(N[s15])),t=15)
  # Very high K values, not likely to be limiting?

parplot(pri$sig_y_p,log(psla$ns),expression(sigma[py]),expression(ln(N[s15])),t=15)

### Other
parplot(goi$tau_mu,log(psla$ns),expression(tau[mu]),expression(ln(N[s15])),t=15)
parplot(log(goi$tau_sig),log(psla$ns),expression(ln(tau[sigma])),expression(ln(N[s15])),t=15)
parplot(goi$rho,log(psla$ns),expression(rho),expression(ln(n[s15])),t=15)

parplot(goi$alpha_G,log(psla$nsK),expression(alpha[G]),expression(ln(N[sK15])),t=15)
parplot(goi$beta_Gz,log(psla$nsK),expression(beta[G]),expression(ln(N[sK15])),t=15)

# Yearly dynamics ---------------------------------------------------------

cbind(mean=apply(psla$z,3,mean),sd=apply(psla$z,3,sd))

parplot(psla$z,log(psla$r),expression(z),expression(r),t=15,xlim=c(-0.5,1.5))
parplot(psla$z,log(psla$r),expression(z),expression(r),
  type="n",xlim=c(-0.5,1.5),ylim=c(-2.5,0.5))
  # doesn't account for zero-reproduction years
with(psla, parplot(z,lYS,expression(z),expression(ln(Y[e])-ln(S[o])),t=15,xlim=c(-1,1),ylim=c(-5,5)))
with(psla, parplot(z,lYS,expression(z),expression(ln(Y[e])-ln(S[o])),type="n",xlim=c(-1,1),ylim=c(-5,5)))

parplot(log(psla$ns),log(psla$r),expression(ln(N[s])),expression(r),type="n",ylim=c(-2.5,0.5))
parplot(psla$G,log(psla$r),expression(G),expression(r),type="n",ylim=c(-2.5,0.5))
parplot(qlogis(psla$G),log(psla$r),expression(logit(G)),expression(r),t=15,xlim=c(-5,5),ylim=c(-2.5,0.5))

parplot(psla$z,log(psla$Y),expression(z),expression(ln(Y)),t=15,xlim=c(-0.5,1.5))
parplot(psla$z,log(psla$Ye),expression(z),expression(ln(Ye)),t=15,xlim=c(-0.5,1.5))
parplot(psla$z,log(psla$Y),expression(z),expression(ln(Y)),xlim=c(-0.5,1.5),type="n",
  ylim=c(-2.5,5))
parplot(psla$z,log(psla$Ye),expression(z),expression(ln(Ye)),xlim=c(-0.5,1.5),type="n",
  ylim=c(-2.5,5))
parplot(log(psla$ng),log(psla$Ye),expression(ln(N[g])),expression(ln(Ye)),t=15)
parplot(qlogis(psla$G),log(psla$Ye),expression(logit(G)),expression(ln(Ye)),t=15,xlim=c(-5,5))

parplot(log(psla$nn),qlogis(psla$Sn),expression(ln(N[n])),expression(S[n]),type="n")
  # Show how K affects inflection point? Higher K -> CC doesn't cross inflection?
parplot(log(psla$Y),qlogis(psla$Sn),expression(ln(Y)),expression(S[n]),type="n")

# Pr(Y>0) ~ gdens

pardensplot(log(psla$nn),qlogis(psla$Sn),expression(ln(N[n])),expression(S[n]),type="n")
pardensplot(log(psla$Y),qlogis(psla$Sn),expression(ln(Y)),expression(S[n]),type="n")

pardensplot(log(psla$nn),log(psla$nnb),expression(ln(N[n])),expression(N[e]),type="n")
parplot(log(psla$nn),log(psla$nnb),expression(ln(N[n])),expression(N[e]),t=15)
log(apply(goi$taun,2,median)*10*nk)

# Individual runs ---------------------------------------------------------

with(psla,matplot(t(log(ns[1:5,,19,1])),type="l",lty=1))
  # not overcompensating DD because doesn't switch direction every year
  # always some variation because of random year terms in reproduction
with(psla,matplot(t(log(Y[1:5,,19,1])),type="l",lty=1))
abline(h=0,lty=3)
  # basically never falls below replacement rate
with(psla,matplot(t(log(Y[1:5,,19,5])),type="l",lty=1))
abline(h=0,lty=3)
  # often below replacement rate

# Population density distributions at given time --------------------------

densplot(log(psla$ns),"ln(N[s])")
densplot(qlogis(psla$Sn),"logit(S[n])",ylim=c(0,1),vab=apply(goi$So3,2,median))
densplot(psla$mnmo,"log(S[n])-log(S[o])",ylim=c(0,1.5))
densplot(log(psla$nn),"ln(N[n])")
densplot(log(psla$Y),"ln(Y)")

densplot(log(psla$ns_p),"ln(N[s]|N[s]>0)",t=50)
densplot(log(psla$ns_e),"ln(N[s]|N[s]=0)",t=50)

j <- 19
par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(gois$rho[,j],log(psla$ns[,15,j,1]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[,j],log(psla$ns[,15,j,1])),col="red")
plot(gois$rho[,j],log(psla$ns[,15,j,6]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[,j],log(psla$ns[,15,j,6])),col="red")

p1 <- psla$ns[,50,j,1]>0
p2 <- psla$ns[,50,j,6]>0
plot(gois$rho[p1,j],log(psla$ns[p1,15,j,1]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[p1,j],log(psla$ns[p1,15,j,1])),col="red")
plot(gois$rho[p2,j],log(psla$ns[p2,15,j,6]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[p2,j],log(psla$ns[p2,15,j,6])),col="red")

plot(gois$rho[!p1,j],log(psla$ns[!p1,15,j,1]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[!p1,j],log(psla$ns[!p1,15,j,1])),col="red")
plot(gois$rho[!p2,j],log(psla$ns[!p2,15,j,6]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[!p2,j],log(psla$ns[!p2,15,j,6])),col="red")

par(mfrow=c(2,1),mar=c(4,4,1,1))
plot(gois$rho[,j],log(psla$ns[,15,j,2]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[,j],log(psla$ns[,15,j,2])),col="red")
plot(gois$rho[,j],log(psla$ns[,15,j,5]),pch="+",ylim=c(0,15))
lines(mysupsmu(gois$rho[,j],log(psla$ns[,15,j,5])),col="red")

# Parameter correlations within a given species ---------------------------

j <- 17
with(goi,plot(log(Kn[,j])~log(hn[,j])))
with(goi,plot(log(alpha_m[,j])~log(beta_m[,j])))

parplot(log(psla$nn),qlogis(psla$Sn),expression(ln(N[n])),expression(S[n]),type="n")

# Changes between climate scenarios ---------------------------------------

### Pairs correlations

pairplot("popchange_pairs_vitalratechange_newenv",rca,2)
pairplot("popchange_pairs_vitalratechange_consenv",cca,2)
pairplot("popchange_pairs_allclimscenarios_",rcaa,3)
pairplot("popchange_pairs_allclimscenarios_alltraits",rcata,3,w=21,h=21)

### Within species

Kn <- goi$Kn*nk/10*tau_s
diffplot(log(Kn),log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(K[n]),yname=expression(N[s]))
diffplot(log(Kn),log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(K[n]),yname=expression(N[s]))

diffplot(goi$alpha_G,log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(alpha[G]),yname=expression(N[s]))
diffplot(goi$alpha_G,log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(alpha[G]),yname=expression(N[s]))

diffplot(goi$alpha_m,log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(alpha[m]),yname=expression(N[s]))
diffplot(goi$alpha_m,log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(alpha[m]),yname=expression(N[s]))
diffplot(qlogis(exp(-exp(goi$alpha_m))),log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(S[o]),yname=expression(N[s]))
diffplot(qlogis(exp(-exp(goi$alpha_m))),log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(S[o]),yname=expression(N[s]))
  # alpha_m determines K, but not particularly strong predictor for dN?

diffplot(qlogis(psla$G),log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(dG),yname=expression(dN[s]),xdiff=T)
diffplot(qlogis(psla$G),log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(dG),yname=expression(dN[s]),xdiff=T)
  # in general, greater reductions in G don't -> greater reductions in N

diffplot(log(psla$Y),log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(dY),yname=expression(dN[s]),xdiff=T)
diffplot(log(psla$Y),log(psla$ns),refcl="mu1_sd1",t=15,
  xname=expression(dY),yname=expression(dN[s]),xdiff=T)
  # y reductions explain population decline due to mean environmental change
  # not an artefact, because y is measured after Ns (?)

diffplot(qlogis(psla$Sn),log(psla$ns),refcl="mu1_sd0",t=15,
  xname=expression(dS[n]),yname=expression(dN[s]),xdiff=T)

### Median population sizes

diffplot(log(psla$Y),log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(dY),yname=expression(dN[s]),xdiff=T)
diffplot(qlogis(psla$G),log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(dG),yname=expression(dN[s]),xdiff=T)
diffplot(log(psla$Y),log(psla$ns),refcl="mu1_sd1",t=NULL,
  xname=expression(dY),yname=expression(dN[s]),xdiff=T)
diffplot(qlogis(psla$G),log(psla$ns),refcl="mu1_sd1",t=NULL,
  xname=expression(dG),yname=expression(dN[s]),xdiff=T)
  # but only 50 timesteps, so could just reflect which sims happened to 
  # get lots of rainfall?

diffplot(pri$beta_p[,,1],log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(beta[p1]),yname=expression(dN[s]))
diffplot(rsi$beta_r[,,1],log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(beta[r1]),yname=expression(dN[s]))
diffplot(goi$alpha_G,log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(alpha[G]),yname=expression(dN[s]))
  # parameters themselves don't seem to have big influence on outcome (!)
diffplot(goi$tau_mu,log(psla$ns),refcl="mu1_sd0",t=NULL,
  xname=expression(tau[mu]),yname=expression(dN[s]),
  xlim=quantile(psla$z,probs=c(0.025,0.975))
  )
  # no obvious benefit of low germination threshold
diffplot(qlogis(exp(-exp(goi$alpha_m))),log(psla$ns),refcl="mu1_sd1",t=NULL,
  xname=expression(S[o]),yname=expression(dN[s]))
  # Seed survival seems to be more important for mean than for variance adaptation

Kn <- goi$Kn*nk/10*tau_s
dns <- with(psla,log(ns[,15,,"mu095_sd1"])-log(ns[,15,,"mu1_sd1"]))
dns[is.nan(dns)] <- NA
mdns <- apply(dns,2,median,na.rm=T)
mKn <- apply(Kn,2,median)
So <- apply(exp(-exp(goi$alpha_m)),2,median)
plot(mdns~log(mKn))
plot(mdns~qlogis(So))
exclsp <- c(18,22)
plot(mdns[-exclsp]~log(mKn[-exclsp]))
  # a few species with low K seem to decline less strongly
plot(mdns[-exclsp]~qlogis(So[-exclsp]))
  # no obvious pattern

Kn <- goi$Kn*nk/10*tau_s
dns <- with(psla,
  log(apply(ns[,,,"mu095_sd1"],c(1,3),median)) 
  - log(apply(ns[,,,"mu1_sd1"],c(1,3),median))
)
mdns <- apply(dns,2,median)
mKn <- apply(Kn,2,median)
plot(mdns~log(mKn))

j <- 19
plot(log(Knarr[,1,j]),dns[,j,6])
lines(mysupsmu(log(Knarr[,1,j]),dns[,j,6]),col="red")

plot(qlogis(psla$So[,1,j,6]),dns[,j,6])
lines(mysupsmu(qlogis(psla$So[,1,j,6]),dns[,j,6]),col="red")

dc_Sn <- qlogis(psla$Sn) - rep(qlogis(psla$Sn[,,,1]),nclim)
dc_ns <- log(psla$ns) - rep(log(psla$ns[,,,1]),nclim)
# dc = difference current; df = difference future

parplot(dc_Sn,dc_ns,expression(dS[n]),expression(dN[s]),t=15)

with(goi,plot(alpha_G[,j]~alpha_m[,j]))
with(goi,plot(apply(alpha_G,2,median)~apply(alpha_m,2,median)))
with(goi,cor.test(apply(alpha_G,2,median),apply(alpha_m,2,median)))

# Climate distributions ---------------------------------------------------

pdf(paste0("Plots/zdists",format(Sys.Date(),"%d%b%Y"),".pdf"),width=4.5,height=4.5)

plot(density(psla$z[,,2]),xlim=c(-2,2),main="",col=cols[2])
for(i in 3:(nclim-1)){
  lines(density(psla$z[,,i]),col=cols[i])
}
abline(v=psla$z[1,1,1],col=cols[1],lty=3)
abline(v=psla$z[1,1,6],col=cols[6],lty=3)

plot(density(tau_p*exp(psla$z[,,2])),xlim=c(0,250),ylim=c(0,0.011),main="",col=cols[2])
for(i in 3:(nclim-1)){
  lines(density(tau_p*exp(psla$z[,,i])),col=cols[i])
}
abline(v=tau_p*exp(psla$z[1,1,1]),col=cols[1],lty=3)
abline(v=tau_p*exp(psla$z[1,1,6]),col=cols[6],lty=3)

dev.off()

