### Iterated simulations of resdient population dynamics and ###
### invasions by new strategies                              ###

library(plyr)
library(reshape2)
library(parallel)
library(abind)
library(RColorBrewer)

### LOAD DATA 

source("Source/invasion_functions.R")
source("Source/figure_functions.R")
source("Source/prediction_functions.R")
source("Source/trait_functions.R")

msy <- read.csv("Output/msy_15Jan2016.csv",header=T)
Tvalues <- read.csv("Output/Tvalues_31Jul2015.csv",header=T)

# Define params -----------------------------------------------------------

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])
# T in years

tau_s <- 100  # adujstment for seeds
tau_d <- 100	# adjustment for density
tau_p <- 100	# adjustment for rainfall
  # NB: these T and tau values don't affect sims, which use default arguments

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
nyear <- length(unique(msy$year))

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
  # go = readRDS("Models/go_pars_tdistpois_naspecies_noerr_noGDD_loglik_RICKER_15Oct2017.rds"),
	gs = readRDS("Models/gnzhh_onhh_pars_medians_26Oct2015.rds"),
		# gs = g site level
		# source script: venable_Stan_GO_descriptive_gnzhh_onhh_26Oct2015
		# uses tau_s = 100
		# but tau actually irrelevant because all multiplicative?
	pr = readRDS("Models/pr_pars_yearhet_squared_pc_02Mar2016.rds"),
	rs = readRDS("Models/rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")
	)

### MEDIAN PARAMETERS

uncertainty <- T
plasticity <- T

pls <- with(pl, list(
  beta_p=pr$beta_p,
  beta_r=rs$beta_r,
  sig_y_p=pr$sig_y_p,
  sig_y_r=rs$sig_y_r,
  sig_s_g=gs$sig_s_g,
  sig_s_p=pr$sig_s_p,
  sig_s_r=rs$sig_s_r,
  sig_o_p=pr$sig_o_p,
  phi_r=rs$phi_r,
  theta_g=gs$theta_g,
  m0=exp(go$alpha_m),
  m1=exp(go$beta_m)
))

if(uncertainty==F){
  for(i in 1:length(pls)){
    ndim <- length(dim(pls[[i]])) 
    if(ndim==1){
      pls[[i]] <- rep(median(pls[[i]]),length(pls[[i]]))
    }
    if(ndim>1){
      keepdim <- 2:ndim # remove first dim, which is always iter
      eachrep <- dim(pls[[i]])[1] # select iterdim
      pls[[i]][] <- rep(apply(pls[[i]],keepdim,median),each=eachrep)
    }
    # if is.null(ndim), do nothing
  }
}

### PARAMS FOR SENSITIVITY ANALYSES
# Creates grid of starting alpha and beta values

nit <- 100 # 100

if(plasticity==F){
  Gsens <- data.frame(
    # alpha_G=qlogis(seq(0.001,0.999,length.out=nit)),
    alpha_G=seq(-5,5,length.out=nit),
    beta_Gz=rep(0,nit)
  )
}

if(plasticity==T){
  Gsens <- expand.grid(
    alpha_G = seq(-5,5,length.out=sqrt(nit)),
    beta_Gz = seq(-5,5,length.out=sqrt(nit))
  )
}

# if(plasticity==T){
#   # assumes positive thresholds
#   tmrange <- c(-1.5,0.5)
#   tsrange <- c(-3,1)
#   # transformed climate range approx. -1 -> 1
#   neach <- 10
#   tau_mu <- seq(tmrange[1],tmrange[2],length.out=neach)
#   tau_sd <- 2^seq(tsrange[1],tsrange[2],length.out=neach)
#   Gsens <- expand.grid(tau_mu=tau_mu,tau_sd=tau_sd)
#   np <- nrow(Gsens) # = neach^2
#   Gsens$alpha_G <- with(Gsens,godalpha_f(tau_mu,tau_sd))
#   Gsens$beta_Gz <- with(Gsens,godbeta_f(tau_sd))
# }

set.seed(1)
Gsens <- Gsens[sample(1:nit,nit,replace=FALSE),]
  # random starting points 
  # -> does't matter how distributed over different climates

pls$am0 <- Gsens$alpha_G
pls$bm0 <- Gsens$beta_Gz

# Sims --------------------------------------------------------------------

# Each core focuses on one climate
# Iterations for one climate can be split up over multiple cores
# (controlled by cpc)

# maml <- as.list(c(1,1,mpam,1,mpam,mpam))
# msdl <- as.list(c(0,1,1,mpsd,mpsd,0))
maml <- as.list(c(1,1,1)) 
msdl <- as.list(c(1/mpsd,1,mpsd)) 
  # scaling mean log rainfall (zamo) only works because sign stays the same

nclim <- length(maml)
cpc <- 25 # CORES per CLIMATE (assumed equal for resident and invader)
ncores <- nclim*cpc
mpos <- rep(1:nclim,each=cpc)

# nstart <- 1
nr <- 100 # 100 # number of repeated invasions
nt <- 1025 # 125 # 10050 
nb <- 25  # number of "burn-in" timesteps to stabilise resident dynamics
nj <- 22
  # min invader iterations per core = nr * nit
nk <- 0  

iseq <- 1:nit

cpos <- rep(1:cpc,times=nclim)
cipos <- rep(1:cpc,each=nit/cpc)

itersetl <- split(iseq,cipos)
  # requires that ni < maxiter
ni <- length(itersetl[[1]]) # iterations PER CORE for RESIDENT simulations

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

# ES G simulations --------------------------------------------------------

system.time({
CL = makeCluster(ncores)
clusterExport(cl=CL, c(
  "BHS","RICKERS",
  "logitnorm","logitmean","logitnormint",
  "nbtmean","nbtnorm","nbtlnmean","fnn","pradj","sprinkle",
  "fixG","ressim","invade_infinite","invade_finite",
  "evolve","multievolve",
  "pls","plasticity",
  "ni","nj","nr","nt","nb","nk",
  "zamo","zsdo","wamo","wsdo",
  "mpos","maml","msdl","cpos",
  "Tvalues","cnames_bycore",
  "itersetl"
  )) 
parLapply(CL, 1:ncores, function(n){
	mam <- maml[[mpos[n]]]
	msd <- msdl[[mpos[n]]]
	iset <- itersetl[[cpos[n]]]
	with(pls, multievolve(
	  ni=ni,nj=nj,nr=nr,nt=nt,nb=nb,nk=nk,
	  zam=zamo+mam*zsdo,zsd=zsdo*msd,
	  wam=wamo+mam*wsdo,wsd=wsdo*msd,
	  beta_p=beta_p[iset,,],beta_r=beta_r[iset,,],
	  sig_y_p=sig_y_p[iset,],sig_y_r=sig_y_r[iset,],
	  sig_s_g=sig_s_g[iset,],sig_s_p=sig_s_p[iset],sig_s_r=sig_s_r[iset],
	  sig_o_p=sig_o_p[iset],phi_r=phi_r[iset],theta_g=theta_g[iset,],
	  m0=m0[iset,],m1=m1[iset,],
	  am0=am0[iset],bm0=bm0[iset],
	  DDFUN=BHS,
	  Sg=1,
	  smut_a=5,smut_b=ifelse(plasticity==TRUE,5,0),
		savefile=paste0("ESS_p",
		  plasticity,
		  "_nk",nk,"_",
		  cnames_bycore[n]
		  ) 
		))
	})
stopCluster(CL)
})
  # OLD METHOD
  # 54 hours (25 cores, long time series [nt=1025])
  # 23 hours (25 cores, short time series with 100 iterations [ni])
    # 46 hours with low-mid-high variability - but most finished in 24 hours
    # 3.3 hours with nt=75, ni=50, nr=50
    # but non-spatial = 2.25 mins!
    # 1.3 hours: ni=100, nr=100, nt=1025, nk=0, cpc=25
    #   (timing very consistent: all finished within 10 mins of each other)
  # 48 hours (25 cores, same but 200 instead of 100 invasions [nr])
  # 100 hours (25 cores, nk=10000, nit=50, nt=125, nr=100; 
  #   2/25 cores didn't finish)
  # 33 hours (10 cores, nk=Inf, nit=100, nr=100, nt=125)
  #
  # FINITE - NEW METHOD
  # 2 hours: 20 cores, nk=100, nit=100, nr=100, nt=125
  # stopped at 100 hours: 10 cores, nk=10000, nit=100, nr=100, nt=125
  #   only 8/30 finished, 7 of which were high-variance (am=1,sd=12)
  #
  # FINITE - NEW METHOD 2 (all  nit=100, nr=100, nt=125)
  # 39 mins: 20 cores, nk=100
  # 50 hours: 25 cores, nk=10000 (but 10 missing, 1 unfinished)
  # A few hours: 1 core, nk=10
  # 18 hours: 2 cores, nk=1000
  # 
  # SHORT SIMS
  # ni=50, nr=100, nt=75, cpc=25
  # nk=1: 53 secs
  # nk=10: 122 secs
  # nk=100: 1010 secs
  # nk=1000, ni=50, nr=50, cpc=2: 3.23 hours
  #   nk=Inf, cpc=2: 1 hour

# write verbose form that allows convergence to be checked?

# Redo missed sims --------------------------------------------------------

# reclim <- c(1,1,1,1,1,2,2,2,2,3,3,3)
# reset <- c(4,6,9,14,25,6,9,14,24,6,9,25)
# # rerun <-  c("mu1_sd081_s4","mu1_sd081_s6","mu1_sd081_s9",
# #   "mu1_sd081_s14","mu1_sd081_s25","mu1_sd1_s6",  
# #   "mu1_sd1_s9","mu1_sd1_s14","mu1_sd1_s24",  
# #   "mu1_sd12_s6","mu1_sd12_s9","mu1_sd12_s25"
# # )
# reit <- unlist(lapply(itersetl[reset],function(x){
#   split(x,1:4)
# }), recursive = TRUE)
# renames <- paste(
#   rep(cnames_unique,table(reclim)),
#   rep(reset,each=4),
#   rep(1:4,times=length(reset)
#   ),sep="_")
# 
# ncoresre <- length(reclim)*4
# 
# system.time({
#   CL = makeCluster(ncoresre)
#   clusterExport(cl=CL, c(
#     "reclim","reset","reit","renames","ncoresre",
#     "BHS","RICKERS",
#     "logitnorm","logitmean","logitnormint",
#     "nbtmean","nbtnorm","nbtlnmean","fnn","pradj","sprinkle",
#     "fixG","ressim","invade_infinite","invade_finite",
#     "evolve","multievolve",
#     "pls", 
#     "ni","nj","nr","nt","nb","nk",
#     "zamo","zsdo","wamo","wsdo",
#     "mpos","maml","msdl","cpos",
#     "Tvalues","cnames_bycore",
#     "itersetl"
#   )) 
#   parLapply(CL, 1:ncoresre, function(n){
#     mam <- maml[[rep(reclim,each=4)[n]]]
#     msd <- msdl[[rep(reclim,each=4)[n]]]
#     iset <- reit[n]
#     with(pls, multievolve(
#       ni=1,nj=nj,nr=nr,nt=nt,nb=nb,nk=nk,
#       zam=zamo+mam*zsdo,zsd=zsdo*msd,
#       wam=wamo+mam*wsdo,wsd=wsdo*msd,
#       beta_p=beta_p[iset,,],beta_r=beta_r[iset,,],
#       sig_y_p=sig_y_p[iset,],sig_y_r=sig_y_r[iset,],
#       sig_s_g=sig_s_g[iset,],sig_s_p=sig_s_p[iset],sig_s_r=sig_s_r[iset],
#       sig_o_p=sig_o_p[iset],phi_r=phi_r[iset],theta_g=theta_g[iset,],
#       m0=m0[iset,],m1=m1[iset,],
#       am0=am0[iset],bm0=bm0[iset],
#       DDFUN=BHS,
#       Sg=1,
#       savefile=paste0("ESS_finite_nk",nk,"_",renames[n]) 
#     ))
#   })
#   stopCluster(CL)
# })

# Read resident simulations back in ---------------------------------------

### Small RAM read-in

# cnames_bycore_small <- cnames_bycore[grep("mu1_sd1",cnames_bycore)]
# ncores_small <- length(cnames_bycore_small)
# psl <- as.list(rep(NA,ncores_small))
# for(n in 1:ncores_small){
#   psl[[n]] <- readRDS(paste0(cnames_bycore_small[n],"_28Apr2016.rds"))
#   }
# names(psl) <- cnames_bycore_small

apptext <- paste0("_p",plasticity,"_nk",nk,"_")
psl <- as.list(rep(NA,ncores))
dir <- paste0(getwd(),"/Sims/")
files <- paste0(dir,list.files(dir))

for(n in 1:ncores){
  curname <- paste0("Sims/ESS",apptext,cnames_bycore[n],"_10Dec2017.rds")
  # curname <- paste0("Sims/ESS_finite_nk",nk,"_",cnames_bycore[n],"_08Dec2017.rds")
  # curname <- paste0("Sims/ESS_infinite_spatial_",cnames_bycore[n],"_28Nov2017.rds") 
    # 06Dec 08Dec
  finished <- grep(curname,files)
  if(length(finished)!=0){
    psl[[n]] <- readRDS(curname)
    }
}
names(psl) <- cnames_bycore

preplace <- function(x){
  dnf <- is.na(x)
  x[dnf] <- x[!dnf][1:sum(dnf)]
  return(x)
}
psls <- unlist(lapply(split(psl,mpos),preplace),recursive=FALSE)
  # sequentially replaces sims that didn't finish
  # requires that for each climate, there are more finished than missing cores

psla <- simcombine(psl,nclim=nclim,cpc=cpc) # psls

psla_out <- list(
  maml=maml,msdl=msdl,
  nk=nk,nt=nt,nb=nb,nr=nr,ni=nit,nj=nj,
  pls=pls,psla=psla
  )
saveRDS(psla_out,paste0("Sims/ESSall",apptext,format(Sys.Date(),"%d%b%Y"),".rds"))

# ES G plots --------------------------------------------------------------

purples <- brewer.pal(9,"Purples")[5] 
blues <- brewer.pal(9,"Blues")[5] 
greens <- brewer.pal(9,"Greens")[5] 
oranges <- brewer.pal(9,"Oranges")[5]
reds <- brewer.pal(9,"Reds")[5] 
greys <- brewer.pal(9,"Greys")[5] 

cols <- c(purples,blues,greens,oranges,reds,greys)
ltys <- c(3,1,3)

colledgetext <- cnames_unique
detledgetext <- c(
  paste0("nstart=",1),
  paste0("ni=",nit),
  paste0("nr=",nr),
  paste0("nb=",nb),
  paste0("nt=",nt),
  paste0("nk=",nk),
  paste0("nj=",nj)
)

colpos <- 1:nclim 

tran <- 25
cols_rgb <- col2rgb(cols)
trancols <- rgb(
  red=cols_rgb[1,],
  green=cols_rgb[2,],
  blue=cols_rgb[3,],
  alpha=tran,
  maxColorValue = 255
)

alpha_G <- as.vector(psla$am)
beta_Gz <- as.vector(psla$bm)
# alpha_G <- with(psl,c(mu1_sd081_s1$am,mu1_sd1_s1$am,mu1_sd12_s1$am))
# beta_Gz <- with(psl,c(mu1_sd081_s1$bm,mu1_sd1_s1$bm,mu1_sd12_s1$bm))

nw <- 100
wseq <- seq(-2,2,length.out=nw)
Gw <- array(dim=c(nw,ni*cpc,nj,nclim))
Gw[] <- plogis(
  matrix(rep(alpha_G,each=nw),nr=nw,nc=ni*cpc*nj*nclim)
  + outer(wseq,beta_Gz,"*")
)
qGw <- aperm(
  apply(Gw,c(1,3,4),quantile,prob=c(0.25,0.50,0.75),na.rm=T), 
  c(2,3,1,4)
)

nob <- nrow(pl$go$alpha_G)
Gwob <- array(dim=c(nw,nob,nj))
for(j in 1:nj){
  Gwob[,,j] <- plogis(
    matrix(rep(pl$go$alpha_G[,j],each=nw),nr=nw,nc=nob)
    + outer(wseq,pl$go$beta_Gz[,j],"*")
  )
}
qGob <- aperm(
  apply(Gwob,c(1,3),quantile,prob=c(0.25,0.50,0.75),na.rm=T), 
  c(2,3,1)
)

wam <- sapply(maml,function(x) zam=wamo+x*wsdo - log(tau_p))
wsd <- sapply(msdl,function(x) zam=wsdo*x)
qw <- rbind(wam-1.96*wsd,wam+1.96*wsd)

trangrey <- rgb(red=190,green=190,blue=190,alpha=0.25,maxColorValue = 255)

# Quantile ES G plots -----------------------------------------------------

pdf(paste0("Plots/ESS",apptext,format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)

plotsetup()

for(j in 1:nspecies){
  matplot(wseq,qGw[,j,,1],type="l",lty=ltys,ylim=c(0,1),col=cols[colpos][1])
  xx <- rep(qw[,colpos[1]],each=2)
  yy <- c(0,1,1,0)
  polygon(xx, yy, col=trancols[colpos][1],border=NA)
  
  for(m in 1:nclim){
    if(m!=1){
      matplot(wseq,qGw[,j,,m],type="l",lty=ltys,add=T,col=cols[colpos][m])
      #[c(1,2,4)][m]
    }
    xx <- rep(qw[,colpos[m]],each=2)
    yy <- c(0,1,1,0)
    polygon(xx, yy, col=trancols[colpos][m],border=NA)
  }
  
  matplot(wseq,qGob[,j,],type="l",lty=ltys,col="black",add=T)
  
  lettlab(j)
  
  if(j %in% 19:23) addxlab("w") 
  if(j %in% seq(1,23,4)) addylab("G") 
}

addledge(ltext=colledgetext[colpos],col=cols[colpos],lty=1) #[c(1,2,4)]
addledge(ltext=detledgetext)

dev.off()

# Example ES G plots ------------------------------------------------------

nd <- 25 # number of draws to plot

pdf(paste0("Plots/Gdraws",apptext,format(Sys.Date(),"%d%b%Y"),".pdf"),
    width=plotwidth,height=plotheight)

for(m in 1:nclim){

  plotsetup()

  for(j in 1:nspecies){
    matplot(wseq,Gw[,1:nd,j,m],type="l",lty=1,ylim=c(0,1),col=cols[colpos][m])
    xx <- rep(qw[,colpos[m]],each=2)
    yy <- c(0,1,1,0)
    polygon(xx, yy, col=trancols[colpos][m],border=NA)
    
    # for(m in 1:nclim){
    #   if(m!=1){
    #     matplot(wseq,Gw[,1:nd,j,m],type="l",lty=1,add=T,col=cols[colpos][m])
    #     #[c(1,2,4)][m]
    #   }
    #   xx <- rep(qw[,colpos[m]],each=2)
    #   yy <- c(0,1,1,0)
    #   polygon(xx, yy, col=trancols[colpos][m],border=NA)
    # }
    
    matplot(wseq,qGob[,j,1],type="l",lty=1,col="black",add=T)
    
    lettlab(j)
    
    if(j %in% 19:23) addxlab("w") 
    if(j %in% seq(1,23,4)) addylab("G") 
  }
  
addledge(ltext=colledgetext[colpos],col=cols[colpos],lty=1) #[c(1,2,4)]
addledge(ltext=detledgetext)

}

dev.off()
