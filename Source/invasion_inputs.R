### Assemble inputs for invasion analyses ###

# Input parameters --------------------------------------------------------

nj <- 22
nit <- 100
nr <- 100 # number of repeated invasions
nt <- 100
nb <- 25 # number of "burn-in" timesteps to stabilise resident dynamics

uncertainty <- T
plasticity <- T
nstart <- 1

cm <- c(0,0,0)
cs <- c(-1,0,1) # powers
nc <- length(cm)

# Basic parameters --------------------------------------------------------

rho <- 0.82

# Climate -----------------------------------------------------------------

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

# Model parameters --------------------------------------------------------

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
  # already permuted

### SELECT MEDIAN PARAMETERS

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

posd <- expand.grid(ipos=1:nit,jpos=1:nj,mpos=1:nc)

itot <- 10^4
ie <- posd$ipos
ce <- with(posd, jpos * itot + ipos)

pard <- with(pl, data.frame(
  beta_p1=pr$beta_p[,,1][ce],
  beta_p2=pr$beta_p[,,2][ce],
  beta_p3=pr$beta_p[,,3][ce],
  beta_p4=pr$beta_p[,,4][ce],
  beta_r1=rs$beta_r[,,1][ce],
  beta_r2=rs$beta_r[,,2][ce],
  beta_r3=rs$beta_r[,,3][ce],
  beta_r4=rs$beta_r[,,4][ce],
  sig_y_p=pr$sig_y_p[ce],
  sig_y_r=rs$sig_y_r[ce],
  sig_s_g=gs$sig_s_g[ce],
  sig_s_p=pr$sig_s_p[ie],
  sig_s_r=rs$sig_s_r[ie],
  sig_o_p=pr$sig_o_p[ie],
  phi_r=rs$phi_r[ie],
  theta_g=gs$theta_g[ce],
  m0=exp(go$alpha_m[ce]),
  m1=exp(go$beta_m[ce])
  ))

pd <- cbind(posd,pard)

### PARAMS FOR SENSITIVITY ANALYSES
# Creates grid of starting alpha and beta values

pd$am0 <- runif(nit,-5,5)[pd$ipos]
if(plasticity==T) pd$bm0 <- runif(nit,-5,5)[pd$ipos]
if(plasticity==F) pd$bm0 <- 0

# Sims --------------------------------------------------------------------

# Each core focuses on one climate
# Iterations for one climate can be split up over multiple cores
# (controlled by cpc)

# maml <- as.list(c(1,1,mpam,1,mpam,mpam))
# msdl <- as.list(c(0,1,1,mpsd,mpsd,0))

# scaling mean log rainfall (zamo) only works because sign stays the same

nclim <- length(maml)
cpc <- 10 # CORES per CLIMATE (assumed equal for resident and invader)
ncores <- nclim*cpc
mpos <- rep(1:nclim,each=cpc)

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


# Simulate normalised climate series --------------------------------------

library(MASS)
set.seed(1)
ntt <- nt * nit
zw <- mvrnorm(n=ntt, mu=c(0,0), Sigma=matrix(c(1,rep(rho,2),1),nr=2,nc=2))
eps_y_p <- rnorm(ntt,0,sig_y_p)
eps_y_r <- rnorm(ntt,0,sig_y_r)

zen <- array(dim=c(nt,nit,4),
             dimnames=list(NULL,NULL,c("z","w","p","r"))
             )
zen[,,1] <- zw[,1]
zen[,,2] <- zw[,2]
zen[,,3] <- eps_y_p
zen[,,4] <- eps_y_r

# Write simulation input RDS files ----------------------------------------

apptext <- paste("j",nj,
                 "i",ni,
                 "r",nr,
                 "t",nt,
                 "b",nb
                 "n",nstart,
                 "cm"
                 
                 
i",nj,
                 sep="_")
                  nit <- 100
                  nr <- 100 # number of repeated invasions
                  nt <- 100
                  nb <- 25 # number of "burn-in" timesteps to stabilise resident dynamics
                  nc <- 10
                  
                  uncertainty <- T
                  plasticity <- T
                  nstart <- 1
                  
                  cm <- c(0,0,0)
                  cs <- c(-1,0,1) # powers
                  nm <- length(maml)
_nk",nk,"_")

saveRDS(list(zen=zen,ipd),paste0("Sims/",zen,"_",cur_date,".rds"))

# PARAMETERS
# CLIMATE VALUES

# ES G simulations --------------------------------------------------------

^msd 

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

# Other -------------------------------------------------------------------

source("Source/invasion_functions.R")
source("Source/figure_functions.R")
source("Source/prediction_functions.R")
source("Source/trait_functions.R")

library(plyr)
library(reshape2)
library(parallel)
library(abind)
library(RColorBrewer)