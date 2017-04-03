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

# Derived parameters ------------------------------------------------------

### From input parameters

goi$tau_mu <- with(goi, godmean_f(alpha_G,beta_Gz) )
goi$tau_sig <- with(goi, godvar_f(beta_Gz) )
goi$rho <- with(goi, alpha_G + beta_Gz*log(zam/tau_p))

goi$m0 <- exp(goi$alpha_m)
goi$m1 <- exp(goi$beta_m)
goi$Kn <- with(goi, Kncalc(m0,m1,T3))
goi$hn <- with(goi, hncalc(m0,m1,T3))

### From simulations

psla$Y <- with(psla,nn/ng)
psla$Ye <- with(psla,nnb/ng)
aNA <- array(dim=c(nit,1,nj,nclim))
psla$r <- log(abind(psla$ns[,-1,,],aNA,along=2)) - log(psla$ns)
  # popualtion growth rates (year t = growth to next year
psla$pY <- apply(psla$nn,2:4,function(x) sum(x>0)/length(x))
  # probability of at least one new seed
  # (could also use to calculate extinction risk)
psla$lYvS <- with(psla,log(Ye)-log(So))

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
q_lYvS <- seriesquant(psla$lYvS)

# Plot results ------------------------------------------------------------

seriesplot(q_ns,"ns",yname=expression(ln(N[s])))
seriesplot(q_nn,"nn",yname=expression(ln(N[n])))
seriesplot(q_G,"G",yname="logit(G)")
seriesplot(q_Sn,"Sn",yname="logit(Sn)")
seriesplot(q_Y,"Y",yname="ln Y")
seriesplot(q_Ye,"Ye",yname="ln Yeff")
seriesplot(q_nnb,"nnb",yname=expression(ln~N[nb]))
seriesplot(q_no,"no",yname=expression(ln~N[o]))
seriesplot(psla$pY,"pY",yname="Pr(Y>0)",quantiles=F)
seriesplot(q_lYvS,"lYvS",yname=expression(ln(Y[e]/S[o])))

# Relative change between scenarios ---------------------------------------

scenbase <- "mu1_cv1"
scennew <- c("mu081_cv1", "mu1_cv12", "mu081_cv12")
tpos <- 15 # averaging over range of z due to different iterations
keepsp <- (spvals %in% c("plpa","vuoc"))==F

qal <- list(ns=q_ns,G=q_G,ng=q_ng,Y=q_Y,Ye=q_Ye,nn=q_nn,Sn=q_Sn)

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

pna <- with(psla,abind(pY,pY,along=4))
pna <- aperm(pna,c(4,1,2,3))
  # hack to stack two copies of pYf along 2nd dimension so that relchange works
rpna <- relchange(qlogis(pna),scenbase="mu1_cv0",scennew="mu1_cv1",keepsp=keepsp)

# Pairs correlations ------------------------------------------------------

pairplot("popchange_pairs_vitalratechange_newenv",rca,2)
pairplot("popchange_pairs_vitalratechange_consenv",cca,2)
pairplot("popchange_pairs_allclimscenarios_",rcaa,3)
pairplot("popchange_pairs_allclimscenarios_alltraits",rcata,3,w=21,h=21)

# Optimal parameters within species ---------------------------------------

parplot(goi$alpha_G,log(psla$ns),expression(alpha[G]),expression(ln(N[s15])),t=15)
parplot(goi$beta_Gz,log(psla$ns),expression(beta[G]),expression(ln(N[s15])),t=15)
parplot(goi$alpha_m,log(psla$ns),expression(alpha[m]),expression(ln(N[s15])),t=15)
parplot(goi$beta_m,log(psla$ns),expression(beta[m]),expression(ln(N[s15])),t=15)
parplot(log(goi$Kn),log(psla$ns),expression(log(K[N])),expression(ln(N[s15])),t=15)
parplot(log(goi$hn),log(psla$ns),expression(log(h[N])),expression(ln(N[s15])),t=15)

parplot(pri$beta_p[,,1],log(psla$ns),expression(beta[p1]),expression(ln(N[s15])),t=15)
parplot(pri$beta_p[,,4],log(psla$ns),expression(beta[p4]),expression(ln(N[s15])),t=15)
parplot(pri$sig_y_p,log(psla$ns),expression(sigma[py]),expression(ln(N[s15])),t=15)
parplot(rsi$beta_r[,,1],log(psla$ns),expression(beta[r1]),expression(ln(N[s15])),t=15)

parplot(goi$tau_mu,log(psla$ns),expression(tau[mu]),expression(ln(N[s15])),t=15)
parplot(log(goi$tau_sig),log(psla$ns),expression(ln(tau[sigma])),expression(ln(N[s15])),t=15)
parplot(goi$rho,log(psla$ns),expression(rho),expression(ln(n[s15])),t=15)

parplot(goi$alpha_G,log(psla$nsK),expression(alpha[G]),expression(ln(N[sK15])),t=15)
parplot(goi$beta_Gz,log(psla$nsK),expression(beta[G]),expression(ln(N[sK15])),t=15)

# Yearly dynamics ---------------------------------------------------------

cbind(apply(psla$z,3,mean),apply(psla$z,3,sd))

parplot(psla$z,log(psla$r),expression(z),expression(r),t=15,xlim=c(-0.5,1.5))
parplot(psla$z,log(psla$r),expression(z),expression(r),
  type="n",xlim=c(-0.5,1.5),ylim=c(-2.5,0.5))
  # doesn't account for zero-reproduction years
with(psla, parplot(z,lYvS,expression(z),expression(ln(Y[e])-ln(S[o])),t=15,xlim=c(-1,1),ylim=c(-5,5)))
with(psla, parplot(z,lYvS,expression(z),expression(ln(Y[e])-ln(S[o])),type="n",xlim=c(-1,1),ylim=c(-5,5)))

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

parplot(log(psla$nn),log(psla$Sn),expression(ln(N[n])),expression(S[n]),type="n")

# Pr(Y>0) ~ gdens


# Parameters correlations within species ----------------------------------

j <- 17
with(goi,plot(log(Kn[,j])~log(hn[,j])))
with(goi,plot(log(alpha_m[,j])~log(beta_m[,j])))

# Climate distributions ---------------------------------------------------

pdf(paste0("Plots/zdists",format(Sys.Date(),"%d%b%Y"),".pdf"),width=4.5,height=4.5)

plot(density(psla$z[,,2]),xlim=c(-2,2),main="",col=cols[2])
for(i in 3:nclim){
  lines(density(psla$z[,,i]),col=cols[i])
}
abline(v=psla$z[,,1],col=cols[1],lty=3)

plot(density(tau_p*exp(psla$z[,,2])),xlim=c(0,250),ylim=c(0,0.011),main="",col=cols[1])
for(i in 3:nclim){
  lines(density(tau_p*exp(psla$z[,,i])),col=cols[i])
}
abline(v=tau_p*exp(psla$z[,,1]),col=cols[1],lty=3)

dev.off()

