#############################################################################
### Simulations of population dynamics for finite area, using species-level #
### median values for all parameters except germination ones                #
#############################################################################

library(plyr)
library(reshape2)
library(parallel)
library(abind)
library(RColorBrewer)

### LOAD DATA 

source("Source/simulation_functions_stochastic.R")
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
	go = readRDS("Models/go_pars_tdistpois_naspecies_noerr_noGDD_loglik_BH_01Mar2017.rds"),
	gs = readRDS("Models/gnzhh_onhh_pars_medians_26Oct2015.rds"),
		# gs = g site level
		# source script: venable_Stan_GO_descriptive_gnzhh_onhh_26Oct2015
		# uses tau_s = 100
		# but tau actually irrelevant because all multiplicative?
	pr = readRDS("Models/pr_pars_yearhet_squared_pc_02Mar2016.rds"),
	rs = readRDS("Models/rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")
	)

### PARAMS FOR SENSITIVITY ANALYSES

# transformed climate range approx. -1 -> 1
tmrange <- c(-1,1)
tsrange <- c(-2,2)
pls <- pl
nsens <- 100
plasticity <- F

# assign median parameters

for(i in 1:length(pl)){
  for(j in 1:length(pl[[i]])){
    ndim <- length(dim(pl[[i]][[j]])) 
    if(ndim==1){
      pls[[i]][[j]] <- rep(median(pl[[i]][[j]]),length(pl[[i]][[j]]))
    }
    if(ndim>1){
      keepdim <- 2:ndim # remove first dim, which is always iter
      eachrep <- dim(pl[[i]][[j]])[1] # select iterdim
      pls[[i]][[j]][] <- rep(apply(pl[[i]][[j]],keepdim,median),each=eachrep)
    }
  # if is.null(ndim), do nothing
  }
}

# assign germination parameters

if(plasticity==F){
  Gsens <- data.frame(
    alpha_G=qlogis(seq(0.001,0.999,length.out=nsens^2)),
    beta_Gz=rep(0,nsens^2)
  )
}

if(plasticity==T){
  Gsens <- expand.grid(
    tau_mu = seq(tmrange[1],tmrange[2],length.out=100),
    tau_sd = exp(seq(tsrange[1],tsrange[2],length.out=100))
  )
  Gsens$alpha_G <- with(Gsens,godalpha_f(tau_mu,tau_sd))
  Gsens$beta_Gz <- with(Gsens,godbeta_f(tau_sd))
}

pls$go$alpha_G[] <- Gsens$alpha_G
pls$go$beta_Gz[] <- Gsens$beta_Gz

nplot <- 4
Gshow <- round(
  rep(seq(1,100,length.out=nplot),times=nplot)
  + rep(seq(0,nsens^2-nsens,length.out=nplot),each=nplot),
  0)
Gplot <- Gsens[Gshow,]

# pdf(paste0("Plots/Germination_functions_", format(Sys.Date(),"%d%b%Y"),".pdf"),
# 		width=7,height=7)
par(mfrow=c(nplot,nplot),mar=c(3,3,1,1),las=1,ann=F,bty="l")
for(i in 1:nrow(Gplot)){
  with(Gplot, curve(plogis(alpha_G[i]+beta_Gz[i]*x),
    xlim=c(-2,2),ylim=c(0,1)) )
  if(plasticity==T){
    legend("topleft",bty="n",
      legend=c(paste0("mu=",signif(Gplot$tau_mu[i],2)),
        paste0("sd=",signif(Gplot$tau_sd[i],2))
        )
      )
  }
}
# dev.off()
  
# Sims --------------------------------------------------------------------

# 1000 per 0.1m^2 - *what does this mean?*

maml <- as.list(c(0,0,mpam,0,mpam,mpam))
msdl <- as.list(c(0,1,1,mpsd,mpsd,0))
  # scaling mean log rainfall (zamo) only works because sign stays the same

### Actual data
# nclim <- length(maml)
# cpc <- 4 # CORES per CLIMATE
# ncores <- nclim*cpc
# mpos <- rep(1:nclim,each=cpc)
# 
# nstart <- rep(10000,nspecies)
# ni <- 250 # iterations PER CORE
# nt <- 50
# nj <- 22
# nk <- 10000

### Sensitivity analyses
nclim <- length(maml)
cpc <- 5 # CORES per CLIMATE
ncores <- nclim*cpc
mpos <- rep(1:nclim,each=cpc)

nstart <- rep(10000,nspecies)
ni <- 2000 # iterations PER CORE
nt <- 50
nj <- 22
nk <- 1000 # 10000

# ni and nk must be >1
nit <- ni*cpc

set.seed(1)
maxiter <- 10000 # max number of iterations in PARAMETERISATION
cpos <- rep(1:cpc,times=nclim)
cipos <- rep(1:cpc,each=ni)
itersetl <- split(1:(ni*cpc),cipos)
  # requires that ni < maxiter

simp <- function(l){
  lapply(l,function(x){
    signif(x,2)
  })
}
  
cnames_unique <- gsub("\\.","",paste0("mu",simp(maml),"_sd",simp(msdl)))
cnames_bycore <- paste0(rep(cnames_unique,each=cpc),"_s",rep(1:cpc,times=nclim))
cnames_merged <- paste(cnames_unique,collapse="_")

# Each core focuses on one climate
# Iterations for one climate can be split up over multiple cores
# (controlled by cpc)
system.time({
CL = makeCluster(ncores)
clusterExport(cl=CL, c("popsim","pls", # was pl
  "zamo","zsdo","wamo","wsdo",
  "mpos","maml","msdl","cpos",
  "nstart","ni","nt","nj","nk",
  "Tvalues","itersetl","cnames_bycore"
  )) 
parLapply(CL, 1:ncores, function(n){
	mam <- maml[[mpos[n]]]
	msd <- msdl[[mpos[n]]]
	popsim(pl=pls,ni=ni,nt=nt,nj=nj,nk=nk, # was pl
	  nstart=nstart,zam=zamo+mam*zsdo,zsd=zsdo*msd,
	  zwy=zwyo,wam=wamo+wam*wsdo,wsd=wsdo*msd,rho=0.82,
		Tvalues=Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
		iterset=itersetl[[cpos[n]]],
		savefile=paste0("s_",cnames_bycore[n]) # added "s_"
		)
	})
stopCluster(CL)
})

# Read back in ------------------------------------------------------------

### Small RAM read-in

# cnames_bycore_small <- cnames_bycore[grep("mu1_sd1",cnames_bycore)]
# ncores_small <- length(cnames_bycore_small)
# psl <- as.list(rep(NA,ncores_small))
# for(n in 1:ncores_small){
#   psl[[n]] <- readRDS(paste0(cnames_bycore_small[n],"_28Apr2016.rds"))
#   }
# names(psl) <- cnames_bycore_small

psl <- as.list(rep(NA,ncores))
for(n in 1:ncores){
  # psl[[n]] <- readRDS(paste0("Sims/",cnames_bycore[n],"_07Apr2017.rds"))
  psl[[n]] <- readRDS(paste0("Sims/s_",cnames_bycore[n],"_19Jun2017.rds"))
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
        else ast <- abind(ast,psl[[firstmpos[j]+(k-1)]][[i]],along=1)
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
psla$P <- ifelse(psla$ns>0,1,0)
  # Binary persistence variable

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
  paste0("nj=",nj),
  paste0("nk=",nk)
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
parplot(goi$alpha_G,log(psla$ns),expression(alpha[G]),expression(ln(N[s15])),t=15)
parplot(goi$beta_Gz,log(psla$ns),expression(beta[G]),expression(ln(N[s15])),t=15)

parplot(gois$alpha_G,log(psla$ns),expression(alpha[G]),expression(ln(N[s15])),t=15,type="n")
parplot(gois$beta_Gz,log(psla$ns),expression(beta[G]),expression(ln(N[s15])),t=15,type="n")

parplot(gois$alpha_G,log(psla$ns_p),expression(alpha[G]),expression(ln(N[s50])),t=50,type="n",ylim=c(6,12))
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

parplot(gois$alpha_G,psla$P,expression(alpha[G50]),expression(P[50]),t=50,type="n")
parplot(plogis(gois$alpha_G),psla$pP,"G",expression(P[50]),t=50,type="n")

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

with(psla,matplot(t(log(ns[1:5,,17,1])),type="l",lty=1))
  # not overcompensating DD because doesn't switch direction every year
  # always some variation because of random year terms in reproduction
with(psla,matplot(t(log(Y[1:5,,17,1])),type="l",lty=1))
abline(h=0,lty=3)
  # basically never falls below replacement rate
with(psla,matplot(t(log(Y[1:5,,17,5])),type="l",lty=1))
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

# Summary stats -----------------------------------------------------------

Yeps <- psla$Y
Yeps[Yeps==0] <- min(psla$Y[!is.nan(psla$Y)],na.rm=T)
pehe_Ysd <- apply(log(Yeps[,,14,2]),1,function(x) sd(x[!is.nan(x) & is.finite(x)],na.rm=T))
scba_Ysd <- apply(log(Yeps[,,19,2]),1,function(x) sd(x[!is.nan(x) & is.finite(x)],na.rm=T))
hist(pehe_Ysd,breaks=100)
hist(scba_Ysd,breaks=100)
  # about 3

pehe_Ymean <- apply(log(Yeps[,,14,2]),1,function(x) mean(x[!is.nan(x) & is.finite(x)],na.rm=T))
scba_Ymean <- apply(log(Yeps[,,19,2]),1,function(x) mean(x[!is.nan(x) & is.finite(x)],na.rm=T))
hist(pehe_Ymean,breaks=100)
hist(scba_Ymean,breaks=100)
  # about 2

  # calculated in current climate (2)
  # but overall magnitude similar for future climate (3)
  # but *both measures biased because ignore zeroes*
