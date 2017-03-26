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
seriesplot(psla$pY,"pY",yname="Pr(Y>0)",quantiles=F)

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

### Pairs plots

pairplot("popchange_pairs_vitalratechange_newenv",rca,2)
pairplot("popchange_pairs_vitalratechange_consenv",cca,2)
pairplot("popchange_pairs_allclimscenarios_",rcaa,3)
pairplot("popchange_pairs_allclimscenarios_alltraits",rcata,3,w=21,h=21)

### Climate graphs

pdf(paste0("Plots/zdists",format(Sys.Date(),"%d%b%Y"),".pdf"),width=4.5,height=4.5)

plot(density(psla$z[,,2]),xlim=c(-2,2),main="",col=cols[2])
for(i in 3:nclim){
  lines(density(psla$z[,,i]),col=cols[i])
  }
abline(v=psla$z[,,1],col="black",lty=3)

plot(density(tau_p*exp(psla$z[,,2])),xlim=c(0,250),ylim=c(0,0.011),main="",col=cols[1])
for(i in 3:nclim){
  lines(density(tau_p*exp(psla$z[,,i])),col=cols[i])
}
abline(v=tau_p*exp(psla$z[,,1]),col="black",lty=3)

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

# Predictions of yearly dynamics ------------------------------------------

### FINISH ANNOTATIONS!

### Y every year, X every year and every species
# y[i,t,j,m] ~ x[i,t,j,m]

# r ~ z
# r ~ G
# r ~ ns
# Y ~ gdens
# Pr(Y>0) ~ gdens

### Y every year, X every year, same for all species
# y[i,t,j,m] ~ x[i,t,m]

# Y ~ z
# Ye ~ z

### Y one year, X one value from sims
# y[i,j,m] ~ x[i,j,m]

# ns_median ~ G_median

### Y one value, X one param
# y[i,j,m] ~ x[i,j]

# ns_median ~ alpha_G

# species index gets ignored

### y:
# - [1000, 50, 22, 5]
# - [1000, 22, 5]

### x : 
# - [1000, 50, 22, 5]
# - [1000, 50, 5]
# - [1000, 22]
# - [1000, 5]

pslm <- pli <- list() # rename?
pslm$ns <- apply(psla$ns,c(1,3,4),median)
pli$alpha_G <- pl$go$alpha_G[unlist(itersetl),]

parplot(pli$alpha_G,log(pslm$ns),t=NULL,xname=expression(alpha[G]),yname=expression(ln(ns[t=15])))

x <- pli$alpha_G
y <- log(pslm$ns)

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
