####################################################
### Population-level models for Pr, Rs, G, and O ###
####################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

library(plyr)
library(reshape2)
library(lme4)
library(rstan)
library(parallel)
library(matrixStats)

setwd(paste0(maindir,"Analyses/Venable"))
# setwd("/mnt/data/home/NIOO/calluml/Source/")
source("venable_figure_functions_23Apr2016.R")

msy <- read.csv("msy_seedests_18Jan2017.csv",header=T)
# no sds attached yet
Tvalues <- read.csv("Tvalues_20Jun2015.csv",header=T)
obserr <- read.csv("observation_error_byspecies_23Jan2017.csv",header=T)

nspecies <- nlevels(msy$species)
spvals <- levels(msy$species)

################
### ADD DATA ###
################

tau_s <- 100		# adjustment for seed density
tau_p <- 100		# adjustment for rainfall

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

# with(msy90,cbind(
#   sig_b=tapply(sig_b,species,median),
#   sig_b=tapply(nu_b,species,median),
#   sig_n=tapply(sig_n,species,median)
#   ))

go_mod <- "
functions {

  real BH(real n_in, real m1, real m2, real T3) {
    real n_out;
    n_out = n_in * exp(-m1*T3) / ( 1 + (m2/m1)*(1-exp(-m1*T3))*n_in );
    return n_out;
    }
  
  real RICKER(real n_in, real m1, real m2, real T3) {
    real n_out;
    n_out = n_in * exp(-(m1+m2*n_in)*T3);
    return n_out;
    }

  }

data {

  int<lower=0> I;
  int<lower=0> J;
  int<lower=0> L;
  vector<lower=0>[L] J0; 	// olsd start year
  int<lower=0> N;	// n positive totseed
  int<lower=0> M;	// n >0 totboth
  
  int<lower=1,upper=L> species[I];
  int<lower=1,upper=J> year[I];
  int<lower=1,upper=I> mpos[M];
  
  vector[I] tprcp;
  int<lower=0> totboth[I];
  int<lower=0> totgerm_p[M];
  int<lower=0> totboth_p[M];
  int<lower=0,upper=1> isseed[I];	
  
  vector[I] larea_g;	
  
  real<lower=0> T1;
  real<lower=0> T2;
  real<lower=0> T3;
  
  int<lower=0,upper=N> npos[I];	// zero because placeholder
  vector[N] lseeddens;

  vector<lower=0>[N] sig_n; // including as data

  }		

parameters {

  real alpha_G_mu;
  real<lower=0> sig_G_alpha;
  vector[L] alpha_G;
  
  real beta_Gz_mu;
  real<lower=0> sig_Gz_beta;
  vector[L] beta_Gz;
  
  real<lower=0>  m1_mu; // for gamma prior specification
  real<lower=0> sig_m_alpha;   	
  vector[L] alpha_m;
  
  real beta_m_mu;
  real<lower=0> sig_m_beta;   		
  vector[L] beta_m;
  
  real lbd0_mu;
  real<lower=0> sig_lbd0;   		
  vector[L] lbd0; 
  // log starting buried density
  
  real lnd_mu;
  real<lower=0> sig_lnd;   		
  vector[N] lnd; //  N not I
  
  real sig_b_mu;
  real<lower=0> sig_b_sig;
  vector<lower=0>[L] sig_b;
  vector<lower=0>[L] nu_b;
  vector[I] eps_b;

  // real<lower=0> sig_n; 
    // including as one param for all species

  }

transformed parameters {

  real alpha_m_mu;
  
  vector<lower=0,upper=1>[I] G;
  vector<lower=0>[I] m1;
  vector<lower=0>[I] m2;
  
  vector[I] lbd;
  vector[I] lgd;
  vector[I] lod;
  vector[I] lgod;
  vector<lower=0>[I] nd;
  vector<lower=0>[I] odb;
  vector<lower=0>[I] ndb;
  
  vector[M] G_adj_lgt;
  vector[I] loglam_b;
  
  vector<lower=0>[I] sig_b_vec;	
  vector<lower=0>[I] nu_b_vec;
  
  alpha_m_mu = log(m1_mu);
  
  for (i in 1:I){
  
    G[i] = inv_logit(alpha_G[species[i]] + beta_Gz[species[i]]*tprcp[i]);
    m1[i] = exp(alpha_m[species[i]]);
    m2[i] = exp(beta_m[species[i]]);
    
    if(year[i]==J0[species[i]]){
      lbd[i] = lbd0[species[i]];
      lod[i] = lbd[i] + log(1-G[i]) - m1[i]*T1; 
      }
  
    if(year[i]>J0[species[i]]){
      lbd[i] = log(odb[i-1] + ndb[i-1]);
      // i-1 always same species, previous year
      lod[i] = lbd[i] + log(1-G[i]) - m1[i-1]*T1;
      // prev year's mortality (so that same mort for one cycle)
      }	
  
    lgd[i] = lbd[i] + log(G[i]);
    
    lgod[i] = log( exp(lgd[i]) + exp(lod[i]) );
    loglam_b[i] = lgod[i] + larea_g[i] + eps_b[i];
    
    odb[i] = exp(lod[i] - m1[i]*(T2+T3));
    
    // NEW SEEDS
    
    if(isseed[i]==0) {	
      nd[i] = 0;
      ndb[i] = 0;
      // no lnd needed because ipos means that never used
      }
    else {
      nd[i] = exp(lnd[npos[i]]);
      ndb[i] = BH(nd[i],m1[i],m2[i],T3);
      }
      // densities all in m2*tau_s
    
    nu_b_vec[i] = nu_b[species[i]];
    sig_b_vec[i] = sig_b[species[i]];
    
    }
  
  for(i in 1:M){
    G_adj_lgt[i] = logit(G[mpos[i]]) + m1[mpos[i]]*T1; 
    // multiplying log-odds by 1/S[1] = exp(m1*T1)
    // totboth already adjusted for area
    // so totgerm and totlive counts comparable
    }

  }

model {

  // PRIORS
  
  alpha_G_mu ~ normal(0,5);
  m1_mu ~ gamma(1,1);
  // moles et al 2003, moriuchi et al 2000?
  
  alpha_G ~ normal(alpha_G_mu,sig_G_alpha);
  beta_Gz ~ normal(beta_Gz_mu,sig_Gz_beta);
  
  alpha_m ~ normal(alpha_m_mu,sig_m_alpha);
  beta_m ~ normal(beta_m_mu,sig_m_beta);
  
  lbd0 ~ normal(lbd0_mu,sig_lbd0);
  lnd ~ normal(lnd_mu,sig_lnd);
  
  sig_b ~ lognormal(sig_b_mu,sig_b_sig);
  nu_b ~ gamma(2,0.1);
  eps_b ~ student_t(nu_b_vec,0,sig_b_vec);
  // eps_b ~ normal(0,sig_b_vec);

  sig_G_alpha ~ cauchy(0,2.5);
  sig_Gz_beta ~ cauchy(0,2.5);
  sig_m_alpha ~ cauchy(0,2.5);
  sig_m_beta ~ cauchy(0,2.5);
  sig_b ~ cauchy(0,2.5);
  sig_b_sig ~ cauchy(0,2.5);
  
  // LIKELIHOOD
  
  totgerm_p ~ binomial_logit(totboth_p,G_adj_lgt);
  totboth ~ poisson_log(loglam_b);
  lseeddens ~ normal(lnd,sig_n); // msy pos
  
  }
"

# setwd("/mnt/data/home/NIOO/calluml/Output/")

go_fit <- stan(model_code=go_mod,data=dl,chains=0)

# WINDOWS
nchains <- 10

system.time({
CL = makeCluster(nchains, outfile=paste0("gofit_poplevel_binomialG_tdistbpar_normndat_",format(Sys.Date(),"%d%b%Y"),".log"))
clusterExport(cl=CL, c("dl","go_fit")) 
go_sflist <- parLapply(CL, 1:nchains, fun = function(cid) {  # number of chains
  require(rstan)
  stan(fit=go_fit, 
       data=dl, 
       chains=1, 
       warmup=2000,
       iter=3000,
       pars=c("alpha_G_mu","beta_Gz_mu"),
       sample_file=paste0("go_fits_chain",cid,"_poplevel_binomialG_tdistbpar_normndat_",format(Sys.Date(),"%d%b%Y"),".csv"),
       chain_id=cid
  )
})
stopCluster(CL)
})

### READ IN

# setwd("/mnt/data/home/NIOO/calluml/Output/")

finchains <- 1:nchains

nchains <- length(finchains)
CL = makeCluster(nchains)
go_sflist <- list()
clusterExport(cl=CL, c("go_sflist","finchains")) 

go_sflist <- parLapply(CL, 1:nchains, function(i){
  require(rstan)
  read_stan_csv(paste0("go_fits_chain",finchains[i],"_poplevel_binomialG_tdistbpar_normndat_13Feb2017.csv"))
  })

go_fit <- sflist2stanfit(go_sflist) 

traceplot(go_fit,pars=c("alpha_G_mu","beta_Gz_mu","alpha_m_mu","beta_m_mu"),nr=2,nc=3)
traceplot(go_fit,pars=c("alpha_G"),nr=6,nc=4)
traceplot(go_fit,pars=c("beta_Gz"),nr=6,nc=4)
traceplot(go_fit,pars=c("alpha_m"),nr=6,nc=4)
traceplot(go_fit,pars=c("beta_m"),nr=6,nc=4)

print(go_fit)
print(go_fit,pars=c("alpha_G","beta_Gz","alpha_m","beta_m"))
traceplot(go_fit,pars=c("alpha_G"))

gopars <- extract(go_fit)
gonames <- names(gopars)
goshow <- sapply(gopars,function(x) (nrow(msy90) %in% dim(x))==F)
print(go_fit,pars=gonames[goshow])
traceplot(go_fit,c("sig_G_alpha","sig_Gz_beta","sig_Gd_beta","sig_m_alpha","sig_m_beta","sig_g_sig","sig_o_sig"))
  # sig_Gd_beta quite autocorrelated

# Likelihood

loglik_all <- with(gopars,rowMedians(log_lik))
( loglik <- sum(loglik_all) )
( deviance <- -2 * loglik ) 

######################
### CALCULATE PARS ###
######################

gopars <- extract(go_fit)

lbd <- with(gopars,colMedians(lbd))
lgd <- with(gopars,colMedians(lgd))
lod <- with(gopars,colMedians(lod))
lnd <- with(gopars,colMedians(lnd))
nd <- with(gopars,colMedians(nd))
odb <- with(gopars,colMedians(odb))
ndb <- with(gopars,colMedians(ndb))

G <- with(gopars,colMedians(G))
m1 <- with(gopars,colMedians(m1))
m2 <- with(gopars,colMedians(m2))

loglam_g <- with(gopars,colMedians(loglam_g))
loglam_o <- with(gopars,colMedians(loglam_o))

loglam_g_marg <- with(gopars,colMedians(loglam_g-eps_g))
loglam_o_marg <- with(gopars,colMedians(loglam_o-eps_o))

sig_G_eps <-  with(gopars,median(sig_G_eps))
sig_m_eps <-  with(gopars,median(sig_m_eps))

eps_m <-  with(gopars,colMedians(eps_m))
eps_G <-  with(gopars,colMedians(eps_G))
eps_g <-  with(gopars,colMedians(eps_g))
eps_o <-  with(gopars,colMedians(eps_o))

nu_g <- with(gopars,colMedians(nu_g))
nu_o <- with(gopars,colMedians(nu_o))

sig_g <- with(gopars,colMedians(sig_g))
sig_o <- with(gopars,colMedians(sig_o))

sdll_g <- with(gopars,colMedians(sdll_g))
sdll_o <- with(gopars,colMedians(sdll_o))

### CALIBRATION

par(mfrow=c(1,2))
hist(eps_G,breaks=50,prob=T)
curve(dnorm(x,0,sig_G_eps),col="red",add=T)
hist(eps_m,breaks=50,prob=T)
curve(dnorm(x,0,sig_m_eps),col="red",add=T)

tgreen <- rgb(red=0,green=1,blue=0,alpha=0.3,maxColorValue=1)
tblue <- rgb(red=0,green=0,blue=1,alpha=0.3,maxColorValue=1)

pdf(paste0("obserror_hists_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=plotwidth,height=plotheight)
par(mfrow=c(6,4))
for(i in 1:nspecies){
  cursp <- msy90$species==spvals[i]
  mybreaks <- hist(eps_g[cursp],breaks=50,prob=T,main=spvals[i])
  hist(rt(10^5,df=nu_g)*msy90$gsd[cursp][1],add=T,col=tgreen,prob=T,breaks=c(-Inf,mybreaks$breaks,Inf))
    }
par(mfrow=c(6,4))
for(i in 1:nspecies){
  cursp <- msy90$species==spvals[i]
  mybreaks <- hist(eps_o[cursp],breaks=50,prob=T,main=spvals[i])
  hist(rt(10^5,df=nu_o)*msy90$osd[cursp][1],add=T,col=tblue,prob=T,breaks=c(-Inf,mybreaks$breaks,Inf))
}
dev.off()

msy90$olsdhat <- with(msy90,totlive/area_o)
msy90$Ghat <- with(msy90,germdhat/(germdhat+olsdhat))

alphaG_med <- with(gopars,colMedians(alpha_G))
beta_Gz_med <- with(gopars,colMedians(beta_Gz))
beta_Gd_med <- with(gopars,colMedians(beta_Gd))

G_med <- plogis(
  alphaG_med[dl$species] 
  + beta_Gz_med[dl$species]*dl$tprcp 
  + beta_Gd_med[dl$species]*(log(msy90$germdhat+msy90$olsdhat)-log(tau_s))
  )

par(mfrow=c(1,1))
plot(qlogis(msy90$Ghat)~qlogis(G_med))
abline(0,1,col="red")
lines(supsmu(qlogis(G_med),qlogis(msy90$Ghat)),col="blue")
table(msy90$Ghat>G_med)/nrow(msy90)

plot(log(msy90$Ghat/G_med) ~ msy90$year)
abline(h=0,col="red",lty=2)

par(mfrow=c(1,2))
plot(tapply(msy90$gsd,msy90$species,median),sig_g)
abline(0,1,col="red")
plot(tapply(msy90$osd,msy90$species,median),sig_o)
abline(0,1,col="red")

par(mfrow=c(1,1))
plot(qlogis(msy90$Ghat)~qlogis(G))
abline(0,1,col="red")
table(msy90$Ghat>G)

with(msy90,plot(qlogis(tapply(Ghat,species,median,na.rm=T))~qlogis(tapply(G,species,median))))
abline(0,1,col="red")

par(mfrow=c(1,2))
plot(log(msy90$germdhat)~I(lgd+log(tau_s)))
abline(0,1,col="red")
plot(log(msy90$olsdhat)~I(lod+log(tau_s)))
abline(0,1,col="red")

npar <- nrow(msy90)
nsim <- 10^3
ntot <- npar*nsim
gsimmat <- osimmat <- area_g_mat <- area_o_mat <- matrix(NA,nr=npar,nc=nsim)
area_g_mat[] <- rep(msy90$area_g,nsim)
area_o_mat[] <- rep(msy90$area_o,nsim)
epsg <- rt(ntot,df=nu_g)*rep(msy90$gsd,nsim)
epso <- rt(ntot,df=nu_o)*rep(msy90$osd,nsim)

gsimmat[] <- rpois(ntot,
  lambda=exp(rep(lgd+log(area_g_mat[,1]*tau_s),times=nsim)+epsg)
  ) / area_g_mat
osimmat[] <- rpois(ntot,
  lambda=exp(rep(lod+log(area_o_mat[,1]*tau_s),times=nsim)+epso)
) / area_o_mat

# make distribution of residuals too

pdf(paste0("hists_lnpois_BH_naspecies_diffnu_GDD_loglik_",format(Sys.Date(),"%d%b%Y"),".pdf"),
  width=12,height=8)
par(mfrow=c(2,3),mar=c(4.5,4.5,1.5,1.5))
mybreaks <- hist(msy90$germdhat,breaks=50,prob=T,density=20,xlab=expression(N[g]),main="")
hist(exp(lgd+log(tau_s)),add=T,col=tgreen,prob=T,breaks=c(0,mybreaks$breaks,Inf))
mybreaks <- hist(msy90$germdhat,breaks=50,prob=T,density=20,xlab=expression(N[g]),main="")
hist(gsimmat,add=T,col=tgreen,prob=T,breaks=c(0,mybreaks$breaks,Inf))
plot(log(msy90$germdhat)~I(lgd+log(tau_s)))
abline(0,1,col="red")

mybreaks <- hist(msy90$olsdhat,breaks=50,prob=T,density=20,xlab=expression(N[o]),main="")
hist(exp(lod+log(tau_s)),add=T,col=tblue,prob=T,breaks=c(0,mybreaks$breaks,Inf))
mybreaks <- hist(msy90$olsdhat,breaks=50,prob=T,density=20,xlab=expression(N[o]),main="")
hist(osimmat,add=T,col=tblue,prob=T,breaks=c(0,mybreaks$breaks,Inf))
plot(log(msy90$olsdhat)~I(lod+log(tau_s)))
abline(0,1,col="red")
dev.off()

par(mfrow=c(1,1))
curve(dgamma(x,shape=2,rate=0.1),xlim=c(0,20))

#################
### SAVE PARS ###
#################

# lapply(gopars,dim)
# lapply(goparl,dim)

gopars <- extract(go_fit)

goparl <- as.list(rep(NA,length(gopars)))
names(goparl) <- names(gopars)
for(i in 1:length(gopars)){
  if(nrow(msy90) %in% dim(gopars[[i]])){
    goparl[[i]] <- colMedians(gopars[[i]])
  }
  else{
    goparl[[i]] <- gopars[[i]]
  }
}

goparl$loglam_b_marg <- with(gopars,colMedians(loglam_b-eps_b))

setwd("/mnt/data/home/NIOO/calluml/Source/")
saveRDS(goparl,
  paste0("go_pars_binomialG_tdistbpar_normndat_",format(Sys.Date(),"%d%b%Y"),".rds")
  )

