#########################################################
### Estimating per-capita reproduction (Y) parameters ###
#########################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
# setwd("/mnt/data/home/NIOO/calluml/Source/")

cdpos <- read.csv("cdpos_15Jan2016.csv",header=T)
csyp <- read.csv("csyp_15Jan2016.csv",header=T)
msy <- read.csv("msy_15Jan2016.csv",header=T)

#keepsp <- c("boin","evmu","scba","pere","plpa")
#keeppl <- sample(levels(csyp$plot),size=10)
#cdpos <- droplevels(subset(cdpos,year>2006 & species %in% keepsp & plot %in% keeppl))
#csyp <- droplevels(subset(csyp,year>2006 & species %in% keepsp & plot %in% keeppl))
#msy <- droplevels(subset(msy,year>2006 & species %in% keepsp))

# library(lme4)
library(rstan)
library(parallel)
library(plyr)
library(matrixStats)

####################
### PREPARE DATA ###
####################

tau_d <- 100 	# adjustment for density
tau_p <- 100 	# adjustment for rainfall

prdat <- subset(csyp, ngerm>0 & is.na(nmeas)==F & is.na(germd/tau_d)==F)
rsdat <- cdpos

spvals <- levels(msy$species)
yrvals <- unique(msy$year)

dat_list <- list()

dat_list$L <- length(spvals)
# dat_list$J <- length(yrvals)

msy$spyear <- with(msy, as.numeric(factor(species:as.factor(year))))
msy$species_num <- as.numeric(msy$species)
msy_merge <- subset(msy,select=c("species","year","spyear"))
prdat <- merge(prdat,msy_merge,by=c("species","year"),all.x=T)
rsdat <- merge(rsdat,msy_merge,by=c("species","year"),all.x=T)
dat_list$J <- length(unique(msy$spyear))

prdat$spsite <- with(prdat, as.numeric(factor(species:plot)))
prdat$species <- match(prdat$species,spvals)
prdat$year <- match(prdat$year,yrvals)

rsdat$spsite <- with(rsdat, as.numeric(factor(species:plot)))
rsdat$species <- match(rsdat$species,spvals)
rsdat$year <- match(rsdat$year,yrvals)

dat_list$I_p <- nrow(prdat)
dat_list$K_p <- max(prdat$spsite)

dat_list$I_r <- nrow(rsdat)
dat_list$K_r <- max(rsdat$spsite)
  # could combine Ks if useful

jointnames <- c("species","spyear","spsite")
names(prdat)[names(prdat) %in% jointnames] <- paste0(names(prdat)[names(prdat) %in% jointnames],"_p")
names(rsdat)[names(rsdat) %in% jointnames] <- paste0(names(rsdat)[names(rsdat) %in% jointnames],"_r")

dat_list <- c(dat_list,
    as.list(subset(prdat,select=c(nfert,nmeas,species_p,spyear_p,spsite_p))),
    as.list(subset(rsdat,select=c(seedsint,species_r,spyear_r,spsite_r)))
    )

dat_list$sp_byyear <- with(msy, species_num[match(unique(spyear),spyear)])
dat_list$sp_bysite_r <- with(rsdat, species_r[match(unique(spsite_r),spsite_r)])

dat_list$x_p <- cbind(
  intfake=rep(1,dat_list$I_p),
  ltprcp=log(prdat$prcp/tau_p),
	ltprcp2=(log(prdat$prcp/tau_p))^2,
  ltgermd=log(prdat$germd/tau_d)
  )

dat_list$x_r <- cbind(
  intfake=rep(1,dat_list$I_r),
  ltprcp=log(rsdat$prcp/tau_p),
	ltprcp2=(log(rsdat$prcp/tau_p))^2,
  ltgermd=log(rsdat$germd/tau_d)
  )

dat_list$M_p <- ncol(dat_list$x_p) # 4 columns in the x matrices
dat_list$M_r <- ncol(dat_list$x_r) # 4 columns in the x matrices
dat_list$u <- cbind(rep(1,dat_list$L))

# dat_list$offset <- c(0,log(tau_p),log(tau_d)) 
# offset multipliers
# first term is intercept -> no offset

### MODEL

prmod <- "
	data {

   	int<lower=1> I_p;   // n data
   	int<lower=1> J; 	  // n years
   	int<lower=1> K_p;   // n species-site combos
   	int<lower=1> L; 	  // n species
   	int<lower=1> M_p; 	  // n betas
   	int<lower=0> nfert[I_p];
   	int<lower=1> nmeas[I_p];
		int<lower=1,upper=J> spyear_p[I_p];
		int<lower=1,upper=K_p> spsite_p[I_p];
		int<lower=1,upper=L> species_p[I_p];
    int<lower=0,upper=L> sp_byyear[J];
   	matrix[I_p,M_p] x_p;
		matrix[L,1] u;			// dummy variable

		}		

	parameters {

		matrix[M_p,L] z_p;
		cholesky_factor_corr[M_p] L_Omega_p;
		vector<lower=0>[M_p] tau_p; 	// prior scale
		matrix[1,M_p] beta_mu_p; 	    
      // combines with u (1s for each species) to create matrix of (identical) beta_mus for each group 

    real sig_y_alpha_p;
		real<lower=0> sig_y_gam_p;

    vector<lower=0>[L] sig_y_p;
    real<lower=0> sig_s_p; 
   	real<lower=0> sig_o_p; 
	
  	vector[J] eps_y_p;
    vector[K_p] eps_s_p;
    vector[I_p] eps_o_p;
	  }

	transformed parameters{

		matrix[L,M_p] beta_p;
		vector[I_p] pi;   // logit scale!

		// COEFFICIENTS
		beta_p <- u * beta_mu_p + (diag_pre_multiply(tau_p,L_Omega_p) * z_p)';
			// u*gamma = M x L matrix (or L x M) = beta_mus for each species

		// LINEAR PREDICTOR
		for (i in 1:I_p)
			pi[i] <- x_p[i] * beta_p[species_p[i]]' // [L,M] or [M,L] matrix
				+ eps_y_p[spyear_p[i]]
				+ eps_s_p[spsite_p[i]]
				+ eps_o_p[i];

		}

	model {

		// PRIORS
		tau_p ~ cauchy(0,2.5);
		to_vector(z_p) ~ normal(0,1);
		L_Omega_p ~ lkj_corr_cholesky(2); // could change to (1)
		to_vector(beta_mu_p) ~ normal(0,5);
			// Change tau/Omega/z/beta_mu priors? Too informative?

    sig_y_p ~ lognormal(sig_y_alpha_p,sig_y_gam_p);
    for(i in 1:J){
      eps_y_p[i] ~ normal(0,sig_y_p[sp_byyear[i]]);
      }

		eps_s_p ~ normal(0,sig_s_p);
		eps_o_p ~ normal(0,sig_o_p);
	
		// LIKELIHOOD
    nfert ~ binomial_logit(nmeas,pi);

	  }

  generated quantities {

		matrix[M_p,M_p] Omega_p;
    vector[I_p] log_lik_p;

		Omega_p <- L_Omega_p * L_Omega_p';

    for (i in 1:I_p){
      log_lik_p[i] <- binomial_logit_log(nfert[i],nmeas[i],pi[i]);
      }

		}
	"

rsmod <- "
	data{

   	int<lower=1> I_r; 			// n data
   	int<lower=1> J; 			// n years
   	int<lower=1> K_r; 			// n species-site combos
   	int<lower=1> L; 			// n species
   	int<lower=1> M_r; 			// n betas
		int<lower=1> seedsint[I_r];
		int<lower=1,upper=J> spyear_r[I_r];
		int<lower=1,upper=K_r> spsite_r[I_r];
		int<lower=1,upper=L> species_r[I_r];
    int<lower=0,upper=L> sp_byyear[J];
    int<lower=0,upper=L> sp_bysite_r[K_r];
   	matrix[I_r,M_r] x_r;
		matrix[L,1] u;				// dummy variable

		}

	parameters{

	  matrix[M_r,L] z_r;
		cholesky_factor_corr[M_r] L_Omega_r;
		vector<lower=0>[M_r] tau_r; 	  // prior scale
		matrix[1,M_r] beta_mu_r;        // combines with u (1s for each species) to create matrix of (identical) beta_mus for each group 
	
    real sig_y_alpha_r;
		real<lower=0> sig_y_gam_r;
	  // real sig_s_alpha_r;
		// real<lower=0> sig_s_gam_r;

    vector<lower=0>[L] sig_y_r;
    // vector<lower=0>[L] sig_s_r;
    real<lower=0> sig_s_r;
		vector[J] eps_y_r;
    vector[K_r] eps_s_r;
		real<lower=0> phi_r;

		}

	transformed parameters{

		matrix[L,M_r] beta_r;
		vector[I_r] loglam_r;

		// COEFFICIENTS
		beta_r <- u * beta_mu_r + (diag_pre_multiply(tau_r,L_Omega_r) * z_r)';
			// u*gamma = M x L matrix (or L x M) = beta_mus for each species

  	// LIKELIHOOD
		for(i in 1:I_r){
			loglam_r[i] <- exp(
				x_r[i] * beta_r[species_r[i]]' // [L,M] or [M,L] matrix
					+ eps_y_r[spyear_r[i]] 
					+ eps_s_r[spsite_r[i]]
				);
				// put exp here if not using neg_binomial_2_log
		  }

		}

	model{

		// PRIORS
		tau_r ~ cauchy(0,2.5);
		to_vector(z_r) ~ normal(0,1);
		L_Omega_r ~ lkj_corr_cholesky(4); // could change to (1)
		to_vector(beta_mu_r) ~ normal(0,5);
			// Change tau/Omega/beta_mu priors? Too informative?

		sig_y_gam_r ~ cauchy(0,2.5);
		// sig_s_gam_r ~ cauchy(0,2.5);
    sig_y_r ~ lognormal(sig_y_alpha_r,sig_y_gam_r);
    // sig_s_r ~ lognormal(sig_s_alpha_r,sig_s_gam_r);
    sig_s_r ~ cauchy(0,2.5);

    for(i in 1:J){
    	eps_y_r[i] ~ normal(0,sig_y_r[sp_byyear[i]]);
    	}
		//for(i in 1:K_r){
	 	//	eps_s_r[i] ~ normal(0,sig_s_r[sp_bysite_r[i]]);
    //	}
		eps_s_r ~ normal(0,sig_s_r);

		for(i in 1:I_r){
			seedsint[i] ~ neg_binomial_2(loglam_r[i],phi_r) T[1,];
			}

		// seedsint ~ neg_binomial_2_log(loglam_r,phi_r);
		// removed truncation T[1,] (and exp around poisson function)

	  }

  generated quantities {

		matrix[M_r,M_r] Omega_r;
    vector[I_r] log_lik_r;

		Omega_r <- L_Omega_r * L_Omega_r';

    for (i in 1:I_r){
      log_lik_r[i] <- neg_binomial_2_log(seedsint[i],loglam_r[i],phi_r)
        / ( 1 - neg_binomial_2_log(0,loglam_r[i],phi_r) );
      }

		}
	"

### COMPILE MODEL

prfit <- stan(model_code=prmod,data=dat_list,chains=0)
rsfit <- stan(model_code=rsmod,data=dat_list,chains=0)
rsfit <- stan(model_code=rsmod,data=dat_list,chains=1,warmup=10,iter=20,init=inits[1])

### SET UP PARALLEL

# setwd("D:/Users/calluml/Desktop/Stan/")
setwd("/mnt/data/home/NIOO/calluml/Output/")

# fit_list <- list(prfit=prfit,rsfit=rsfit)
# par_list <- list("beta_mu_p","beta_mu_r")
# apptext <- c("yearhet_nositehet_squared_pc","yearhet_sitehet_linear_pc_trunc")
# warm_list <- c(500,25)
# iter_list <- c(1500,75)

# PR

fit_list <- list(prfit=prfit)
par_list <- list("beta_mu_p")
apptext <- c("yearhet_nositehet_squared_loglik")
warm_list <- 500
iter_list <- 1500

nmods <- length(fit_list)
cpm <- 10 # chains per mod # earlier 10
fit_ind <- rep(1:nmods,each=cpm)
chain_ind <- rep(1:cpm,times=nmods)
nchains <- nmods*cpm
mod_names <- paste(names(fit_list),apptext,sep="_")

# RS

fit_list <- list(rsfit=rsfit)
par_list <- list("beta_mu_r")
apptext <- c("yearhet_nositehet_squared_trunc_loglik")
warm_list <- 500
iter_list <- 1500

nmods <- length(fit_list)
cpm <- 12 # chains per mod # earlier 10
fit_ind <- rep(1:nmods,each=cpm)
chain_ind <- rep(1:cpm,times=nmods)
nchains <- nmods*cpm
mod_names <- paste(names(fit_list),apptext,sep="_")

inits <- readRDS("rs_inits_sq_b_03Mar2016.rds")
inits <- inits[1:cpm]
inits_sq <- lapply(inits,function(x){
	x[names(x) %in% c("z_r","L_Omega_r","tau_r","beta_mu_r","beta_r","Omega_r")==F]
	})

system.time({
  CL = makeCluster(nchains, outfile=paste0(paste(mod_names,collapse="_"),"_",format(Sys.Date(),"%d%b%Y"),".log"))
  clusterExport(cl=CL, c("dat_list","fit_list","par_list","fit_ind","chain_ind","mod_names","warm_list","iter_list","inits")) 
  p_r_list <- parLapply(CL, 1:nchains, fun=function(i) {  # number of chains
    require(rstan)
    stan(fit=fit_list[[fit_ind[i]]], 
         data=dat_list, 
         chains=1, 
         warmup=warm_list[fit_ind[i]], 	
         iter=iter_list[fit_ind[i]],
         pars=par_list[[fit_ind[i]]],
         chain_id=i,
    		 init=inits[i], # OMITTED FOR PR!
         sample_file=paste0(mod_names[fit_ind[i]],"_chain",chain_ind[i],"_",format(Sys.Date(),"%d%b%Y"),".csv")
    )
  })
  stopCluster(CL)
})
	# was warmup=500, iter=1000

### READ BACK IN

# PR

pchains <- 10
CL = makeCluster(pchains)
p_sflist <- list()
clusterExport(cl=CL, c("p_sflist","pchains")) 
p_sflist <- parLapply(CL, 1:pchains, function(i){
  require(rstan)
  if(i %in% 1:pchains){
    read_stan_csv(paste0("prfit_yearhet_nositehet_squared_loglik_chain",i,"_19May2016.csv"))
    }
  })
stopCluster(CL)

# RS

finchains <- 1:10 # cutting off spare chains 11,12
maxchains <- max(finchains)   # !!! only runs loop up to here !!!
CL = makeCluster(nchains) # load from earlier
r_sflist <- list()
clusterExport(cl=CL, c("r_sflist","fit_ind","chain_ind","mod_names","finchains")) 

r_sflist <- parLapply(CL, 1:maxchains, function(i){
  require(rstan)
  if(i %in% finchains){
    read_stan_csv(paste0(mod_names[fit_ind[i]],"_chain",chain_ind[i],"_19May2016.csv"))
    }
  })
stopCluster(CL)

# TRACEPLOTS

prfit <- sflist2stanfit(p_sflist)
traceplot(prfit)

# rsfit <- sflist2stanfit(p_r_sflist[rsfin])
rsfit <- sflist2stanfit(r_sflist)  # 03 Mar
traceplot(rsfit)
traceplot(rsfit,pars="beta_mu_r")
traceplot(rsfit,pars="phi_r")

traceplot(sflist2stanfit(r_sflist[finchains[10]]))
rsconv <- sflist2stanfit(r_sflist[finchains[c(5,10)]])
print(rsfit,pars=("beta_mu_r"))

### INITIAL VALUES

# setwd("/mnt/data/home/NIOO/calluml/Source/")
# rs <- readRDS("rs_pars_siteyearhet_trunc_18Nov2015.rds")
# str(rs)
# with(rs, apply(beta_mu_r,c(2,3),mean))

initextract <- function(m,n,endits=T){
	initlist <- list()
	for(i in 1:n){
		niter <- m@stan_args[[1]]$iter - m@stan_args[[1]]$warmup
		if(endits==F){
			itnum <- sample(x=1:niter,size=1)
			set.seed(i)
			}
		if(endits==T){
			itnum <- niter - n + i # won't work without permute=F
			}
		allpars <- extract(m)
		initlist[[i]] <- lapply(allpars,function(p){
			pdim <- length(dim(p))
			if(pdim==1) return(p[itnum])
			if(pdim==2) return(p[itnum,])
			if(pdim==3) return(p[itnum,,])
			if(pdim==4) return(p[itnum,,,])
			})
		}
	return(initlist)
	}

inits <- initextract(m=rsconv,n=12,endits=F)
inits <- lapply(inits,function(x){
	pos <- which(names(x)=="beta_mu_r")
	x[[pos]] <- matrix(x[[pos]],nr=1,nc=4) # changed to 4 for squared
	return(x)
	})
	# turning beta_mu_r into row vector
	
print(prfit,pars="lp__")
print(rsfit,pars="lp__")

###############################
### EXTRACT AND SAVE PARAMS ###
###############################

prpars <- extract(prfit,permuted=T) 
rspars <- extract(rsfit,permuted=T) 

# params with uncertainty
# predictions and marginal predictions without uncertainty

### LIKELIHOODS

log_lik_p <- with(prpars,rowSums(log_lik_p))
median(log_lik_p)
log_lik_r <- with(rspars,rowSums(log_lik_r))
median(log_lik_r)

prsq <- prpars$beta_mu_p[,,1]
table(prsq<0)

### ALL PARAMS

prparl <- rsparl <- list()

### CALCULATE LOGLAM FOR RS
# didn't save loglam so need to calculate for each row of cd

nbtmean <- function(mu,phi){
	mu / ( 1 - (phi/(mu+phi))^phi )
	}
	# calculates mean of truncated negative binomial distribution

niter <- length(rspars$phi_r)
loglam_r <- mu_r_all <- matrix(NA,nr=niter,nc=dat_list$I_r)
for(i in 1:niter){
	loglam_r[i,] <- rowSums(dat_list$x_r * rspars$beta_r[i,dat_list$species_r,])
			+ rspars$eps_y_r[i,dat_list$spyear_r] 
			+ rspars$eps_s_r[i,dat_list$spsite_r]
	mu_r_all[i,] <- nbtmean(exp(loglam_r[i,]),rspars$phi_r[i])
	}

### CALCULATE MEDIANS

library(matrixStats)

pi <- colMedians(prpars$pi)
pi_marg <- with(prpars, colMedians(pi - (eps_s_p[,prdat$spsite] + eps_o_p)))
eps_y_p <- with(prpars,colMedians(eps_y_p))
eps_s_p <- with(prpars,colMedians(eps_s_p))
eps_o_p <- with(prpars,colMedians(eps_o_p))
Omega_p <- apply(prpars$Omega_p,c(2,3),median)

mu_r <- colMedians(mu_r_all)
# loglam_r_marg <- with(rspars, colMedians(loglam_r - eps_s_r[,rsdat$spsite]))
eps_y_r <- with(rspars,colMedians(eps_y_r))
eps_s_r <- with(rspars,colMedians(eps_s_r))
Omega_r <- apply(rspars$Omega_r,c(2,3),median)
  # marginal predictions *still include year effects*

liknb <- function(x,...){
  dnbinom(x,
    mu=rspars$loglam_r,
    size=rep(rspars$phi_r,times=dat_list$I_r),
    ...
  )
}

yobsmat <- matrix(rep(dat_list$seedsint,each=10^4),nr=10^4,nc=dat_list$I_r)

log_lik_r_full <- matrix(NA,nr=10^4,ncol=dat_list$I_r)
log_lik_r_full[] <- 
  liknb(yobsmat,log=T) - log( 1 - liknb(rep(0,dat_list$I_r*10^4),log=F) ) 

log_lik_r <- with(rspars,rowSums(log_lik_r_full))
median(log_lik_r)

library(loo)
loo1 <- loo(log_lik_r_full)
waic1 <- waic(log_lik_r_full)

cbind(
  rspars$loglam_r[1,log_lik_r_full[1,]==-Inf],
  yobsmat[1,log_lik_r_full[1,]==-Inf]
  )

curve(dnbinom(x,mu=rspars$loglam_r[1],size=rspars$phi_r[1]),xlim=c(0,100))

### ASSESS SPECIES HETEROGENEITY

nspecies <- length(unique(msy$species))
nsim <- 10^4

# SITE

eps_s_species_p <- prdat$species_p[match(unique(prdat$spsite_p),prdat$spsite_p)]
eps_s_p_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_s_p_sd[i,] <- tapply(prpars$eps_s_p[i,],eps_s_species_p,sd)
  }

eps_s_species_r <- rsdat$species_r[match(unique(rsdat$spsite_r),rsdat$spsite_r)]
eps_s_r_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_s_r_sd[i,] <- tapply(rspars$eps_s_r[i,],eps_s_species_r,sd)
  }

boxplot(eps_s_p_sd,range=0)
abline(h=median(prpars$sig_s_p),col="red",lty=3)
boxplot(eps_s_r_sd,range=0)
abline(h=median(rspars$sig_s_r),col="red",lty=3)

# YEAR
msy$species_num <- as.numeric(msy$species)

eps_y_species_p <- msy$species_num[match(unique(prdat$spyear_p),msy$spyear)]
eps_y_p_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_y_p_sd[i,] <- tapply(prpars$eps_y_p[i,unique(prdat$spyear_p)],eps_y_species_p,sd)
  }

eps_y_species_r <- msy$species_num[match(unique(rsdat$spyear_r),msy$spyear)]
eps_y_r_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_y_r_sd[i,] <- tapply(rspars$eps_y_r[i,unique(rsdat$spyear_r)],eps_y_species_r,sd)
  }

boxplot(eps_y_p_sd,range=0)
abline(h=median(prpars$sig_y_p),col="red",lty=3)
boxplot(eps_y_r_sd,range=0)
abline(h=median(rspars$sig_y_r),col="red",lty=3)

### SAVE

pr_keep <- c("beta_mu_p","sig_y_p","sig_s_p","sig_o_p","beta_p")  
rs_keep <- c("beta_mu_r","sig_y_r","sig_s_r","phi_r","beta_r")   

prparl <- c(prpars[names(prpars) %in% pr_keep],
            list(pi=pi,
                 pi_marg=pi_marg,
                 eps_y_p=eps_y_p,
                 eps_s_p=eps_s_p,
                 eps_o_p=eps_o_p,
                 Omega_p=Omega_p
                 )
            )
rsparl <- c(rspars[names(rspars) %in% rs_keep],
            list(#loglam_r=loglam_r,
                 #loglam_r_marg=loglam_r_marg,
            		 mu_r=mu_r,
                 eps_y_r=eps_y_r,
                 eps_s_r=eps_s_r,
                 Omega_r=Omega_r
                 )
            )

saveRDS(prparl,paste0("pr_pars_yearhet_squared_pc_",format(Sys.Date(),"%d%b%Y"),".rds"))
saveRDS(rsparl,paste0("rs_pars_yearhet_squared_pc_trunc_",format(Sys.Date(),"%d%b%Y"),".rds"))
saveRDS(inits,paste0("rs_inits_sq_b_",format(Sys.Date(),"%d%b%Y"),".rds"))

