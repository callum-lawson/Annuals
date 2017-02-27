################################################################################
### Model to describe G and O values as function of species, year, and site, ###
### comparing negative binomial and lognormal-Poisson models                 ###
################################################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
# setwd("/mnt/data/home/NIOO/calluml/Source/")
source("venable_figure_functions_10Mar2016.R")

csyp <- read.csv("csyp_15Jan2016.csv",header=T)
ssyp <- read.csv("ssyp_15Jan2016.csv",header=T)
msy <- read.csv("msy_15Jan2016.csv",header=T)

library(rstan)
# library(lme4)
library(parallel)
library(reshape2)
library(plyr)
library(matrixStats) # for colMedians

############
### DATA ###
############

exlsp <- with(subset(msy,year<2001),
  names(which(tapply(totlive,species,sum,na.rm=T)<5))
)
  # species to be omitteds

tau_s <- 100 # Density per 10 cm square

spvals <- levels(msy$species)
keepvars <- c("species","year","plot","area")

gdata <- subset(csyp,is.na(ngerm)==F,select=c("species","year","plot","ngerm","area"))
odata <- subset(ssyp,is.na(totlive)==F
  & (
    ((species %in% exlsp)==T & year>=2001)
    |
    ((species %in% exlsp)==F & year>=1990)
    ),
  select=c("species","year","plot","totlive","area_o")
  )

gkeepnames <- names(gdata) %in% c("ngerm") == F
okeepnames <- names(odata) %in% c("totlive","area_o") == F

names(gdata)[gkeepnames] <- paste(names(gdata)[gkeepnames],"g",sep="_")
names(odata)[okeepnames] <- paste(names(odata)[okeepnames],"o",sep="_")

gdata$year_g <- as.factor(gdata$year_g)
odata$year_o <- as.factor(odata$year_o)

gdata$spyear_g <- with(gdata, as.numeric(factor(species_g:year_g)))
odata$spyear_o <- with(odata, as.numeric(factor(species_o:year_o)))
	# these have different numbers, so will have to match later

gdata$spsite_g <- with(gdata, as.numeric(factor(species_g:plot_g)))
odata$spsite_o <- with(odata, as.numeric(factor(species_o:plot_o)))

gdata$species_g <- match(gdata$species_g,spvals)
odata$species_o <- match(odata$species_o,spvals)

gdata$larea_g <- log(gdata$area_g * tau_s) # tau_s adjusts units from m2
odata$larea_o <- log(odata$area_o * tau_s)
	# fill in areas for unsurveyed qudrats

#gdata$ispos_g <- gdata$ngerm>0
#odata$ispos_o <- odata$totlive>0

#gdatap <- subset(gdata, ispos_g)
#odatap <- subset(odata, ispos_o)
   # DOESN'T FIT

glist <- as.list(gdata[,names(gdata) %in% c("year_g","plot_g")==F])
olist <- as.list(odata[,names(odata) %in% c("year_o","plot_o")==F])

dat_list <- c(glist,olist)

dat_list$L <- length(spvals)
dat_list$I_g <- nrow(gdata)
dat_list$I_o <- nrow(odata)
dat_list$J_g <- length(unique(gdata$spyear_g))
dat_list$J_o <- length(unique(odata$spyear_o))
dat_list$K_g <- length(unique(gdata$spsite_g))
dat_list$K_o <- length(unique(odata$spsite_o))

dat_list$sp_bysite_g <- with(gdata, species_g[match(unique(spsite_g),spsite_g)])
dat_list$sp_bysite_o <- with(odata, species_o[match(unique(spsite_o),spsite_o)])
dat_list$sp_byyear_g <- with(gdata, species_g[match(unique(spyear_g),spyear_g)])
dat_list$sp_byyear_o <- with(odata, species_o[match(unique(spyear_o),spyear_o)])

#dat_list$I_gp <- sum(gdata$ngerm_pos)
#dat_list$I_op <- sum(odata$totlive_pos)

############
### STAN ###
############

gmod_nzhh <- "
	data {

		int<lower=0> L;			// n species

		int<lower=0> I_g;			// n obs in g
		int<lower=0> J_g; 		// n years	
		int<lower=0> K_g;			// n sites in g
   	int<lower=0> ngerm[I_g];
   	vector[I_g] larea_g;	
		int<lower=0,upper=L> species_g[I_g];
		int<lower=0,upper=J_g> spyear_g[I_g];
		int<lower=0,upper=K_g> spsite_g[I_g];
    int<lower=0,upper=L> sp_byyear_g[J_g];
    int<lower=0,upper=L> sp_bysite_g[K_g];

		}		

	parameters {

		real alpha_mu;
		real<lower=0> sig_alpha_g;
		real sig_y_alpha_g;
		real<lower=0> sig_y_gam_g;
    real sig_s_alpha_g;
		real<lower=0> sig_s_gam_g;
		real alpha_psi_g;
		real<lower=0> sig_psi_g;
		real<lower=0> phi_g;

    vector[L] alpha_g;
    vector<lower=0>[L] sig_y_g;
    vector<lower=0>[L] sig_s_g;
    vector[J_g] eps_y_g;
 		vector[K_g] eps_s_g;
 		vector[L] psi_g;

		}
		
	transformed parameters {

		vector[L] theta_g;
		vector[I_g] loglam_g;

    for(i in 1:L){
      theta_g[i] <- inv_logit(psi_g[i]);
      }

		for (i in 1:I_g){
			loglam_g[i] <- larea_g[i] 
				+ alpha_g[species_g[i]]
				+ eps_y_g[spyear_g[i]]
				+ eps_s_g[spsite_g[i]];
			} 
		
		}

	model {

		alpha_g ~ normal(alpha_mu,sig_alpha_g);
    sig_y_g ~ lognormal(sig_y_alpha_g,sig_y_gam_g);
    sig_s_g ~ lognormal(sig_s_alpha_g,sig_s_gam_g);
    for(i in 1:J_g){
	    eps_y_g[i] ~ normal(0,sig_y_g[sp_byyear_g[i]]);
      }    
    for(i in 1:K_g){
	    eps_s_g[i] ~ normal(0,sig_s_g[sp_bysite_g[i]]);
      }
		psi_g ~ normal(alpha_psi_g,sig_psi_g);

    for(i in 1:I_g){
      if (ngerm[i] == 0){
        increment_log_prob(
          log_sum_exp(
            bernoulli_log(1,theta_g[species_g[i]]),
            bernoulli_log(0,theta_g[species_g[i]]) + neg_binomial_2_log_log(ngerm[i],loglam_g[i],phi_g)
            )
          );
        }
      else{
        increment_log_prob(
          bernoulli_log(0,theta_g[species_g[i]]) + neg_binomial_2_log_log(ngerm[i],loglam_g[i],phi_g)
          );
        }
      }

		}
	"

omod_nhh <- "
  data {

    int<lower=0> L;			// n species

    int<lower=0> I_o;			// n obs in o
    int<lower=0> J_o; 		// n years	
    int<lower=0> K_o;			// n sites in o
    int<lower=0> totlive[I_o];
    vector[I_o] larea_o;	
    int<lower=0,upper=L> species_o[I_o];
    int<lower=0,upper=J_o> spyear_o[I_o];
    int<lower=0,upper=K_o> spsite_o[I_o];
    int<lower=0,upper=L> sp_byyear_o[J_o];
    int<lower=0,upper=L> sp_bysite_o[K_o];

  }		

  parameters {

    real alpha_mu;
    real<lower=0> sig_alpha_o;
    real sig_y_alpha_o;
    real<lower=0> sig_y_gam_o;
    real sig_s_alpha_o;
    real<lower=0> sig_s_gam_o;

    vector[L] alpha_o;
    vector<lower=0>[L] sig_y_o;
    vector<lower=0>[L] sig_s_o;
    vector[J_o] eps_y_o;
    vector[K_o] eps_s_o;
    real<lower=0> phi_o;

  }

  transformed parameters {

    vector[I_o] loglam_o;

    for (i in 1:I_o){
      loglam_o[i] <- larea_o[i] 
      + alpha_o[species_o[i]]
      + eps_y_o[spyear_o[i]]
      + eps_s_o[spsite_o[i]];
      } 

  }

  model {

    alpha_o ~ normal(alpha_mu,sig_alpha_o);
    sig_y_o ~ lognormal(sig_y_alpha_o,sig_y_gam_o);
    sig_s_o ~ lognormal(sig_s_alpha_o,sig_s_gam_o);
    for(i in 1:J_o){
      eps_y_o[i] ~ normal(0,sig_y_o[sp_byyear_o[i]]);
      }    
    for(i in 1:K_o){
      eps_s_o[i] ~ normal(0,sig_s_o[sp_bysite_o[i]]);
      }

    totlive ~ neg_binomial_2_log(loglam_o,phi_o);

    }
"

### INITIATE MODEL

# setwd("D:/Users/calluml/Desktop/Stan/")
# setwd("C:/Users/Callum/Desktop/Stan/")
setwd("/mnt/data/home/NIOO/calluml/Output/")

gfit_nzhh <- stan(model_code=gmod_nzhh,data=dat_list,chains=0)
ofit_nhh <- stan(model_code=omod_nhh,data=dat_list,chains=0,warmup=10,iter=20)

### WINDOWS

fit_list <- list(gfit_nzhh=gfit_nzhh,ofit_nhh=ofit_nhh)
# fit_list <- list(gfit_nzhh=NA,ofit_nhh=NA)

nmods <- length(fit_list)
cpm <- 10 # chains per mod
fit_ind <- rep(1:nmods,each=cpm)
chain_ind <- rep(1:cpm,times=nmods)
nchains <- nmods*cpm
mod_names <- names(fit_list)

system.time({
CL = makeCluster(nchains, outfile=paste0("g_o_fit_nzhh_",format(Sys.Date(),"%d%b%Y"),".log"))
clusterExport(cl=CL, c("dat_list","fit_list","fit_ind","chain_ind","mod_names")) 
g_o_list <- parLapply(CL, 1:nchains, fun=function(i) {  # number of chains
	require(rstan)
 	stan(fit=fit_list[[fit_ind[i]]], 
		data=dat_list, 
		chains=1, 
		warmup=500, 	
		iter=1500,
		pars=c("alpha_mu"), # arbitrary
		chain_id=i,
		sample_file=paste0(mod_names[fit_ind[i]],"_chain",chain_ind[i],"_",format(Sys.Date(),"%d%b%Y"),".csv")
		)
	})
stopCluster(CL)
})

### READ BACK IN

CL = makeCluster(nchains) # load from earlier
go_sflist <- list()
clusterExport(cl=CL, c("go_sflist","fit_ind","chain_ind","mod_names")) 

go_sflist <- parLapply(CL, 1:nchains, function(i){
  require(rstan)
  read_stan_csv(paste0(mod_names[fit_ind[i]],"_chain",chain_ind[i],"_28Mar2016.csv"))
  })

gfit_nzhh <- sflist2stanfit(go_sflist[1:10])
ofit_nhh <- sflist2stanfit(go_sflist[11:20])

print(gfit_nzhh,pars="lp__")
print(ofit_nhh,pars="lp__")

traceplot(gfit_nzhh,ask=T)
traceplot(ofit_nhh,ask=T)
print(gfit_nzhh,pars=c("alpha_mu","sig_alpha_g","sig_s_alpha_g","sig_s_gam_g","sig_y_g","sig_s_g","phi_g","theta_g"))
print(ofit_nhh,pars=c("alpha_mu","sig_alpha_o","sig_s_alpha_o","sig_s_gam_o","sig_y_o","sig_s_o","phi_o","theta_o"))

### EXTRACT AND SAVE

gfit_nzhh_pars <- extract(gfit_nzhh)
ofit_nhh_pars <- extract(ofit_nhh)

lam_g <- exp(colMedians(gfit_nzhh_pars$loglam_g))
lam_o <- exp(colMedians(ofit_nhh_pars$loglam_o))

loglam_g_marg <- with(gfit_nzhh_pars, colMedians(loglam_g - eps_s_g[,gdata$spsite_g]))
loglam_o_marg <- with(ofit_nhh_pars, colMedians(loglam_o - eps_s_o[,odata$spsite_o]))

eps_y_g <- colMedians(gfit_nzhh_pars$eps_y_g)
eps_y_o <- colMedians(ofit_nhh_pars$eps_y_o)
eps_s_g <- colMedians(gfit_nzhh_pars$eps_s_g)
eps_s_o <- colMedians(ofit_nhh_pars$eps_s_o)

### HETEROSCEDASTICITY

nspecies <- length(unique(msy$species))
nsim <- 10^4

# YEARS

eps_y_g <- gfit_nzhh_pars$eps_y_g
eps_y_species_g <- gdata$species_g[match(unique(gdata$spyear_g),gdata$spyear_g)]
eps_y_g_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_y_g_sd[i,] <- tapply(eps_y_g[i,],eps_y_species_g,sd)
  }

eps_y_o <- ofit_nhh_pars$eps_y_o
eps_y_species_o <- odata$species_o[match(unique(odata$spyear_o),odata$spyear_o)]
eps_y_o_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_y_o_sd[i,] <- tapply(eps_y_o[i,],eps_y_species_o,sd)
  }

boxplot(eps_y_g_sd,range=0)
boxplot(eps_y_o_sd,range=0)

# SITES

eps_s_g <- gfit_nzhh_pars$eps_s_g
eps_s_species_g <- gdata$species_g[match(unique(gdata$spsite_g),gdata$spsite_g)]
eps_s_g_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_s_g_sd[i,] <- tapply(eps_s_g[i,],eps_s_species_g,sd)
  }

eps_s_o <- ofit_nhh_pars$eps_s_o
eps_s_species_o <- odata$species_o[match(unique(odata$spsite_o),odata$spsite_o)]
eps_s_o_sd <- matrix(nr=nsim,nc=nspecies)
for(i in 1:nsim){
  eps_s_o_sd[i,] <- tapply(eps_s_o[i,],eps_s_species_o,sd)
  }

boxplot(eps_s_g_sd,range=0)
boxplot(eps_s_o_sd,range=0)

### SAVE PARAMS

#g_keep <- c("alpha_g","sig_y_g","sig_s_g","phi_g","theta_g")
o_keep <- c("alpha_o","sig_y_o","sig_s_o","phi_o")

#gparl <- gfit_nzhh_pars[names(gfit_nzhh_pars) %in% g_keep]
oparl <- ofit_nhh_pars[names(ofit_nhh_pars) %in% o_keep]
goparl <- c(
    #gparl,
    oparl,
    list(
      #lam_g=lam_g,
      lam_o=lam_o,
      #loglam_g_marg=loglam_g_marg,
      loglam_o_marg=loglam_o_marg,
      #eps_y_g=eps_y_g,
      eps_y_o=eps_y_o,
      #eps_s_g=eps_s_g,
      eps_s_o=eps_s_o
      )
    )

saveRDS(goparl,paste0("onhh_pars_medians_naspecies_",format(Sys.Date(),"%d%b%Y"),".rds"))

#saveRDS(goparl,paste0("gnzhh_onhh_pars_medians_",format(Sys.Date(),"%d%b%Y"),".rds"))

