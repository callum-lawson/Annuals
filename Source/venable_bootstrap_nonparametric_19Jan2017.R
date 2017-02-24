################################################
### Non-parametric bootstrap of total counts ###
################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"
# setwd("/mnt/data/home/NIOO/calluml/Source/")

library(plyr)
library(reshape2)
library(rstan)
library(parallel)
library(lme4)
library(matrixStats)

# setwd(paste0(maindir,"Analyses/Venable"))
# source("venable_figure_functions_23Apr2016.R")

msy <- read.csv("msy_seedests_18Jan2017.csv",header=T)
  # seedests not actually necessary?
csyp <- read.csv("csyp_seedests_19Jan2017.csv",header=T) 
ssyp <- read.csv("ssyp_15Jan2016.csv",header=T) 

####################
### PREPARE DATA ###
####################

spvals <- levels(msy$species)

idvars <- c("species","year","plot")
csyp_s <- subset(csyp,
  is.na(germd)==F,
  select=c(idvars,"ngerm","area")
  )
  # for ng

exlsp <- with(subset(msy,year<2001),
  names(which(tapply(totlive,species,sum,na.rm=T)<5))
  )
  # species to be omitted due to weird seed data
ssyp_s <- subset(ssyp,
  is.na(totlive)==F & (
    ((species %in% exlsp)==T & year>=2001)
    |
    ((species %in% exlsp)==F & year>=1990)
    ),
  select=c(idvars,"totlive","area_o")
  )
  # for no
nsyp_s <- subset(csyp,
  is.na(nsd)==F,
  select=c(idvars,"seed_hat","area")
  )
  # for nn
nsyp_s$seed_hat <- round(nsyp_s$seed_hat,0)
  # expected values were given as decimals

names(csyp_s) <- names(ssyp_s) <- names(nsyp_s) <- c(idvars,"count","area")
csyp_s <- droplevels(csyp_s)
ssyp_s <- droplevels(ssyp_s)
nsyp_s <- droplevels(nsyp_s)
  # needed?

tau_d <- 100 	# adjustment for density
tau_p <- 100 	# adjustment for rainfall
tau_s <- 100 # Density per 10 cm square

##################################################
### BOOTSTRAP SAMPLING OF TOTALS FROM RAW DATA ###
##################################################

bootstrap <- function(indat,nboot){
  nrow <- nrow(indat)
  lplot <- levels(indat$plot)
  nplot <- nlevels(indat$plot)
  plotvals <- sample(lplot,size=nplot*nboot,replace=T)
  
  rowvals <- lapply(plotvals,function(x) which(indat$plot==x))
  sampdat <- indat[unlist(rowvals),]
  nrowvals <- sapply(rowvals,length)
  sumnrowvals <- tapply(nrowvals,rep(1:nboot,each=nplot),sum)
  sampdat$boot <- rep(1:nboot,sumnrowvals)
  
  totdat <- ddply(sampdat,.(boot,species,year),summarise,
    area = sum(area),
    count = sum(count)
    )
  return(totdat)
  }
  
nboot <- 100
cbd <- bootstrap(csyp_s,nboot=nboot)
sbd <- bootstrap(ssyp_s,nboot=nboot)
nbd <- bootstrap(nsyp_s,nboot=nboot)
  # sums over plots

bbd <- merge(cbd,sbd,by=c("boot","species","year"),sort=F)
  # dropping cdb for which sbd species-year combos are missing
  # (weird seed data; see above)
bbd$count <- with(bbd, round(count.x + count.y*(area.x/area.y),0))
  # scaling seed counts to match germ counts
bbd$area <- with(bbd, area.x)
  # using germ area (same as seed area once scaled)
bbd <- subset(bbd,select=c("species","year","boot","count","area"))

### DATA EXPLORATION

hist(bbd$count,breaks=10^4,xlim=c(0,500))
hist(nbd$count,breaks=10^5,xlim=c(0,500))
table(nbd$count>0)
with(nbd, hist(count[count>0],breaks=10^5,xlim=c(0,5000)))
with(nbd, hist(log(count[count>0]),breaks=10^5))

with(bbd,barplot(tapply(count,list(year,species),mean)))
with(indat,barplot(tapply(count/area,list(year,species),mean)))
with(subset(indat,year %in% 1990:1999),barplot(tapply(count/area,list(year,species),mean)))

with(indat,tapply(count/area,list(year,species),mean))
sumtab <- with(indat,log(tapply(count/area,list(year,species),sum)))
barplot(t(sumtab))

spmeans <- with(indat, tapply(count/area,species,mean))
resid <- with(indat,(count/area) / spmeans[species])
tapply(resid,indat$species,var)
  # variances between years too different for different species?
  
with(subset(bbd, species==spvals[20] & year=="2014"),hist(count/area,breaks=10^4,xlim=c(0,500)))

#######################################
### MODEL-FITTING (T-DISTRIBUTIONS) ###
#######################################

### PREPARE DATA

# bigsp <- with(bbd, names(sort(tapply(count,species,sum),decreasing=T)))
# indat <- droplevels(subset(bbd,species %in% bigsp[1:10]))
  
# oldsp <- names(which(!is.na(with(subset(bbd,year=="1990"), tapply(count,species,sum)))))
# indat <- droplevels(subset(bbd,species %in% oldsp))

# indat <- subset(bbd,species=="scba")
# indat <- droplevels(indat)

indat <- bbd

indat$year <- as.factor(indat$year)
indat$spyear <- with(indat, as.numeric(factor(species:year)))
  # these have different numbers, so will have to match later
indat$species <- as.numeric(as.factor(indat$species)) # match(indat$species,spvals)
  # to make sure that species goes 1:new species number
indat$larea <- log(indat$area * tau_d) 
  # tau_d = tau_s = tau_p

dl <- as.list(indat)
dl$L <- length(unique(indat$species))
dl$I <- nrow(indat)
dl$J <- length(unique(indat$spyear))

### PARAMETERISE MODEL

tdist_mod <- "
  data {
    int<lower=1> L;     // n species (for nu)
    int<lower=1> I;     // n obs
    int<lower=1> J;     // n spyears
    int<lower=0> count[I];
    vector[I] larea;
    int<lower=1,upper=L> species[I];
    int<lower=1,upper=J> spyear[I];
    }		

  parameters {
    real alpha_mu;
    real<lower=0> sig_alpha;
    vector[J] alpha;

    real sig_o_mu;
    real<lower=0> sig_o_sig;
		vector<lower=0>[L] sig_o;
		vector<lower=0>[L] nu_o;
		vector[I] eps_o;
    }
  
  transformed parameters {
    vector[I] loglam;
		vector<lower=0>[I] nu_o_vec;
 		vector<lower=0>[I] sig_o_vec;	

    for (i in 1:I){
      nu_o_vec[i] = nu_o[species[i]];
      sig_o_vec[i] = sig_o[species[i]];
      loglam[i] = larea[i] + alpha[spyear[i]] + eps_o[i];
      }
    }
  
  model {
    alpha ~ normal(alpha_mu,sig_alpha);

    sig_o ~ lognormal(sig_o_mu,sig_o_sig);
    nu_o ~ gamma(2,0.1);
  	eps_o ~ student_t(nu_o_vec,0,sig_o_vec);

    sig_alpha ~ cauchy(0,2.5);
    sig_o_sig ~ cauchy(0,2.5);

    count ~ poisson_log(loglam);
    }
"

# tdist_mod <- "
#   data {
#     int<lower=1> I;     // n obs
#     int<lower=1> J;     // n spyears
#     int<lower=0> count[I];
#     vector[I] larea;
#     int<lower=0,upper=J> spyear[I];
#     }		
#   
#   parameters {
#     real alpha_mu;
#     real<lower=0> sig_y;
#     vector[J] eps_y;
#     }
#   
#   transformed parameters {
#     vector[I] loglam;
#     
#     for (i in 1:I){
#       loglam[i] = larea[i] + alpha_mu + eps_y[spyear[i]];
#       }
#     }
#   
#   model {
#     eps_y ~ normal(0,sig_y);
#     sig_y ~ cauchy(0,2.5);
# 
#     count ~ poisson_log(loglam);
#     }
# "

tdist_fit <- stan(model_code=tdist_mod,data=dl,chains=0)

# setwd("/mnt/data/home/NIOO/calluml/Output/")

# WINDOWS
nchains <- 10

system.time({
  CL = makeCluster(nchains, outfile=paste0("tdist_descriptive_",format(Sys.Date(),"%d%b%Y"),".log"))
  clusterExport(cl=CL, c("dl","tdist_fit")) 
  tdist_sflist <- parLapply(CL, 1:nchains, fun = function(cid) {
    require(rstan)
    stan(fit=tdist_fit, 
      data=dl, 
      chains=1, 
      warmup=2000,
      iter=3000,
      pars=c("alpha_mu"),
      sample_file=paste0("tdist_descriptive_fits_chain",cid,"_",format(Sys.Date(),"%d%b%Y"),".csv"),
      chain_id=cid
      )
    })
  stopCluster(CL)
  })

### READ IN

# setwd("/mnt/data/home/NIOO/calluml/Output/")

# finchains <- c(4,6,8,10)
finchains <- 1:10

nchains <- length(finchains)
CL = makeCluster(nchains)
tdist_sflist <- list()
clusterExport(cl=CL, c("tdist_sflist","finchains")) 

tdist_sflist <- parLapply(CL, 1:nchains, function(i){
  require(rstan)
  read_stan_csv(paste0("tdist_descriptive_fits_chain",finchains[i],"_19Jan2017.csv"))
  })

tdist_fit <- sflist2stanfit(tdist_sflist) 

traceplot(tdist_fit)
mypars <- c("sig_o","nu_o")
traceplot(tdist_fit,pars=mypars)
print(tdist_fit,pars=mypars)

pairs(tdist_fit,pars=c("alpha_mu","sig_alpha","sig_y"))

tpars <- rstan::extract(tdist_fit) # extract function also from dplyr
sig_o_t <- with(tpars,colMedians(sig_o))
nu_o_t <- with(tpars,colMedians(nu_o))

### New seed densities

nbdagg <- ddply(nbd,.(species,year),summarise,
  totcount = sum(count)
  )
  # only using years for which total count was >0
library(dplyr)
nbda <- left_join(nbd,nbdagg,by=c("species","year"))
nbdap <- subset(nbda,totcount>0)

stddevs <- function(obj){
  unlist(lapply(VarCorr(obj), function(m) sqrt(diag(m))))
	}

ssds <- sapply(split(nbdap,nbdap$species),function(d){
    tsm_r <- glmer(count ~ offset(log(area*tau_s))
      + (1|year/boot),
      family="poisson",
      data=d
    )
    obssd <- stddevs(tsm_r)["boot:year.(Intercept)"]
    return(obssd)
  })

### Save outputs

obserrors <- data.frame(
  species = spvals,
  sig_b = sig_o_t, 
  nu_b = nu_o_t,
  sig_n = ssds,
  row.names=NULL
  )

write.csv(obserrors,
  paste0("observation_error_byspecies_",format(Sys.Date(),"%d%b%Y"),".csv"),
  row.names=F
  )

