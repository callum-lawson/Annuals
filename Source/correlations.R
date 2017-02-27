#################################################################
### Look at among-species correlations across years and plots ###
#################################################################

library(plyr)
library(lme4)

source("Source/figure_functions.R")

csy <- read.csv("Output/csy_15Jan2016.csv",header=T)
msy <- read.csv("Output/msy_15Jan2016.csv",header=T)

csyp <- read.csv("Output/csyp_15Jan2016.csv",header=T)
ssyp <- read.csv("Output/ssyp_15Jan2016.csv",header=T)

sb <- read.csv("Output/sb_15Jan2016.csv",header=T)

# csyp <- read.csv("csyp_20Jun2015.csv",header=T)
# ssyp <- read.csv("ssyp_20Jun2015.csv",header=T)

#############################
### VARIANCE PARTITIONING ###
#############################

sb$obslevel <- 1:nrow(sb)
ssyp$obslevel <- 1:nrow(ssyp)
sb$area <- rep(pi*(0.054/2)^2,nrow(sb))	# diameter is 5.4cm
csyp$obslevel <- as.factor(1:nrow(csyp))

sb$spyear <- with(sb, factor(species:as.factor(year)))
csyp$spyear <- with(csyp, factor(species:as.factor(year)))

model <- glmer(nlive~offset(log(area*1000))+(1|species/plot)+(1|year)+(1|obslevel),family=poisson,data=sb)
model2 <- glmer(totlive~offset(log(area*1000))+(1|species/plot)+(1|year)+(1|obslevel),family=poisson,data=ssyp)

modelb <- glmer(nlive~offset(log(area*1000))+(1|species/plot)+(1|spyear)+(1|obslevel),family=poisson,data=sb)
	# didn't converge
model2b <- glmer(totlive~offset(log(area*1000))+(1|species/plot)+(1|spyear)+(1|obslevel),family=poisson,data=ssyp)

summary(model)
summary(model2)
overdisp_fun(model)
overdisp_fun(model2)

### each year then calculate means

gmod <- glmer(ngerm~offset(log(area*1000))+(1|species/plot)+(1|year)+(1|obslevel),
	family=poisson,data=csyp
	)
omod <- glmer(nlive~offset(log(area*1000))+(1|species/plot)+(1|year)+(1|obslevel),
	family=poisson,data=sb
	)
omod2 <-  glmer(totlive~offset(log(area*1000))+(1|species/plot)+(1|year)+(1|obslevel),
	family=poisson,data=ssyp
	)

summary(gmod)
summary(omod)
summary(omod2)

gmeans <- coef(gmod)$species[[1]]
omeans <- coef(omod)$species[[1]]
omeans2 <- coef(omod2)$species[[1]]
gvar <- sqrt(4.274 + 3.089)
ovar <- sqrt(5.5996 + 0.9191)
ovar2 <- sqrt(3.906 + 0.608)

n <- 10^4
	# 10^8 made no difference to estimates

	# could try without 2006

nspecies <- nlevels(csy$species)
gsum <- osum <- osum2 <- gcalc <- gcalc2 <- vector()

for(i in 1:nspecies){
	gsim <- rpois(n,lambda=rlnorm(n,meanlog=gmeans[i],sdlog=gvar))
	osim <- rpois(n,lambda=rlnorm(n,meanlog=omeans[i],sdlog=ovar))
	osim2 <- rpois(n,lambda=rlnorm(n,meanlog=omeans2[i],sdlog=ovar2))
	gsum[i] <- sum(gsim)
	osum[i] <- sum(osim)
	osum2[i] <- sum(osim2)
	}

gcalc <- gsum/(gsum+osum)
gcalc2 <- gsum/(gsum+osum2)

names(gcalc) <- names(gcalc2) <- levels(csy$species)

msy$G <- with(msy,germd/(germd+olsd))
gtrue <- with(msy,tapply(G,list(species,year),mean,na.rm=T))
gtruemarg <- apply(gtrue,1,mean,na.rm=T)

plot(gcalc,gcalc2,xlim=c(0,1))
	# estimates consistent regardless of level
plot(gcalc,gtruemarg,xlim=c(0,1))
plot(gcalc2,gtruemarg,xlim=c(0,1))

# Calibration

ng_true <- with(msy, tapply(germd/10^3,species,mean,na.rm=T))
plot(log(ng_true)~log(gsum/n))
abline(0,1,col="red",lty=2)
abline(lm(log(ng_true)~log(gsum/n)),col="blue",lty=3)

no_true <- with(msy, tapply(olsd/10^3,species,mean,na.rm=T))
plot(log(no_true)~log(osum/n))
abline(0,1,col="red",lty=2)
abline(lm(log(no_true)~log(osum/n)),col="blue",lty=3)

# GLMM ADMB

library(glmmADMB)
csyp$year <- as.factor(csyp$year)
csyp$stupid <- rnbinom(nrow(csyp),size=1,mu=1)
csyp$ltarea <- with(csyp,log(area*1000))

model1 <- glmmadmb(ngerm~(1|species/plot)+(1|year)+offset(ltarea),
	data=subset(csyp,is.na(ngerm)==F & is.na(area)==F),
	zeroInflation=FALSE,family="nbinom",link="log")

model1 <- glmmadmb(stupid~1,
	data=csyp,
	zeroInflation=FALSE,family="nbinom",link="log")
	# won't even start running

### SAMPLE SIZE SIMULATION

N <- 10^(1:6)
I <- length(N)
R <- 10^3

estmean <- matrix(NA,nr=R,nc=I)

for(i in 1:I){
	for(j in 1:R){
		estmean[j,i] <- mean(rpois(N,lambda=rlnorm(N,meanlog=0,sdlog=1)))
		}
	}

boxplot(estmean)
	# no obvious effect of sample size on estimates

############################
### CORRELATION FUNCTION ###
############################

paircor_ij <- function(i,j,veclist){
	if(i==j) return(NA)
	if(i!=j){
		fakeadd <- veclist[[i]] + veclist[[j]]
			# arbitrary, to check if both non-NA
		nfinite <- sum(is.na(fakeadd)==F)
		if(nfinite < 3) return(NA) 
		if(nfinite >= 3) cor.test(veclist[[i]],veclist[[j]])$estimate 
		}
	}
paircor <- Vectorize(paircor_ij,vectorize.args=list("i","j"))

nspecies <- nlevels(csy$species)

### experimenting

plot(log(germd)~log(olsd),data=msy)
plot(log((germd+olsd)/prevcsd)~log(prevcsd),data=msy,col=species,pch="+")
plot((germd)/(olsd+germd)~log(prevcsd),data=msy,col=species,pch="+")
plot((germd)/(olsd+germd)~log(olsd+germd),data=msy,col=species,pch="+")

nlevels(csyp$plot) *  nlevels(csy$species) *  length(unique(csy$year)) * 8000
	# all = 5 million
	# used to be 3.2 million, so not a big issue

with(subset(csyp,species=="scba"),plot(log((olsdbar+germd)/(prevnsd+prevolsdbar))~log(prevnsd+prevolsdbar))) 
with(subset(ssyp,species=="scba"),plot(log((olsd+germdbar)/(prevnsdbar+prevolsd))~log(prevnsdbar+prevolsd))) 

with(subset(ssyp,species=="scba"),plot(log(olsd/prevolsd)~log(prevnsdbar+prevolsd))) 
	# low prevolsd -> high prevgermd -> high prevnsd, so can't separate

hist(csyp$nsd[csyp$year=="1990"],breaks=1000)

plot(log(quantile(csyp$nsd,prob=seq(0,1,0.001),na.rm=T)),xlim=c(800,1000))
points(log(quantile(ssyp$olsd,prob=seq(0,1,0.001),na.rm=T)),col="red")

plot(log(quantile(csyp$nsdbar,prob=seq(0,1,0.001),na.rm=T)))
points(log(quantile(ssyp$olsdbar,prob=seq(0,1,0.001),na.rm=T)),col="red")

plot(log(nsd)~year,data=msy)
points(log(olsd)~year,data=msy,col="red")

par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(log(olsd)~log(prevolsd),data=ssyp)
with(ssyp,lines(supsmu(log(prevolsd),log(olsd)),col="red"))
with(ssyp,lines(supsmu(log(prevolsdbar),log(olsd)),col="blue"))
abline(0,1,col="gray",lty=2)
plot(log(nsd)~log(prevnsd),data=csyp)
with(csyp,lines(supsmu(log(prevnsd),log(nsd)),col="red"))
with(csyp,lines(supsmu(log(prevnsdbar),log(nsd)),col="blue"))
abline(0,1,col="gray",lty=2)

par(mfrow=c(1,2),mar=c(2,2,2,2))
plot(log(olsd)~log(prevolsd),data=subset(ssyp,species=="scba"))
with(subset(ssyp,species=="scba"),lines(supsmu(log(prevolsd),log(olsd)),col="red"))
with(subset(ssyp,species=="scba"),lines(supsmu(log(prevolsdbar),log(olsd)),col="blue"))
abline(0,1,col="gray",lty=2)

with(ssyp, cor.test(olsd,prevolsd))
with(ssyp, cor.test(olsd,prevolsdbar))

with(subset(ssyp,olsd>0 & prevolsd>0), cor.test( log(olsd),log(prevolsd) ) )
with(subset(ssyp,olsd>0 & prevolsdbar>0), cor.test( log(olsd),log(prevolsdbar) ) )
	# cor slightly better if using local for seed plots (but not much)

with(csyp, cor.test(nsd,prevnsd))
with(csyp, cor.test(nsd,prevnsdbar))

with(subset(csyp,nsd>0 & prevnsd>0), cor.test( log(nsd),log(prevnsd) ) )
with(subset(csyp,nsd>0 & prevnsdbar>0), cor.test( log(nsd),log(prevnsdbar) ) )
	# cor better if using local for germ plots

	# Suggests either seed global, germ local, or both global

	# (for all species and years pooled)
	# definitely quite strong correlation, but not huge relative to other noise sources
	# -> some movement to other sites...?
	# correlation for olsd stronger than nsd in both cases
	# clearly substantial sampling error because of regression to the mean

par(mfrow=c(1,3))
hist(msy$germd,breaks=100)
hist(msy$nsd,breaks=100)
hist(msy$olsd,breaks=100)
hist(log(msy$germd),breaks=100,xlim=c(0,10))
hist(log(msy$nsd),breaks=100,xlim=c(0,10))
hist(log(msy$olsd),breaks=100,xlim=c(0,10))

plot(log((germd+olsd)/(prevnsd+prevolsd))~log(prevnsd+prevolsd),data=msy)

nsd_mean <- with(csyp,tapply(nsd,list(species,year),mean,na.rm=T))
nsd_quant <- with(csyp,tapply(nsd,list(species,year),quantile,prob=0.8,na.rm=T))
plot(log(nsd_mean)~log(nsd_quant))

with(msy,hist(germd/(prevolsd+prevnsd),breaks=100))
	# a few have WAY more germinants than seeds form last year
with(msy,hist(germd/(germd+olsd),breaks=100))

with(msy,qlogis(mean(germd/(germd+olsd),na.rm=T)))

with(msy,hist(log(germd/(prevolsd+prevnsd)),breaks=100))
with(msy,hist(log(germd/(germd+olsd)),breaks=100))

par(mfrow=c(1,2))
with(csyp,hist(log(germd/(prevolsdbar+prevnsdbar)),breaks=100))
with(csyp,hist(log(germd/(prevolsdbar+germd)),breaks=100))

toobig <- subset(csyp,germd/(prevolsdbar+prevnsdbar)>1)
toobigtab <- with(toobig,table(species,year))
rowSums(toobigtab)
colSums(toobigtab)

toobig2 <- subset(msy,(germd+olsd)/(prevolsd+prevnsd)>1)
toobigtab2 <- with(toobig2,table(species))

with(csyp,table(species,year,is.na(pr)==F))

par(mfrow=c(1,2))
with(csyp,hist(log(germd/(prevolsdbar+prevnsdbar)),breaks=100))
with(csyp,hist(log(germd/(prevolsdbar+prevnsd)),breaks=100))

qlog1 <- with(msy,qlogis(germd/(prevolsd+prevnsd)))
qlog2 <- with(msy,qlogis(germd/(olsd+germd)))
par(mfrow=c(1,2))
hist(qlog1[is.finite(qlog1) & is.na(qlog1)==F],breaks=100)
hist(qlog2[is.finite(qlog2) & is.na(qlog2)==F],breaks=100)
mean(qlog1[is.finite(qlog1) & is.na(qlog1)==F])
mean(qlog2[is.finite(qlog2) & is.na(qlog2)==F])

plot(log((germd+olsd)/(prevnsd+prevolsd))~log(prevnsd+prevolsd),data=msy)
plot(log((germd+olsd)/(prevnsd+mean(prevolsd,na.rm=T)))~log(prevnsd),data=msy)
plot(log((germd+olsd)/(mean(prevnsd,na.rm=T)+prevolsd))~log(prevolsd),data=msy)

with(subset(msy,year==2007),hist((germd+olsd)/prevolsd,breaks=100))
subset(msy,year==2007,select=c("species","germd","olsd","prevolsd","prevnsd"))

subset(ssyp,year==2006,select=c("species","olsd","nplot","totlive"))
subset(ssyp,year==2007,select=c("species","olsd","nplot","totlive"))

# sampling idea

jnames <- c("species","year","plot")

c_s <- merge(
	subset(csyp,select=c("species","year","plot","germd","prevnsd")),
	subset(ssyp,select=c("species","year","prevolsd")),
	by=c("species","year")
	)

with(c_s,hist(qlogis(germd/(prevnsd+prevolsd)),breaks=100))
	# very small, but makes sense because no seed mortality

### weird model thing yo

library(lme4)
csyp$obslevel <- as.factor(1:nrow(csyp))
ssyp$obslevel <- as.factor(1:nrow(ssyp))
sb$obslevel <- as.factor(1:nrow(sb))

m_g <- glmer(ngerm~offset(log(area*1000))+(1|obslevel),family=poisson,data=csyp)
m_o <- glmer(totlive~offset(log(area*1000))+(1|obslevel),family=poisson,data=ssyp)
m_o2 <- glmer(nlive~offset(log(corearea*1000))+(1|obslevel),family=poisson,data=sb)

m_oB <- glmer(nlive~offset(log(area*1000))
	+(1|species/plot/replicate)
	+(1|year)
	+(1|obslevel),
	family=poisson,	
	data=sb)

summary(m_g)
summary(m_o)
summary(m_o2)

stddevs <- function(obj) unlist(lapply(VarCorr(obj), function(m) sqrt(diag(m))))

nrep <- 10^4

g_mu <- fixef(m_g)
g_sig <- stddevs(m_g)

o_mu <- fixef(m_o)
o_sig <- stddevs(m_o)

g_mu_fixed <- exp(g_mu)*1000 # adjusting for area
g_var <- rpois(nrep,lambda=exp(rnorm(nrep,mean=g_mu,sd=g_sig))*1000)
g_mu_var <- mean(g_var)

o_mu_fixed <- exp(o_mu)*1000 # adjusting for area
o_var <- rpois(nrep,lambda=exp(rnorm(nrep,mean=o_mu,sd=o_sig))*1000)
o_mu_var <- mean(o_var)

g_mu_fixed 
g_mu_var
o_mu_fixed
o_mu_var

# difference increases with sample size

g_mu_fixed / (o_mu_fixed + g_mu_fixed)
per <- g_var / (o_var + g_var); mean(per[is.nan(per)==F])
g_mu_var / (o_mu_var + g_mu_var)

# var bigger for plot level but smaller for rep level

m_o <- glmer(totlive~(1|obslevel),family=poisson,data=ssyp)
o_mu <- fixef(m_o)
o_sig <- stddevs(m_o)
with(subset(sb,is.na(nlive)==F),lines(density(nlive),col="red"))
with(subset(sb,is.na(nlive)==F),hist(nlive,breaks=1000,prob=T,xlim=c(0,10)))
curve(dlnorm(x,meanlog=o_mu,sdlog=o_sig),add=T,col="red")

curve(dlnorm(x,meanlog=g_mu,sdlog=g_sig),col="darkgreen")
curve(dlnorm(x,meanlog=o_mu,sdlog=o_sig),add=T,col="black")

### raw survival rates

plot(log10(olsd/prevolsd)~log10(prevolsd),data=subset(msy,prevnsd==0))
plot(log10(olsd/prevolsd)~log10(prevolsd),data=subset(ssyp,prevnsdbar==0))

plot(log10(olsd/prevolsd)~log10(prevolsd),data=subset(ssyp,prevnsdbar==0))

### difference in germ and olsd dists

with(subset(csyp,is.na(germd)==F & germd<200),plot(density(germd)))
with(subset(ssyp,is.na(olsd)==F & olsd<200),lines(density(olsd),col="blue"))
sb$olsd <- with(sb,nlive/area)
with(subset(sb,is.na(olsd)==F),plot(density(olsd),col="red"))

####################
### SEED SPATIAL ###
####################

i <- 1
cursp <- levels(sb$species)[i]
curdat <- subset(sb,species==cursp)

plot(nlive~plot,border = "white",data=curdat)
points(jitter(nlive,amount=0)~plot,data=curdat)

#################################
### GERMINATION SAMPLING PROB ###
#################################

myyear <- 2005
with(subset(sb,species=="pere" & year==myyear),mean(nlive,na.rm=T))
with(subset(sb,species=="pere" & year==myyear),mean(nlive[-which(nlive==max(nlive,na.rm=T))],na.rm=T))

all <- merge(csyp,ssyp,by=c("species","year"))
all$G <- with(all, germd/(olsd+germd))
hist(all$G,breaks=100)
hist(qlogis(all$G),breaks=100)
with(all,tapply(G,species,mean,na.rm=T))
with(subset(all,is.nan(G)==F & is.na(G)==F & G>0 & G<1),
	plogis(tapply(qlogis(G),species,mean,na.rm=T))	
	)

all$G <- with(all, germd/(prevolsd+prevnsd))
hist(all$G,breaks=100)
hist(qlogis(all$G),breaks=100)
with(all,tapply(G,species,mean,na.rm=T))
with(subset(all,is.nan(G)==F & is.na(G)==F & G>0 & G<1),
	plogis(tapply(qlogis(G),species,mean,na.rm=T))	
	)

###################
### SIMULATIONS ###
###################

n_rep <- 10^2
n_sim <- 10^2
G_true <- 0.5

store1 <- store2 <- vector(length=n_rep)

for (i in 1:n_rep){
	n_g_true <- exp(rnorm(n_sim,mean=o_mu,sd=o_sig*1.5))*1000*G_true
	n_o_true <- exp(rnorm(n_sim,mean=o_mu,sd=o_sig))*1000*(1-G_true)
	n_g_obs <- rpois(n_sim,n_g_true)
	n_o_obs <- rpois(n_sim,n_o_true)

	store1[i] <- mean(n_g_obs)/(mean(n_o_obs)+mean(n_g_obs))
	hi <- n_g_obs/(n_o_obs+n_g_obs); store2[i] <- mean(hi[is.nan(hi)==F])

	}

plot(log(n_g_obs) ~ log(n_o_obs))

par(mfrow=c(1,2))
hist(store1,breaks=100,xlim=c(0,1))
abline(v=G_true,lty=2,col="red")
hist(store2,breaks=100,xlim=c(0,1))
abline(v=G_true,lty=2,col="red")

G_true - mean(store1,na.rm=T)
G_true - mean(store2,na.rm=T)

	# using exp(expected o value) is biased
	# but better estimates if randomly pair
	# however, this assumes that G and O follow same distribution

################
### TEMPORAL ###
################

### 	NEW VARIABLES

ggrow <- with(msy,germd/(prevcsd))
ggrow[is.finite(ggrow)==F] <- NA
ogrow <- with(msy,olsd/(prevcsd))
ogrow[is.finite(ogrow)==F] <- NA

# sgrow <- with(msy,germd/(germd+olsd))
# sgrow[is.finite(sgrow)==F] <- NA
# fgerm <- with(msy,(nsd+olsd)/(prevcsd))
# fgerm[is.finite(fgerm)==F] <- NA

### LISTS

pr_year_list <- split(csy$pr,csy$species)
rs_year_list <- split(csy$rs,csy$species)
pcr_year_list <- split(csy$pcr,csy$species)

ggrow_year_list <- split(ggrow,msy$species)
ogrow_year_list <- split(ogrow,msy$species)

# sgrow_year_list <- split(sgrow,msy$species)
# fgerm_year_list <- split(fgerm,msy$species)

### MATRICES 

cor_pr_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=pr_year_list)
cor_rs_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=rs_year_list)
cor_pcr_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=pcr_year_list)
cor_ggrow_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=ggrow_year_list)
cor_ogrow_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=ogrow_year_list)

# cor_sgrow_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=sgrow_year_list)
# cor_fgerm_year_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=fgerm_year_list)

par(mfrow=c(1,2))
hist(cor_ggrow_year_mat[20,],xlim=c(-1,1))
abline(v=0,lty=2,col="red")
hist(cor_ggrow_year_mat[16,],xlim=c(-1,1))
abline(v=0,lty=2,col="red")

###############
### SPATIAL ###
###############

### NEW VARIABLES AND DATAFRAMES

csyp$ggrow <- with(csyp,germd/(prevnsdbar+prevolsdbar))
csyp$ggrow[is.finite(csyp$ggrow)==F] <- NA
ssyp$ogrow <- with(ssyp,olsd/(prevnsdbar+prevolsdbar))
ssyp$ogrow[is.finite(ssyp$ogrow)==F] <- NA

csp <- ddply(csyp, .(species,plot), summarize,
	pr = mean(pr,na.rm=T),
	rseed = mean(rseed,na.rm=T),
	pcr = mean(pcr,na.rm=T),
	ggrow = mean(ggrow,na.rm=T)
	)

ssp <- ddply(ssyp, .(species,plot), summarize,
	ogrow = mean(ogrow,na.rm=T)
	)

### LISTS

pr_plot_list <- split(csp$pr,csp$species)
rs_plot_list <- split(csp$rs,csp$species)
pcr_plot_list <- split(csp$pcr,csp$species)

ggrow_plot_list <- split(csp$ggrow,csp$species)
ogrow_plot_list <- split(ssp$ogrow,ssp$species)

### MATRICES 

cor_pr_plot_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=pr_plot_list)
cor_rs_plot_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=rs_plot_list)
cor_pcr_plot_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=pcr_plot_list)
cor_ggrow_plot_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=ggrow_plot_list)
cor_ogrow_plot_mat <- outer(1:nspecies,1:nspecies,paircor,veclist=ogrow_plot_list)

############
### PLOT ###
############

pdf(paste0("Plots/correlations_year_plot_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=12,height=5)

par(mfrow=c(2,5),mar=c(2,2,2,2),oma=c(6,6,2,2),las=1,bty="l")

corhist(cor_pr_year_mat,xlab="pr")
addylab("Probability density")
lettlab2(1)
corhist(cor_rs_year_mat,xlab="rs")
lettlab2(2)
corhist(cor_pcr_year_mat,xlab="pcr")
lettlab2(3)
corhist(cor_ggrow_year_mat,xlab="ggrow")
lettlab2(4)
corhist(cor_ogrow_year_mat,xlab="ogrow")
lettlab2(5)

corhist(cor_pr_plot_mat,xlab="pr")
addylab("Probability density")
flexxlab(expression(cor~bgroup("(","Pr(Y>1)",")")))
lettlab2(6)
corhist(cor_rs_plot_mat,xlab="rs")
flexxlab(expression(cor~bgroup("(","Y|Y>1",")")))
lettlab2(7)
corhist(cor_pcr_plot_mat,xlab="pcr")
flexxlab(expression(cor~bgroup("(","Y",")")))
lettlab2(8)
corhist(cor_ggrow_plot_mat,xlab="ggrow")
flexxlab(expression(cor~bgroup("(",frac(N[g](t),N[o](t-1)+N[n](t-1)),")")),myline=6)
lettlab2(9)
corhist(cor_ogrow_plot_mat,xlab="ogrow")
flexxlab(expression(cor~bgroup("(",frac(N[o](t),N[o](t-1)+N[n](t-1)),")")),myline=6)
lettlab2(10)

dev.off()

	# Suggests that spatial should be separate for each species

	# A good year for one species is a good year for another,
	# but a good site for one can be a bad site for another

### CLEAN

remove(csp,ssp) # and all the lists
gc()

