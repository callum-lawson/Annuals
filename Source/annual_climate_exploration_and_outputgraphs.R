##############################################################
### Expolore distribution and time series of rainfall data ###
##############################################################

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Dropbox/NIOO/"
# maindir <- "D:/Dropbox/NIOO/"

setwd(paste0(maindir,"Analyses/Venable"))
source("venable_figure_functions_25Nov2015.R")

library(bbmle)
library(reshape2)

####################
### CLIMATE DATA ###
####################

ncy <- read.csv("ncy_15Jan2016.csv",header=T)
ncy <- subset(ncy,is.na(seasprcp)==F)
	# removes first value (missing because no previous winter)

### QUICK EXPLORE

plot(seasprcp~year,data=ncy,type="b")
points(germprcp~year,data=ncy,col="red")

hist(ncy$seasprcp,breaks=30)
hist(ncy$germprcp,breaks=30)
hist(log(ncy$seasprcp),breaks=30)
hist(log(ncy$germprcp),breaks=30)

ncy$yearcat <- cut(ncy$year,5)
par(mfrow=c(2,1),mar=c(1,1,1,1))
boxplot(seasprcp~yearcat,data=ncy,las=3)
boxplot(log(seasprcp)~yearcat,data=ncy,las=3)

#############################
### GROWING SEASON MODELS ###
#############################

m <- mean(ncy$seasprcp)
v <- var(ncy$seasprcp)
shape <- m^2/v
scale <- v/m

lm <- mean(log(ncy$seasprcp))
lsd <- sd(log(ncy$seasprcp))

curve(dgamma(x,shape=shape,scale=scale),xlim=c(0,50))
curve(dlnorm(x,meanlog=lm,sdlog=lsd),xlim=c(0,50))

mg <- mle2(
	seasprcp ~ dgamma(shape=a,scale=b),
	start=list(a=m^2/v,b=v/m),
	data=ncy
	)

ml <- mle2(
	seasprcp ~ dlnorm(meanlog=mu,sdlog=sigma),
	start=list(mu=lm,sigma=lsd),
	data=ncy
	)

AIC(mg, ml)
	# indistinguishable

### AUTOCORRELATION

par(mar=c(4,4,2,2))
ar(ncy$seasprcp)
	# no autocorrelation?
acf(ncy$seasprcp)
	# blue = significance bounds

fake <- rnorm(nrow(ncy),0,1)
ar(fake)
acf(fake)

#################################
### GERMINATION SEASON MODELS ###
#################################

m <- mean(ncy$germprcp)
v <- var(ncy$germprcp)
shape <- m^2/v
scale <- v/m

lm <- mean(log(ncy$germprcp))
lsd <- sd(log(ncy$germprcp))

curve(dgamma(x,shape=shape,scale=scale),xlim=c(0,50))
curve(dlnorm(x,meanlog=lm,sdlog=lsd),xlim=c(0,50))

mg <- mle2(
	germprcp ~ dgamma(shape=a,scale=b),
	start=list(a=m^2/v,b=v/m),
	data=ncy
	)

ml <- mle2(
	germprcp ~ dlnorm(meanlog=mu,sdlog=sigma),
	start=list(mu=lm,sigma=lsd),
	data=ncy
	)

AIC(mg, ml)
	# indistinguishable

### AUTOCORRELATION

par(mar=c(4,4,2,2))
ar(ncy$germprcp)
	# no autocorrelation?
acf(ncy$germprcp)
	# blue = significance bounds

fake <- rnorm(nrow(ncy),0,1)
ar(fake)
acf(fake)
# acf(ncy$germprcp)

#####################
### OUTPUT GRAPHS ###
#####################

pdf(paste0("precipitation_summaries_",format(Sys.Date(),"%d%b%Y"),".pdf"),
		width=7,height=7)

layout(matrix(c(rep(1,2),2,3),byrow=T,nc=2,nr=2))
par(mar=c(5,5,2,2),oma=c(0,0,1,1),las=1,bty="l")

### TIME SERIES
plot(seasprcp~year,data=ncy,type="b",pch=16,ylim=c(0,400),ylab="Winter precipitation (mm)")
yrange <- c(1983,2014)
srange <- c(1990,2014)
yheight <- 350
sheight <- 410
arrows(x0=yrange[1],x1=yrange[2],y0=yheight,y1=yheight,code=3,angle=90,length=0.025)
text(x=median(yrange),y=yheight,labels="seedling censuses",pos=3,xpd=T)
arrows(x0=srange[1],x1=srange[2],y0=sheight,y1=sheight,code=3,angle=90,length=0.025)
text(x=median(srange),y=sheight,labels="seed censuses",pos=3,xpd=T)
mtext(bquote(bold((.(letters[1])))),side=3,line=0.5,xpd=T,adj=0.05)

### AUTOCORRELATION
acf(ncy$seasprcp,
	main="",
	xlab="Lag (years)",
	ylab="Autocorrelation"
	)
mtext(bquote(bold((.(letters[2])))),side=3,line=0.5,xpd=T,adj=0.05)

### CORRELATION WITH GERM PRCP
plot(log(ncy$seasprcp)~log(ncy$germprcp),
	xlab="ln(Germination precipitation)",
	ylab="ln(Winter precipitation)",
	pch=16
	)
corval <- with(ncy, round(cor.test(log(seasprcp),log(germprcp))$estimate,2))
legend("bottomright",
	legend=substitute(paste(rho,"=",rhohat),list(rhohat=corval)),
	bty="n"
	)
mtext(bquote(bold((.(letters[3])))),side=3,line=0,xpd=T,adj=0.05)

dev.off()

###################
### SIMULATIONS ###
###################

lncv <- function(m,s,n=10^6){
	x <- rlnorm(n,meanlog=m,sdlog=s)
	sd(x)/mean(x)	
	}

mseq <- -5:5
sseq <- seq(0,1,length.out=11)

d <- expand.grid(mseq,sseq)
d$cv <- with(d, mapply(lncv,Var1,Var2))
a <- acast(d,Var1 ~ Var2)
matplot(a,type="l",col= colorRampPalette(c("blue","red"))(11))

### MEAN AND CV

sdcalc <- function(CV) sqrt(log(CV^2+1))
meancalc <- function(mu,CV) exp(mu + log(CV^2+1)/2)

sdcalc(1)
meancalc(mu=1,CV=2)




