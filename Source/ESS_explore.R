### Develop explanations for patterns in ESS germination results ###

source("Source/invasion_functions_iterative.R")

# Empirical simulations ---------------------------------------------------

outlist <- evolve(
  nr=1000,nt=10100,nb=100,
  zam=zamo+mam*zsdo,wam=wamo+mam*wsdo,
  zsd=zsdo*msd,wsd=wsdo*msd,rho=0.82,
  beta_p=pl$pr$beta_p[2,19,],beta_r=pl$rs$beta_r[2,19,],
  sig_y_p=pl$pr$sig_y_p[2,19],sig_y_r=pl$rs$sig_y_r[2,19],
  sig_o_p=pl$pr$sig_o_p[2],phi_r=pl$rs$phi_r[1],
  m0=exp(pl$go$alpha_m[2,19]),m1=exp(pl$go$beta_m[2,19]),
  am0=1,bm0=1,
  as0=1,bs0=1,
  abr0=0.1,
  smut_m=1,smut_s=0.1,smut_r=0.1,
  savefile=NULL,
  nsmin=10^-50
)

zam=zamo+mam*zsdo
wam=wamo+mam*wsdo
zsd=zsdo*msd
wsd=wsdo*msd
rho=0.82
zw_mu <- c(zam,wam) - log(tau_p)
zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
zw <- mvrnorm(n=nt, mu=zw_mu, Sigma=zw_sig)

eps_y_p <- rnorm(nt,0,sig_y_p)
eps_y_r <- rnorm(nt,0,sig_y_r)

rd <- ressim(w=zw[,2],x_z,am=1,bm=1,as=1,bs=1,abr=1,
  beta_p=pl$pr$beta_p[1,19,],beta_r=pl$rs$beta_r[1,19,],
  eps_y_p=rep(0,nt),eps_y_r=rep(0,nt),
  sig_o_p=pl$pr$sig_o_p[1],phi_r=pl$rs$phi_r[1],
  So=0.1,m0=exp(pl$go$alpha_m[1,19]),m1=exp(pl$go$beta_m[1,19]),
  nsmin=10^-50
)

plot(outlist$es$am,type="l")
plot(outlist$es$bm,type="l")
plot(outlist$es$as,type="l")
plot(outlist$es$bs,type="l")
plot(outlist$es$abr,type="l")

plot(outlist2$es$am)
plot(outlist2$es$bm)
plot(outlist2$es$as)

with(outlist$es[nr,],
  curve(coaG(w=x,am=am,bm=bm,as=as,bs=bs,abr=abr),xlim=c(-5,5),ylim=c(0,1))
)
with(outlist$es[nr,],curve(fixG(x,am,bm),xlim=c(-5,5),col="red",ylim=c(0,1)))
with(outlist2$es[nr,],curve(fixG(x,am,bm),add=T,col="red"))
with(outlist3$es[nr,],curve(fixG(x,am,bm),add=T,col="red"))
with(outlist4$es[nr,],curve(fixG(x,am,bm),add=T,col="red"))
with(outlist5$es[nr,],curve(fixG(x,am,bm),add=T,col="red"))

with(outlistB$es[nr,],curve(fixG(x,am,bm),add=T,col="blue"))
with(outlistB2$es[nr,],curve(fixG(x,am,bm),add=T,col="blue"))
with(outlistB3$es[nr,],curve(fixG(x,am,bm),add=T,col="blue"))
with(outlistB4$es[nr,],curve(fixG(x,am,bm),add=T,col="blue"))
with(outlistB5$es[nr,],curve(fixG(x,am,bm),add=T,col="blue"))

with(outlistC$es[nr,],curve(fixG(x,am,bm),add=T,col="green"))
with(outlistC2$es[nr,],curve(fixG(x,am,bm),add=T,col="green"))
with(outlistC3$es[nr,],curve(fixG(x,am,bm),add=T,col="green"))
with(outlistC4$es[nr,],curve(fixG(x,am,bm),add=T,col="green"))
with(outlistC5$es[nr,],curve(fixG(x,am,bm),add=T,col="green"))

abline(v=quantile(outlist$zw[,1],probs=c(0.05,0.95)),lty=3)

rd1 <- with(es[1,], ressim(zw[,2],x_z,am,bm,as,bs,abr,
  beta_p,beta_r,
  eps_y_p,eps_y_r,
  sig_o_p,phi_r,
  So,m0,m1,
  nt,nsmin,
  full=T
) )
rd2 <- with(es[nr,], ressim(zw[,2],x_z,am,bm,as,bs,abr,
  beta_p,beta_r,
  eps_y_p,eps_y_r,
  sig_o_p,phi_r,
  So,m0,m1,
  nt,nsmin,
  full=T
) )

with(rd2,plot(log(Ye)~zw[,2]))
with(rd2,lines(supsmu(zw[,2],log(Ye)),col="red"))
abline(h=0,col="blue",lty=3)
abline(v=0,col="blue",lty=3)

G1 <- rd1$Gres
G2 <- rd2$Gres

d1 <- log((1-G2)*So + G2*rd1$Ye) - log((1-G1)*So + G1*rd1$Ye)
d2 <- log((1-G2)*So + G2*rd2$Ye) - log((1-G1)*So + G1*rd2$Ye)

par(mfrow=c(1,1))
plot(d1 ~ w)
lines(supsmu(w,d1),col="orange")
lines(supsmu(w,d2),col="red")
abline(h=0,col="blue",lty=3)
curve(log(fixG(x,es[nr,]$am,es[nr,]$bm)/fixG(x,es[1,]$am,es[1,]$bm)),add=T,col="purple")
curve(log((1-fixG(x,es[nr,]$am,es[nr,]$bm))/(1-fixG(x,es[1,]$am,es[1,]$bm))),add=T,col="purple")
curve(log(
  (So*(1-fixG(x,es[nr,]$am,es[nr,]$bm))+fixG(x,es[nr,]$am,es[nr,]$bm))
  /(So*(1-fixG(x,es[1,]$am,es[1,]$bm))+fixG(x,es[1,]$am,es[1,]$bm))
),add=T,col="purple")
# curve for Y=1

rd1med <- median(rd1$ns)
c1 <- rd1$ns<quantile(rd1$ns,0.05) # rd1$ns<=rd1med
c2 <- rd1$ns>quantile(rd1$ns,0.95) # rd1$ns>rd1med
d1a <- d1[c1]
d1b <- d1[c2]
w1 <- w[c1]
w2 <- w[c2]
plot(d1a ~ w1,col="blue")
points(d1b ~ w2,col="purple")
lines(supsmu(w1,d1a),col="black")
lines(supsmu(w2,d1b),col="black")
abline(h=0,col="blue",lty=3)

par(mfrow=c(2,2),mar=c(2,2,2,2))
with(es[1,],curve(fixG(x,am,bm),xlim=c(-5,5),ylim=c(0,1),col="orange"))
with(es[nr,],curve(fixG(x,am,bm),add=T,col="red"))
abline(v=quantile(zw[,2],probs=c(0.05,0.95)),lty=3)
with(rd2,plot(log(Ye)~zw[,2],
  type="n",
  xlim=quantile(zw[,2],probs=c(0.005,0.995)),
  ylim=c(-2,2)))
with(rd1, lines(supsmu(zw[,2],log(Ye)),col="orange"))
with(rd2, lines(supsmu(zw[,2],log(Ye)),col="red"))
abline(h=0,col="blue",lty=3)
abline(v=0,col="blue",lty=3)
plot(density(log(rd1$ns)),xlim=c(1.5,5.5),col="orange")
lines(density(log(rd2$ns)),col="red")
plot(density(log(rd1$ns*rd1$G*rd1$Ye),n=2048),xlim=c(4,5),col="orange")
lines(density(log(rd2$ns*rd2$G*rd2$Ye),n=2048),col="red")

plot(density(log(rd1$ns*rd1$G*rd1$Ye),n=2048),xlim=c(-4,4),col="orange",ylim=c(0,0.1))
lines(density(log(rd2$ns*rd2$G*rd2$Ye),n=2048),col="red")

K <- exp(-m0*T3) / ( (m1/m0)*(1-exp(-m0*T3))/(tau_s/10) ) # = K_Y
# dN/dt = K - (So+G)N
G <- 1
K/(So+G) # overall K if env always very favourable
quantile(rd1$ns,0.95) # even at this value, Y>=1

par(mfrow=c(2,1))
ws <- w < median(w)
wl <- w >= median(w)
with(rd1[ws,], hist(log(Ye/So),breaks=1000,col=rgb(0,0,1,alpha=0.25),border=NA,xlim=c(-5,5)))
with(rd1[wl,], hist(log(Ye/So),breaks=1000,add=T,col=rgb(1,0,0,alpha=0.25),border=NA))
abline(v=0)
ws <- w < median(w)
wl <- w >= median(w)
with(rd2[ws,], hist(log(Ye/So),breaks=1000,col=rgb(0,0,1,alpha=0.25),border=NA,xlim=c(-5,5)))
with(rd2[wl,], hist(log(Ye/So),breaks=1000,add=T,col=rgb(1,0,0,alpha=0.25),border=NA))
abline(v=0)

with(rd1[ws,], mean(log(Ye/So)))
with(rd1[wl,], mean(log(Ye/So)))
with(rd2[ws,], mean(log(Ye/So)))
with(rd2[wl,], mean(log(Ye/So)))

plot(log(Ye/So)~w,data=rd1,type="n",ylim=c(-5,1))
with(rd1, lines(supsmu(w,log(Ye/So)),col="orange"))
with(rd2, lines(supsmu(w,log(Ye/So)),col="red"))
abline(v=mean(w),lty=3)

h1 <- log((1-G1)*So + G1*rd2$Ye) - log((1-G1)*So + G1*rd1$Ye)
h2 <- log((1-G2)*So + G2*rd2$Ye) - log((1-G2)*So + G2*rd1$Ye)
plot(h1~w,type="n")
lines(supsmu(w,h1),col="orange")
lines(supsmu(w,h2),col="red")
mean(h2-h1)

for(t in 1:nt){
  x_t <- c(x_z[t,],-10)
  pi_bar_t <- sum(beta_p * x_t) + eps_y_p[t]
  mu_bar_t <- sum(beta_r * x_t) + eps_y_r[t]
  pr_t <- logitnormint(mu=pi_bar_t,sigma=sig_o_p)
  rs_t <- nbtmean(exp(mu_bar_t),phi_r)
  Ye2[t] <- pr_t * rs_t * BHS(pr_t * rs_t,m0,m1) 
}
par(new=F)
plot(density(log(Ye)))
par(new=T)
plot(density(log(Ye2)),col="blue")

# nr=1000;nt=250;nb=50;
# zam=zamo+mam*zsdo;wam=wamo+mam*wsdo;
# zsd=zsdo*msd;wsd=wsdo*msd;rho=0.82;
# beta_p=pl$pr$beta_p[1,19,];beta_r=pl$rs$beta_r[1,19,];
# sig_y_p=pl$pr$sig_y_p[1,19];sig_y_r=pl$rs$sig_y_r[1,19];
# sig_o_p=pl$pr$sig_o_p[1];phi_r=pl$rs$phi_r[1];
# m0=exp(pl$go$alpha_m[1,19]);m1=exp(pl$go$beta_m[1,19]);
# am0=1;bm0=1;
# as0=1;bs0=1;
# abr0=0.1;
# smut_m=1;smut_s=0.1;smut_r=0.1;
# savefile=NULL;
# nsmin=10^-50
# 
# w=zw[,2]
# attach(es[1,])

# Simple simulations ------------------------------------------------------

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

Gres <- Ginv <- plogis(-5,5,length.out=100)
gd <- expand.grid(Gres,Ginv)

i <- 1
j <- 19
So <- exp(-exp(pl$go$m1[i,j]))


# Ellner 1985 exploration -------------------------------------------------

msy <- read.csv("Output/msy_seedests_18Jan2017.csv",header=T)

msy$Y <- with(msy, nsdbar/germdbar)
msy$S <- 0.5
msy$lambda <- with(msy,csdbar/prevcsdbar) # this probably wrong
msy$X <- msy$csd

msy$wc <- with(msy,ifelse(gprcp>median(gprcp,na.rm=T),0,1))

Glo <- with(subset(msy,
                   wc==0 & !is.na(Y) & !is.na(lambda) & !is.na(X) & lambda>0),
     mean(S/lambda) / mean(Y/lambda)
)
Ghi <- with(subset(msy,
                   wc==1 & !is.na(Y) & !is.na(lambda) & !is.na(X) & lambda>0),
            mean(S/lambda) / mean(Y/lambda)
)

with(subset(msy,wc==0),hist(log(Y/lambda),breaks=100))
with(subset(msy,wc==1),hist(log(Y/lambda),breaks=100))


sims <- readRDS("Sims/res_mu1_sd1_s7_16Aug2017.rds")
sims$wc <- with(sims,ifelse(w>median(w,na.rm=T),0,1))
sims$Ye <- with(sims, Y*Sn)
sims$lambda <- with(sims, G*Ye + (1-G)*So)
relY <- with(sims, Ye/lambda)
relS <- with(sims, So/lambda)
i <- 1
j <- 19

with(sims,mean(relS[i,30:250,j]) / mean(relY[i,30:250,j]))


1/mean(1/exp(rnorm(10^3,0,1)))
1/mean(1/exp(rnorm(10^3,0,3)))
hist(1/exp(rnorm(10^3,0,1)),breaks=1000)


