nt <- nb + 100
nk <- 100
n <- 1
mam <- maml[[mpos[n]]]
msd <- msdl[[mpos[n]]]
iset <- itersetl[[cpos[n]]]

par1 <- with(pls,list(
ni=ni,nj=nj,nr=nr,nt=nt,nb=nb,nk=nk,
zam=zamo+mam*zsdo,zsd=zsdo*msd,
wam=wamo+mam*wsdo,wsd=wsdo*msd,
beta_p=beta_p[iset,,],beta_r=beta_r[iset,,],
sig_y_p=sig_y_p[iset,],sig_y_r=sig_y_r[iset,],
sig_s_g=sig_s_g[iset,],sig_s_p=sig_s_p[iset],sig_s_r=sig_s_r[iset],
sig_o_p=sig_o_p[iset],phi_r=phi_r[iset],theta_g=theta_g[iset,],
m0=m0[iset,],m1=m1[iset,],
am0=am0[iset],bm0=bm0[iset],
DDFUN=BHS,
Sg=1
))

i <- 1; j <- 19;

par2 <- with(par1,list(
nr=nr,nt=nt,nb=nb,nk=nk,
zam=zam,wam=wam,zsd=zsd,wsd=wsd,rho=0.82,
beta_p=beta_p[i,j,],beta_r=beta_r[i,j,],
sig_y_p=sig_y_p[i,j],sig_y_r=sig_y_r[i,j],
sig_s_g=sig_s_g[i,j],sig_s_p=sig_s_p[i],sig_s_r=sig_s_r[i],
sig_o_p=sig_o_p[i],phi_r=phi_r[i],
theta_g=ifelse(is.null(theta_g),NULL,theta_g[i,j]),
m0=m0[i,j],m1=m1[i,j],
am0=am0[i],bm0=bm0[i],
DDFUN=BHS,
Sg=1
))

attach(par1)
attach(par2)

require(MASS)

zw_mu <- c(zam,wam) - log(tau_p)
zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
zw <- mvrnorm(n=nt, mu=zw_mu, Sigma=zw_sig)

eps_y_p <- rnorm(nt,0,sig_y_p)
eps_y_r <- rnorm(nt,0,sig_y_r)
So <- exp(-m0)

x_z <- matrix(nr=nt,nc=3)
x_z[,1] <- 1 # intercept
x_z[,2] <- zw[,1]
x_z[,3] <- zw[,1]^2 

es <- data.frame(am=rep(NA,times=nr),bm=rep(NA,times=nr))
# as=rep(NA,times=nr), bs=rep(NA,times=nr),abr=rep(NA,times=nr)

es[1,] <- c(am0=0,bm0=1) # ,as0,bs0,abr0)

am <- am0
bm <- bm0
ami <- es$am[i] + rnorm(1,0,smut_m)
bmi <- es$bm[i] + rnorm(1,0,smut_m)
w=zw[,2]

nc=5   # n consecutive t that ns must be < nsmin
nstart=1
intsd=10
tau_d=100

set.seed(1)
library("profvis")
profvis({
  invade_finite(w,x_z,am,bm,ami,bmi,
                beta_p,beta_r,
                eps_y_p,eps_y_r,
                sig_s_g,sig_s_p,sig_s_r,
                sig_o_p,phi_r,theta_g,
                So,m0,m1,
                nt,nb,nk,nsmin,ngmin,
                DDFUN,
                Sg
  )
  })
?profvis

### Testing truncated negbin distribution

pradj <- function(pr,mu,phi){
  q <- dnbinom(0, mu=mu, size=phi) # Pr(Y>0)
  return(pr / (1-q)) # zero-inflated
}

pradj2 <- function(pr,mu,phi){
  q <- dnbinom(0, mu=mu, size=phi) # Pr(Y>0)
  return((1-pr) / q) # proportion of required zeroes provided by pr
}

pr_obs <- 3/4
pr_lef <- pradj(pr_obs,mu,phi)
pr_fra <- pradj2(pr_obs,mu,phi)

nsim <- 10^4
nrep <- 100
mu <- 2
phi <- 1
prob <- phi/(phi+mu)
y1 <- y2 <- y3 <- rep(NA,nsim)
require(countreg)
for(i in 1:nsim){
  y1[i] <- sum(rbinom(nrep, size=1, prob=pr_obs) * rtrunc(n=nrep,spec="nbinom",prob=prob,size=phi,a=0))
  # y2[i] <- rnbinom(1, prob=prob, size=phi*rbinom(1, size=nrep, prob=pr_lef))
  nneg <- rbinom(1, size=nrep, prob=pr_fra)
  y3[i] <- rnbinom(1, prob=prob, size=phi*nneg) 
  + sum(rtrunc(n=nrep-nneg,spec="nbinom",prob=prob,size=phi,a=0))
}

hist(y1,breaks=1000)
# hist(y2,breaks=1000,add=T,border="red")
hist(y3,breaks=1000,add=T,border="blue")

mean(y1)
mean(y2)
sd(y1)
sd(y2)
# Testing re-parameterisation to standard negative binomial distribution

x1 <- x2 <- rep(NA,nsim)
require(countreg)
for(i in 1:nsim){
  x1[i] <- sum(rnbinom(nrep, mu=mu, size=phi))
  x2[i] <- sum(rbinom(nrep, size=1, prob=2/3) * rztnbinom(n=nrep, mu=mu, size=phi))
}

hist(x1,breaks=1000)
hist(x2,breaks=1000,add=T,border="red")

mean(x1)
mean(x2)
sd(x1)
sd(x2)

dnbinom(0, mu=mu, size=phi)
dnbinom(0, mu=nrep*mu, size=phi)
1/3^nsim

t1 <- rnbinom(10^6, mu=mu, size=phi) + rnbinom(10^6, mu=mu, size=phi)
t2 <- rbinom(10^6, size=1, prob=2/3) * rztnbinom(n=nrep, mu=mu, size=phi) + rbinom(10^6, size=1, prob=2/3) * rztnbinom(n=nrep, mu=mu, size=phi)

hist(t1,breaks=1000)
hist(t2,breaks=1000,add=T,border="red")

mean(t1)
mean(t2)
sd(t1)
sd(t2)

# sum negbin

a1 <- a2 <- rep(NA,nsim)
require(countreg)
for(i in 1:nsim){
  a1[i] <- sum(rnbinom(nrep, prob=prob, size=phi))
  a2[i] <- rnbinom(1, prob=prob, size=phi*nrep)
}

hist(a1,breaks=1000)
hist(a2,breaks=1000,add=T,border="red")

# negbin lognormal

sum(negbin(lognormal(mu,sig),phi))
sum( negbin(exp(mu+sig),phi) + negbin(exp(mu,-sig),phi) )
