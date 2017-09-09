### Calculate ES alpha_G and beta_G by iteratively perturbing parameters  ###
### and attempting re-invasion                                            ###

BHS <- function(n,m0,m1){
  exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*n/(tau_s/10) )
}
# original model in m^2
# DD functions use density in 0.01 m^2 = 10 x 10 cm plots
# But we want to use 0.1m^2 plots (to match scale of quadrats)
# Therefore, tau_s set to 10 instead of 100
  
logitnorm <- function(x,mu,sigma){
  plogis(x) * dnorm(x, mean=mu, sd=sigma)
}
# code borrowed from logitnorm package
  
logitnormint <- Vectorize(function(mu,sigma,intsd=10,...){
  integrate(logitnorm,
            mu=mu,sigma=sigma,
            lower=mu-intsd*sigma,
            upper=mu+intsd*sigma,
            ...)$value
})
  
nbtmean <- function(mu,phi){
  mu / ( 1 - (phi/(mu+phi))^phi )
}
# mean for hurdle model

fixG <- function(w,a,b){
  plogis(a+b*w)
}

coaG_w <- function(w,am,bm,as,bs,abr,nint=10,intsd=2){
  require(mvtnorm)
  xa <- seq(am-intsd*as,am+intsd*as,length.out=nint)
  xb <- seq(bm-intsd*bs,bm+intsd*bs,length.out=nint)
  ab <- expand.grid(xa=xa,xb=xb)
  abm <- c(am,bm)
  abs <- matrix(c(as^2,rep(abr*as*bs,2),bs^2),nr=2,nc=2)
  pg <- dmvnorm(ab, mean=abm, sigma=abs, log=F) 
  sum(pg * fixG(w,ab$xa,ab$xb))/sum(pg)
}
  
coaG <- Vectorize(coaG_w,vectorize.args="w")
# coaG(w=zw[,2],am,bm,as,bs,abr,nint=10,intsd=2)

ressim <- function(w,x_z,am,bm,as,bs,abr,
                   beta_p,beta_r,
                   eps_y_p,eps_y_r,
                   sig_o_p,phi_r,
                   So,m0,m1,
                   nt,nsmin
                   ){
  ns <- ng <- nn <- Ye <- rep(NA,nt)
  #Gres <- coaG(w,am,bm,as,bs,abr)
  Gres <- fixG(w,am,bm)
  ns[1] <- nstart
  t <- 1
  while(ns[t] > nsmin & t <= nt){
    ng[t] <- Gres[t] * ns[t]
    x_t <- c(x_z[t,],log(ng[t])-log(tau_d/10))
    pi_bar_t <- sum(beta_p * x_t) + eps_y_p[t]
    mu_bar_t <- sum(beta_r * x_t) + eps_y_r[t]
    pr_t <- logitnormint(mu=pi_bar_t,sigma=sig_o_p)
    rs_t <- nbtmean(exp(mu_bar_t),phi_r)
    nn[t] <- ng[t] * pr_t * rs_t
    Ye[t] <- nn[t] * BHS(nn[t],m0,m1) / ng[t]
    if(t<nt) ns[t+1] <- ns[t] * ( (1-Gres[t])*So + Gres[t]*Ye[t] )
    t <- t + 1
  }
  return(data.frame(Gres=Gres,Ye=Ye))
}

invade <- function(w,ami,bmi,asi,bsi,abri,Gres,Ye,So,nt,nb){
  if(!NA %in% Ye){ # t = final value at which loop stopped 
    #Ginv <- coaG(w,ami,bmi,asi,bsi,abri)
    Ginv <- fixG(w,ami,bmi)
    delta_r <- log((1-Ginv)*So + Ginv*Ye) - log((1-Gres)*So + Gres*Ye)
    invaded <- mean(delta_r[(nb+1):nt]) > 0
  }
  if(NA %in% Ye){   
    invaded <- TRUE 
    # if resident goes extinct, invader establishes immediately
  }
  return(invaded)
}

evolve <- function(
  nr,nt,nb,
  zam,wam,zsd,wsd,rho=0.82,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_o_p,phi_r,
  m0,m1,
  am0,bm0,
  as0,bs0,
  abr0,
  smut_m=0.5,smut_s=0.1,smut_r=0.1,
  savefile=NULL,
  nsmin=10^-50
  ){
  
  require(MASS)
  cur_date <- format(Sys.Date(),"%d%b%Y")
  
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
  
  es <- data.frame(am=rep(NA,times=nr),
                   bm=rep(NA,times=nr),
                   as=rep(NA,times=nr),
                   bs=rep(NA,times=nr),
                   abr=rep(NA,times=nr)
                   )
  es[1,] <- c(am0,bm0,as0,bs0,abr0)

  for(i in 1:nr){
    
    if(i==1){
      rd <- with(es[i,], ressim(zw[,2],x_z,am,bm,as,bs,abr,
                                beta_p,beta_r,
                                eps_y_p,eps_y_r,
                                sig_o_p,phi_r,
                                So,m0,m1,
                                nt,nsmin
                                ) )
        # simulate starting resident dynamics
    }
    ami <- es$am[i] + rnorm(1,0,smut_m)
    bmi <- es$bm[i] + rnorm(1,0,smut_m)
    asi <- es$as[i] * exp(rnorm(1,0,smut_s))
    bsi <- es$bs[i] * exp(rnorm(1,0,smut_s))
    abri <- plogis( (qlogis((es$abr[i]+1)/2) + rnorm(1,0,smut_r)) )*2 - 1 
      # transform [0,1] to correlation range of [-1,1]
    if(abri==-1) abri <- -0.99
    if(abri==+1) abri <- +0.99
    invaded <- invade(zw[,2],ami,bmi,asi,bsi,abri,rd$Gres,rd$Ye,So,nt,nb)
    if(i < nr){
      if(invaded==TRUE){
        es[i+1,] <- c(ami,bmi,asi,bsi,abri)
        rd <- ressim(zw[,2],x_z,ami,bmi,asi,bsi,abri,
                     beta_p,beta_r,
                     eps_y_p,eps_y_r,
                     sig_o_p,phi_r,
                     So,m0,m1,
                     nt,nsmin) # simulate new resident dynamics
      } 
      if(invaded==FALSE){
        es[i+1,] <- es[i,]
      } 
    }
  } # close i loop
  
  outlist <- list(zw=zw,es=es)  
  
  if(is.null(savefile)){
    return(outlist)
  }
  
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
}

outlist <- evolve(
  nr=1000,nt=10100,nb=100,
  zam=zamo+mam*zsdo,wam=wamo+mam*wsdo,
  zsd=zsdo*msd,wsd=wsdo*msd,rho=0.82,
  beta_p=pl$pr$beta_p[1,19,],beta_r=pl$rs$beta_r[1,19,],
  sig_y_p=pl$pr$sig_y_p[1,19],sig_y_r=pl$rs$sig_y_r[1,19],
  sig_o_p=pl$pr$sig_o_p[1],phi_r=pl$rs$phi_r[1],
  m0=exp(pl$go$alpha_m[1,19]),m1=exp(pl$go$beta_m[1,19]),
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

with(rd,plot(log(Ye)~zw[,2]))
with(rd, lines(supsmu(zw[,2],log(Ye)),col="red"))
abline(h=0,col="blue",lty=3)
abline(v=0,col="blue",lty=3)

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



