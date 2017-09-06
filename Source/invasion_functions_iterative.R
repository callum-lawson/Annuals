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
  
logitnormint <- Vectorize(function(mu,sigma){
  integrate(logitnorm,
            mu=mu,sigma=sigma,
            lower=mu-intsd*sigma,
            upper=mu+intsd*sigma,
            rel.tol=rel.tol,
            abs.tol=abs.tol)$value
})
  
nbtmean <- function(mu,phi){
  mu / ( 1 - (phi/(mu+phi))^phi )
}
# mean for hurdle model
  
ressim <- function(alpha_G,beta_G){
  ns <- ng <- nn <- Ye <- rep(NA,nt)
  G <- plogis(alpha_G + beta_G*zw[,2])
  ns[1] <- nstart
  t <- 1
  while(ns[t] > nsmin & t <= nt){
    ng[t] <- G[t] * ns[t]
    x_t <- c(x_z[t,],log(ng[t])-log(tau_d/10))
    pi_bar_t <- sum(beta_p * x_t)
    mu_bar_t <- sum(beta_r * x_t)
    pr_t <- logitnormint(mu=pi_bar_t+eps_y_p[t],sigma=sig_o_p)
    rs_t <- nbtmean(exp(mu_bar_t),phi_r)
    nn[t] <- ng[t] * pr_t * rs_t
    Ye[t] <- nn[t] * BHS(nn[t],m0,m1) / ng[t]
    ns[t+1] <- ns[t] * ( (1-G[t])*So + G[t]*Ye[t] )
    t <- t + 1
  }
  return(data.frame(G=G,Ye=Ye))
}
  
invade <- function(alpha_G_inv,beta_G_inv,rd){
  if(!NA %in% rd$Ye){ # t = final value at which loop stopped  
    Ginv <- plogis(alpha_G_inv,beta_G_inv*zw[,2])
    delta_r <- with(rd, 
      log( (1-Ginv)*So + Ginv*Ye ) - log( (1-G)*So + G*Ye )
      )[(nb+1):nt] 
    invaded <- mean(delta_r) > 0
  }
  if(NA %in% rd$Ye){   
    invaded <- TRUE # if old resident goes extinct, invader establishes immediately
  }
  return(invaded)
}
	
evolve <- function(
  ni,nr,nt,nb,
  zam,wam,zsd,wsd,rho=0.82,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_o_p,phi_r,
  m0,m1,
  alpha_G0,beta_G0,
  sig_alpha_G,sig_beta_G,
  nstart=1,
  iterset=NULL,
  savefile=NULL,
  rel.tol=10^-5,
  abs.tol=0, # .Machine$double.eps^0.25,
  intsd=10,
  nsmin=10^-50
  ){
  
  require(MASS)
  cur_date <- format(Sys.Date(),"%d%b%Y")
  
  T1 <- with(Tvalues,duration[period=="T1"])
  T2 <- with(Tvalues,duration[period=="T2"])
  T3 <- with(Tvalues,duration[period=="T3"])
  # T in years
  
  tau_s <- 100  # adujstment for seeds
  tau_d <- 100	# adjustment for density
  tau_p <- 100	# adjustment for rainfall
  
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
  
  es <- data.frame(alpha_G=rep(NA,times=nr),beta_G=rep(NA,times=nr))
  es[1,] <- c(alpha_G0,beta_G0)

  for(i in 1:nr){
    if(i==1){
      rd <- ressim(alpha_G0, beta_G0) # simulate starting resident dynamics
    }
    alpha_G_inv <- es$alpha_G[i] + rnorm(1,0,sig_alpha_G)
    beta_G_inv <- es$beta_G[i] + rnorm(1,0,sig_beta_G)
    invaded <- invade(alpha_G_inv,beta_G_inv,rd)
    if(i < nr){
      if(invaded==TRUE){
        es$alpha_G[i+1] <- alpha_G_inv
        es$beta_G[i+1] <- beta_G_inv
        rd <- ressim(alpha_G_inv,beta_G_inv) # simulate new resident dynamics
      } 
      if(invaded==FALSE){
        es$alpha_G[i+1] <- es$alpha_G[i]
        es$beta_G[i+1] <- es$beta_G[i]
      } 
    }
  }
  
  outlist <- list(zw=zw,es=es)  
  
  if(is.null(savefile)){
    return(outlist)
  }
  
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
}

outlist <- evolve(
  ni,nr,nt,nb,
  zam,wam,zsd,wsd,rho=0.82,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_o_p,phi_r,
  m0,m1,
  alpha_G0,beta_G0,
  sig_alpha_G,sig_beta_G,
  nstart=1,
  iterset=NULL,
  savefile=NULL,
  rel.tol=10^-5,
  abs.tol=0, # .Machine$double.eps^0.25,
  intsd=10,
  nsmin=10^-50
)

plot(beta_G~alpha_G,data=outlist$es,type="b",col="red")
with(outlist$es, text(alpha_G[1],beta_G[1],labels="start",pos=3) )
with(outlist$es, text(alpha_G[nr],beta_G[nr],labels="end",pos=3) )

with(outlist$es[nr,],curve(plogis(alpha_G + beta_G*x),xlim=c(-2*zsd,2*zsd)))

ni;nr=1000;nt=280;nb;20
zam=zamo+mam*zsdo;zsd=zsdo*msd;
wam=wamo+mam*wsdo;wsd=wsdo*msd;
rho=0.82;
beta_p=pl$pr$beta_p[1,19,];beta_r=pl$rs$beta_r[1,19,];
sig_y_p=pl$pr$sig_y_p[1,19];sig_y_r=pl$rs$sig_y_r[1,19];
sig_o_p=pl$pr$sig_o_p[1];phi_r=pl$rs$phi_r[1];
m0=exp(pl$go$alpha_m[1,19]);m1=exp(pl$go$beta_m[1,19]);
alpha_G0=0;beta_G0=0;
sig_alpha_G=0.1;sig_beta_G=0.1;
nstart=1;
iterset=NULL;
savefile=NULL;
rel.tol=10^-5;
abs.tol=0;
intsd=10;
nsmin=10^-50
