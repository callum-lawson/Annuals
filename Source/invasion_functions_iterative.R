### Calculate ES alpha_G and beta_G by iteratively perturbing parameters  ###
### and attempting re-invasion                                            ###

# - generate invader parameters
# - calculate whether invader invades resident, return Y/N
# - if Y, take new params; if N, keep old ones
# - loop this over n invasions
# - loop this over n scenarios (clim, params / species)

cur_date <- format(Sys.Date(),"%d%b%Y")

T1 <- with(Tvalues,duration[period=="T1"])
T2 <- with(Tvalues,duration[period=="T2"])
T3 <- with(Tvalues,duration[period=="T3"])

tau_p <- 10^2
tau_d <- 10^2
tau_s <- 10^2
nstart <- 1
nburn <- 1 
nt <- 1 # time series length basically = number of replicates
eps_alpha <- 1
eps_beta <- 1
savefile <- TRUE
  
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
  
invade_f <- function(
  zam,wam,tau_p,rho,zsd,wsd,sig_y_p,sig_r_p,
  m0,tau_d,beta_p,beta_r,sig_o_p,phi_r,m1,
  alpha_G_res,beta_G_res,alpha_G_inv,beta_G_inv
  )  {
  
  zw_mu <- c(zam,wam) - log(tau_p)
  zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
  zw <- mvrnorm(n=nt, mu=zw_mu, Sigma=zw_sig)
  eps_y_p <- rnorm(nt,0,sig_y_p)
  eps_y_p <- rnorm(nt,0,sig_r_p)
  So <- rep(exp(-m0),nt)
  Gres <- plogis(alpha_G_res,beta_G_res*w)
  ns[1] <- nstart
  t <- 1
  
  while(ns[t] > nsmin){
    ng[t] <- Gres[t] * ns[t]
    x_t <- c(1,z[t],z[t]^2,log(ng[t])-log(tau_d/10))
    pi_bar_t <- sum(beta_p * x_t)
    mu_bar_t <- sum(beta_r * x_t)
    pr_t <- logitnormint(mu=pi_bar_t+eps_y_p[t], sigma=sig_o_p)
    rs_t <- nbtmean(exp(mu_bar_t),phi_r)
    nn[t] <- ng[t] * pr_t * rs_t
    Sn[t] <- BHS(nn[t],m0[t],m1[t])
    ns[t+1] <- (ns[t]*(1-G[t])*So[t] + Sn[t]*nn[t])
    t <- t + 1
  }
  
  if(ns[t] > nsmin){   
    # t = final value at which loop stopped
    Y <- nn/ng
    Ginv <- plogis(alpha_G_inv,beta_G_inv*w)
    invade <- mean(log(Ginv*Y*Sn + (1-Ginv)*So)) > 0
  }
  if(ns[t] <= nsmin){   
    invade <- NA
  }
  return(invade)
}
	
  
  outlist <- list(
    zam=zam,zsd=zsd,
    ni=ni,nt=nt,nj=nj,
    G=G,So=So,Sn=Sn,Y=Y,
    ns=ns
  )
  
  if(is.null(savefile)){
    return(outlist)
  }
  
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
}