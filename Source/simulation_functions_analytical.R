##########################################################################
### Functions for simulations of population dynamics for infinite area ###
##########################################################################

# check <- rtrunc(n=10^6,spec="nbinom",
#        mu=exp(rnorm(10^6,eta_bar_t[l,j]+eps_y_r[i,t,j],sig_s_r[i])),
#        size=rep(phi_r[i],10^6),
#        a=0)
# 
# ni <- 2
# nt <- 2
# nj <- 22
# zam <- zamo
# zsd <- zsdo
# wam <- wamo
# wsd <- wsdo
# rho <- 0.82
# iterset <- NULL
# savefile <- NULL
# rel.tol <- .Machine$double.eps^0.25/100
# intsd <- 10
# ngmin <- 0
# rho=0.82;
# Tvalues;tau_p=10^2;tau_d=10^2;tau_s=10^2;
# iterset=NULL;
# savefile=NULL;
# rel.tol=10^-5; # .Machine$double.eps^0.25
# intsd=10;
# ngmin=10^-50 # still necessary?
# i <- 1
# j <- 1
# t <- 1

# 
# set.seed(1)
# 
# nk <- 10^3
# nstart <- 10^3
# # run full code of popsim and save as outlist1
# 
# nstart <- 10^(5-3)
# # run i loop of popana using same parameters and save as outlist 2
# 
# plot(outlist1$nn[1,1,]/nk ~ outlist2$nn[1,1,],
#   xlab=expression(N[n]~infinite), ylab=expression(N[n]~finite))
# plot(outlist1$nn[2,1,]/nk ~ outlist2$nn[2,1,])
# abline(0,1,col="red")

popana <- function(pl,ni,nt,nj=22,nstart,
  zwy,zam,zsd,wam,wsd,
  Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
  iterset=NULL,
  savefile=NULL,
  rel.tol=10^-5,
  abs.tol=0, # .Machine$double.eps^0.25,
  intsd=10,
  nsmin=10^-50
  ){
  # ni = runs; nt = years; nj = species; nk = 0.1 m^2 sites
  # nk must be >1
  
  if(ni*nt*nj > 10^9 | ni*nt*nj > 10^9) stop("matrices are too large")
  
  cur_date <- format(Sys.Date(),"%d%b%Y")
  
  T1 <- with(Tvalues,duration[period=="T1"])
  T2 <- with(Tvalues,duration[period=="T2"])
  T3 <- with(Tvalues,duration[period=="T3"])
  
  go <- pl$go
  gs <- pl$gs
  pr <- pl$pr
  rs <- pl$rs
  
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
  
  ### DEFINE PARAMETERS ###
    
  iter_go <- iter_pr <- iter_rs <- iterset
  alpha_G <- go$alpha_G[iterset,]
  beta_Gz <- go$beta_Gz[iterset,]

  alpha_G <- go$alpha_G[iterset,]
  beta_Gz <- go$beta_Gz[iterset,]
  alpha_m <- go$alpha_m[iterset,]
  beta_m <- go$alpha_m[iterset,]
  
  beta_p <- beta_r <- matrix(nr=ni*nj,nc=4) 
  beta_p[] <- as.vector(pr$beta_p[iterset,,]) # iterations vary faster than species
  sig_o_p <- rep(pr$sig_o_p[iterset],times=nj)  # doesn't vary by species
  
  beta_r[] <- as.vector(rs$beta_r[iterset,,])
  phi_r <- rep(rs$phi_r[iterset],times=nj)     # doesn't vary by species
  
  zwy <- zwy[iterset,]
  
  ### RANDOM YEAR EFFECTS ###
  
  eps_y_p <- eps_y_r <- array(NA,c(ni,nt,nj))
  eps_y_p[] <- rep(zwy$eps_y_pn,times=nj) * rep(pr$sig_y_p[iterset,],each=nt)
  eps_y_r[] <- rep(zwy$eps_y_rn,times=nj) * rep(rs$sig_y_r[iterset,],each=nt)
  
  ### RAINFALL	###
  
  z <- matrix(zwy$zn*zsd + zam - log(tau_p),nr=ni,nc=nt)
  w <- matrix(zwy$wn*wsd + wam - log(tau_p),nr=ni,nc=nt)

  ### POPULATION DYNAMICS ###
  
  ns <- ng <- nn <- array(NA,c(ni,nt,nj))
  G <- m0 <- m1 <- So <- Sn <- array(NA,c(ni,nt,nj))
  ns[,1,] <- rep(nstart,each=ni)
  # for each i, replicate vector of starting densities for all species
  
  for(i in 1:ni){
    
    m0[i,,] <- exp(rep(alpha_m[i,],each=nt))
    m1[i,,] <- exp(rep(beta_m[i,],each=nt))
    So[i,,] <- exp(-m0[i,,])
    # annual DI seed mortality (T1+T2+T3)
    
    G[i,,] <- plogis(
      matrix(rep(alpha_G[i,],nt),nr=nt,nc=nj,byrow=T)
      + outer(w[i,],beta_Gz[i,],"*")
    )
    # for each i
    #   make two [t,j] matrices and sum, by:
    #     multiplying each t of w by each j of beta_Gz
    #     repeating each j of alpha_G t times
    # G not DD, so can calculate all beforehand
    # tau-adjusted rainfall included in Gmod
  }
    
  for(t in 1:nt){
      
    ng[,t,] <- G[,t,] * ns[,t,]
      
    x_t <- matrix(nr=ni*nj,nc=4)
    x_t[,1] <- 1
    x_t[,2] <- rep(z[,t],times=nj)
    x_t[,3] <- rep(z[,t]^2,times=nj)
    x_t[,4] <- log(ng[,t,]) - log(tau_d/10)
      
    pi_bar_t <- rowSums(beta_p * x_t)
    mu_bar_t <- rowSums(beta_r * x_t)
    
    al <- ns[,t,] > nsmin # alive
    pr_t <- rep(0,ni*nj) # no reproduction if dead
    pr_t[al] <- logitnormint(
      mu = (pi_bar_t + as.vector(eps_y_p[,t,]))[al],
      sigma = sig_o_p[al]
      )
    
    rs_t <- nbtmean(exp(mu_bar_t),phi_r)
      
    nn[,t,] <- ng[,t,] * pr_t * rs_t
      
    Sn[,t,] <- BHS(nn[,t,],m0[,t,],m1[,t,])
      
    if(t<nt){
      ns[,t+1,][al] <- (ns[,t,]*(1-G[,t,])*So[,t,] + Sn[,t,]*nn[,t,])[al]
      ns[,t+1,][!al] <- 0
    }
      
  } # t loop
  
  Y <- nn/ng
  
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

# popana <- function(pl,ni,nt,nj=22,nstart,
# 	zam,zsd,wam,wsd,rho=0.82,
# 	Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
# 	iterset=NULL,
#   savefile=NULL,
#   rel.tol=10^-5,
#   abs.tol=0, # .Machine$double.eps^0.25
#   intsd=10,
#   ngmin=10^-50 # still necessary?
# 	){
# 	# ni = runs; nt = years; nj = species; nk = 0.1 m^2 sites
# 	# nk must be >1
# 
# 	if(ni*nt*nj > 10^9 | ni*nt*nj > 10^9) stop("matrices are too large")
# 
# 	require(MASS)
# 
#   cur_date <- format(Sys.Date(),"%d%b%Y")
# 
# 	T1 <- with(Tvalues,duration[period=="T1"])
# 	T2 <- with(Tvalues,duration[period=="T2"])
# 	T3 <- with(Tvalues,duration[period=="T3"])
# 
# 	go <- pl$go
# 	gs <- pl$gs
# 	pr <- pl$pr
# 	rs <- pl$rs
# 
# 	BHS <- function(n,m0,m1){
# 		exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*n/(tau_s/10) )
# 	}
# 	  # original model in m^2
# 	  # DD functions use density in 0.01 m^2 = 10 x 10 cm plots
# 	  # But we want to use 0.1m^2 plots (to match scale of quadrats)
# 	  # Therefore, tau_s set to 10 instead of 100
# 
# 	logitmean <- function(mu,sigma,...){
# 	  flogit <- function(x){
# 	    plogis(x) * dnorm(x, mean=mu, sd=sigma)
# 	  }
# 	  integrate(flogit, ...)$value
# 	}
# 	  # code borrowed from logitnorm package
# 
# 	nbtmean <- function(mu,phi){
# 	  mu / ( 1 - (phi/(mu+phi))^phi )
# 	}
# 	  # mean for hurdle model
# 
# 	nbtlnmean <- function(eta,sigma,phi,...){
# 	  fnbt <- function(x){
# 	    nbtmean(mu=exp(x),phi) * dnorm(x, mean=eta, sd=sigma)
# 	  }
# 	  integrate(fnbt, ...)$value
# 	}
#   	# finite limits required to stop integration from crashing
# 
#   	# - calculate probability of each value from lognormal distribution
#   	# - each of these values produces a mean from a trunc negbin distribution
#   	# - then integrate to calculate the mean of these means
#   	# - can be done because:
#   	# sum(negbin(lognormal(mu,sig),phi))
#   	# = sum( negbin(exp(mu+sig),phi) + negbin(exp(mu,-sig),phi) )
#   	# Ref: Econometric Analysis of Count Data - Rainer Winkelmann
# 	  # (checked by simulation that still works with zero-truncation)
# 
# 	### SAMPLE ITERATIONS ###
# 
# 	if(is.null(iterset)){
# 	  iter_go <- sample(1:length(go$alpha_G_mu),ni)
# 	  iter_gs <- sample(1:length(gs$phi_g),ni)
# 	  iter_pr <- sample(1:length(pr$sig_o),ni)
# 	  iter_rs <- sample(1:length(rs$phi_r),ni) # change to phi_g later
# 	  }
# 	if(!is.null(iterset)){
# 	  iter_go <- iter_gs <- iter_pr <- iter_rs <- iterset
# 	  alpha_G <- go$alpha_G[iter_go,]
# 	  beta_Gz <- go$beta_Gz[iter_go,]
# 	  }
# 
# 	### LOAD PARAMS ###
# 
# 	alpha_G <- go$alpha_G[iter_go,]
# 	beta_Gz <- go$beta_Gz[iter_go,]
# 	alpha_m <- go$alpha_m[iter_go,]
# 	beta_m <- go$alpha_m[iter_go,]
# 
# 	sig_s_g <- gs$sig_s_g[iter_gs,]
# 	# theta_g <- gs$theta_g[iter_gs,]
# 		# zero-inflated, NOT hurdle
# 		# theta_g is prob(0)
# 	  # but ignoring here because zero-germinant sites don't contribute
# 	  # (i.e. lognormal portion of G distribution can be considered in isolation)
# 
# 	beta_p <- pr$beta_p[iter_pr,,]
# 	sig_y_p <- pr$sig_y_p[iter_pr,]
# 	sig_s_p <- pr$sig_s_p[iter_pr]  # doesn't vary by species
# 	sig_o_p <- pr$sig_o_p[iter_pr]  # doesn't vary by species
# 
# 	beta_r <- rs$beta_r[iter_rs,,]
# 	sig_y_r <- rs$sig_y_r[iter_rs,]
# 	sig_s_r <- rs$sig_s_r[iter_rs]  # doesn't vary by species
# 	phi_r <- rs$phi_r[iter_rs]          # doesn't vary by species
# 
# 	### RESCALE VARIANCE TERMS ###
# 	  # relative to first species (boin)
# 
# 	sig_y_p1 <- sig_y_p[,1]
# 	sig_y_r1 <- sig_y_r[,1]
# 	  # sds for boin
# 	scale_y_p <- as.vector(sig_y_p / sig_y_p1)
# 	scale_y_r <- as.vector(sig_y_r / sig_y_r1)
#   	# for each iter and species, i*j vector of sd scales relative to boin
#   	# (i varies faster than j)
# 
# 	### SIMULATE YEAR / SITE EFFECTS ###
# 
# 	# SIMULATE AS VECTORS
#   	# same-scale eps values for all species
#   	# i.e. simulating for boin and will re-scale later
# 	eps_y_p1 <- rnorm(ni*nt,0,rep(sig_y_p1,times=nt))
# 	eps_y_r1 <- rnorm(ni*nt,0,rep(sig_y_r1,times=nt))
# 
# 	# CONVERT TO SPECIES-SPECIFIC VALUES AS ARRAYS
# 	eps_y_p <- eps_y_r <- array(NA,c(ni,nt,nj))
# 	eps_y_p[] <- rep(eps_y_p1,times=nj) * rep(scale_y_p,each=nt)
# 	eps_y_r[] <- rep(eps_y_r1,times=nj) * rep(scale_y_r,each=nt)
#   	# relying on automatic (reverse) filling order:
#   	# rows vary fastest, then cols, then gridno
# 
# 	### SIMULATE RAINFALL	###
# 
# 	z <- w <- array(NA,c(ni,nt))
# 
# 	zw_mu <- c(zam,wam)
# 	zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
# 
# 	zw <- mvrnorm(n=ni*nt, mu=zw_mu, Sigma=zw_sig)
# 
# 	z[] <- zw[,1] - log(tau_p)
# 	# transformed log winter rainfall
# 	w[] <- zw[,2] - log(tau_p)
# 	# transformed log germination rainfall
# 
# 	### CREATE DENSITY STORAGE OBJECTS ###
# 
# 	ns <- ng <- nn <- array(NA,c(ni,nt,nj))
# 		# summed counts - permanently saved
# 	G <- m0 <- m1 <- So <- Sn <- array(NA,c(ni,nt,nj))
# 		# derived params - permanently saved
# 	ns[,1,] <- rep(nstart,each=ni)
# 		# for each i, replicate vector of starting densities for all species
# 
# 	### BEGIN CALCULATIONS ###
# 
# 	for(i in 1:ni){
# 
# 	  m0[i,,] <- exp(rep(alpha_m[i,],each=nt))
# 	  m1[i,,] <- exp(rep(beta_m[i,],each=nt))
# 	  So[i,,] <- exp(-m0[i,,])
# 			# annual DI seed mortality (T1+T2+T3)
# 
# 	  G[i,,] <- plogis(
# 	    matrix(rep(alpha_G[i,],nt),nr=nt,nc=nj,byrow=T)
# 	    + outer(w[i,],beta_Gz[i,],"*")
# 	    )
# 	    # for each i
# 	    #   make two [t,j] matrices and sum, by:
# 	    #     multiplying each t of w by each j of beta_Gz
# 	    #     repeating each j of alpha_G t times
#   	  # G not DD, so can calculate all beforehand
#   	  # tau-adjusted rainfall included in Gmod
# 
# 		for(t in 1:nt){
# 
# 			### GERMINATION ###
# 
# 			ng[i,t,] <- G[i,t,] * ns[i,t,]
# 
# 			### REPRODUCTION ###
# 
# 			xvec <- c(1,z[i,t],z[i,t]^2)
# 			  # climate data into model matrix
# 
# 			for(j in 1:nj){
# 
# 			  # running all this within a loop because integration has to be run
# 			  # one-at-a-time
# 
# 			  if(ng[i,t,j] < ngmin){
# 			    nn[i,t,j] <- 0
# 			  }
# 
# 			  if(ng[i,t,j] >= ngmin){
# 
# 			    lgmu <- log(ng[i,t,j]) - (sig_s_g[i,j]^2 / 2)
#     			  # arithmetic mean = ng[i,t,]
#     			  # logarithmic sd = sig_s_g[i,,j]
#     			  # mean of lognormal distribution = log(am) - sig^2 / 2
# 
# 			    intlo <- lgmu - intsd * sig_s_g[i,j]
# 			    inthi <- lgmu + intsd * sig_s_g[i,j]
#   			    # setting range to 10 sds to improve convergence
#   			    # (outside this range, ng=0 -> nn=0)
# 
# 			    fnn <- function(g){
# 			      # g = log(N[g]) for a given plot
# 
# 			      dg <- dnorm(g,mean=lgmu,sd=sig_s_g[i,j])
# 
# 			      nl <- length(g)
# 			      x_t <- matrix(nr=nl,nc=4)
# 			      x_t[,1:3] <- rep(xvec,each=nl)
# 			      x_t[,4] <- g - log(tau_d/10)
# 
# 			        # tau_d/10 density adjustment explained above
# 			      pi_bar_t <- beta_p[i,j,] %*% t(x_t)
# 			      eta_bar_t <- beta_r[i,j,] %*% t(x_t)
#   			      # each density (lng) has own associated world of sites
#   			      # but spatial aspects of pr(Y>0) and pr(Y|Y>0) considered independent,
#   			      # so can be simply added together
#   			      # can't average across years in this way because non-independent
# 
# 			      pr_t <- rs_t <- rep(NA,nl)
# 			      for(l in 1:nl){
# 			        pr_t[l] <- logitmean(
# 			          mu = pi_bar_t[l] + eps_y_p[i,t,j],
# 			          sigma = sqrt(sig_s_p[i]^2 + sig_o_p[i]^2),
# 			          lower=-Inf,
# 			          upper=Inf,
# 			          rel.tol=rel.tol,
# 			          abs.tol=abs.tol
# 			          )
# 
# 			        eta_t <- eta_bar_t[l] + eps_y_r[i,t,j]
# 			        sigma_t <- sig_s_r[i]
# 			        rs_t[l] <- nbtlnmean(
# 			          eta = eta_t,
# 			          sigma = sigma_t,
# 			          phi = phi_r[i],
# 			          lower = eta_t - intsd * sigma_t,
# 			          upper = eta_t + intsd * sigma_t,
# 			          rel.tol=rel.tol,
# 			          abs.tol=abs.tol
# 			          )
# 			        }
# 
# 			      lnY_t <- g + log(pr_t) + log(rs_t)
# 			      # expected log density of new seeds for each possible germinant density
# 			      # log-transforming to try and improve numerical stability
# 
# 			      return(exp(log(dg) + lnY_t))
# 			        # expected overall mean density of seeds
# 			      }
# 
# 			    nn[i,t,j] <- integrate(fnn, intlo, inthi,
# 			      rel.tol=rel.tol,abs.tol=abs.tol)$value
# 
# 			    }
# 			  } # close j loop
# 
# 			Sn[i,t,] <- BHS(nn[i,t,],m0[i,t,],m1[i,t,])
# 
# 			### BURIED SEED DENSITY ###
# 
# 			if(t<nt) ns[i,t+1,] <- So[i,t,]*(ns[i,t,]-ng[i,t,]) + Sn[i,t,]*nn[i,t,]
# 
# 			} # t loop
# 
# 		} # i loop
# 
# 	
# 	Y <- nn/ng
# 	
# 	outlist <- list(
# 	  zam=zam,zsd=zsd,
# 	  ni=ni,nt=nt,nj=nj,
# 	  z=z,w=w,
# 	  G=G,So=So,Sn=Sn,Y=Y,
# 	  ns=ns
# 	)
# 
# 	if(is.null(savefile)){
# 		return(outlist)
# 	}
# 
# 	if(!is.null(savefile)){
# 		saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
# 	}
# }

# popana <- function(pl,ni,nt,nj=22,nstart,
#   zam,zsd,wam,wsd,rho=0.82,
#   Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
#   iterset=NULL,
#   savefile=NULL,
#   rel.tol=10^-5, # .Machine$double.eps^0.25
#   abs.tol=0, 
#   intsd=10,
#   lgmin=-100 
#   ){
#   # ni = runs; nt = years; nj = species; nk = 0.1 m^2 sites 
#   # nk must be >1
#   
#   pr_f <- function(lp,lpmu,lpsd){
#     plogis(lp) * dnorm(lp,lpmu,lpsd)
#   }
#   
#   rs_f <- function(lr,lrmu,lrsd,phi_r){
#     nbtmean(mu=exp(lr),phi_r) * dnorm(lr,lrmu,lrsd)
#   }
#   
#   pr_int <- Vectorize(function(lpmu,lpsd){
#     integrate(pr_f,
#       lpmu=lpmu,lpsd=lpsd,
#       lower=-Inf,
#       upper=Inf,
#       rel.tol=rel.tol,
#       abs.tol=abs.tol)$value
#   })
#   
#   rs_int <- Vectorize(function(lrmu,lrsd,phi_r){
#     integrate(rs_f,
#       lrmu=lrmu,lrsd=lrsd,phi_r=phi_r,
#       lower=lrmu-intsd*lrsd,
#       upper=lrmu+intsd*lrsd,
#       rel.tol=rel.tol,
#       abs.tol=abs.tol)$value
#   })
#   
#   nn_f <- function(lg,lgmu,lgsd,zt,
#     beta_p1,beta_r1,
#     beta_p2,beta_r2,
#     beta_p3,beta_r3,
#     beta_p4,beta_r4,
#     eps_y_p,eps_y_r,
#     sig_a_p,sig_s_r,phi_r
#   ){
#     
#     # lg = log(N[g]) for a given plot
#     
#     dg <- dnorm(lg,mean=lgmu,sd=lgsd)
#     # tau_d/10 density adjustment explained above
#     
#     xvec <- c(1,zt,zt^2)
#     # climate data into model matrix
#     
#     nlg = length(lg)
#     x_t <- matrix(nr=nlg,nc=4)
#     x_t[,1:3] <- rep(xvec,each=nlg) 
#     x_t[,4] <- lg - log(tau_d/10)
#     # tau_d/10 density adjustment explained above
#     
#     beta_p <- cbind(beta_p1,beta_p2,beta_p3,beta_p4)
#     beta_r <- cbind(beta_r1,beta_r2,beta_r3,beta_r4)
#     
#     pi_bar_t <- beta_p %*% t(x_t)
#     eta_bar_t <- beta_r %*% t(x_t)
#     # each density (lng) has own associated world of sites
#     # but spatial aspects of pr(Y>0) and pr(Y|Y>0) considered independent,
#     # so can be simply added together
#     # can't average across years in this way because non-independent
#     
#     pr_t <- pr_int(lpmu=pi_bar_t + eps_y_p, lpsd=sig_a_p)
#     rs_t <- rs_int(lrmu=eta_bar_t + eps_y_r, lrsd=sig_s_r, phi_r=phi_r)
#     # - calculate probability of each value from lognormal dist
#     # - each of these values produces a mean from a trunc negbin dist
#     # - then integrate to calculate the mean of these means
#     # - can be done because:
#     # sum(negbin(lognormal(mu,sig),phi_r)) 
#     # = sum( negbin(exp(mu+sig),phi_r) + negbin(exp(mu,-sig),phi_r) )
#     # Ref: Econometric Analysis of Count Data - Rainer Winkelmann 
#     # (checked by simulation that still works with zero-truncation)
#     
#     lnY_t <- lg + log(pr_t) + log(rs_t)
#     # expected log density of new seeds for each possible germinant density
#     # log-transforming to try and improve numerical stability
#     
#     exp(log(dg) + lnY_t)
#     # expected overall mean density of seeds
#   } 
#   
#   nn_int <- Vectorize(function(lgmu,lgsd,zt,
#     beta_p1,beta_r1,beta_p2,beta_r2,
#     beta_p3,beta_r3,beta_p4,beta_r4,
#     eps_y_p,eps_y_r,sig_a_p,sig_s_r,phi_r){
#     integrate(nn_f,
#       lgmu=lgmu,lgsd=lgsd,zt=zt,
#       beta_p1=beta_p1,beta_r1=beta_r1,
#       beta_p2=beta_p2,beta_r2=beta_r2,
#       beta_p3=beta_p3,beta_r3=beta_r3,
#       beta_p4=beta_p4,beta_r4=beta_r4,
#       eps_y_p=eps_y_p,eps_y_r=eps_y_r,
#       sig_a_p=sig_a_p,sig_s_r=sig_s_r,phi_r=phi_r,
#       lower=lgmu-intsd*lgsd,
#       upper=lgmu+intsd*lgsd,
#       rel.tol=rel.tol,
#       abs.tol=abs.tol)$value
#   })
#   
#   if(ni*nt*nj > 10^9 | ni*nt*nj > 10^9) stop("matrices are too large")
#   
#   require(MASS) # for negative binomial distribution
#   
#   cur_date <- format(Sys.Date(),"%d%b%Y")
#   
#   T1 <- with(Tvalues,duration[period=="T1"])
#   T2 <- with(Tvalues,duration[period=="T2"])
#   T3 <- with(Tvalues,duration[period=="T3"])
#   
#   go <- pl$go
#   gs <- pl$gs
#   pr <- pl$pr
#   rs <- pl$rs
#   
#   BHS <- function(n,m0,m1){
#     exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*n/(tau_s/10) )	
#   }
#   # original model in m^2
#   # DD functions use density in 0.01 m^2 = 10 x 10 cm plots
#   # But we want to use 0.1m^2 plots (to match scale of quadrats)
#   # Therefore, tau_s set to 10 instead of 100
#   
#   nbtmean <- function(mu,phi){
#     mu / ( 1 - (phi/(mu+phi))^phi )
#   } 
#   # mean for hurdle model
#   
#   ### SAMPLE ITERATIONS ###
#   
#   if(is.null(iterset)){
#     iter_go <- sample(1:length(go$alpha_G_mu),ni)
#     iter_gs <- sample(1:length(gs$phi_g),ni)
#     iter_pr <- sample(1:length(pr$sig_o),ni)
#     iter_rs <- sample(1:length(rs$phi_r),ni) # change to phi_g later
#   }
#   if(!is.null(iterset)){
#     iter_go <- iter_gs <- iter_pr <- iter_rs <- iterset
#     alpha_G <- go$alpha_G[iter_go,]
#     beta_Gz <- go$beta_Gz[iter_go,]
#   }
#   
#   ### LOAD PARAMS ###
#   
#   alpha_G <- go$alpha_G[iter_go,]
#   beta_Gz <- go$beta_Gz[iter_go,]  
#   alpha_m <- go$alpha_m[iter_go,]
#   beta_m <- go$alpha_m[iter_go,]
#   
#   sig_s_g <- gs$sig_s_g[iter_gs,]
#   # theta_g <- gs$theta_g[iter_gs,]
#   # zero-inflated, NOT hurdle
#   # theta_g is prob(0)
#   # but ignoring here because zero-germinant sites don't contribute
#   # (i.e. lognormal portion of G distribution can be considered in isolation)
#   
#   beta_p1 <- as.vector(pr$beta_p[iter_pr,,1])
#   beta_p2 <- as.vector(pr$beta_p[iter_pr,,2])
#   beta_p3 <- as.vector(pr$beta_p[iter_pr,,3])
#   beta_p4 <- as.vector(pr$beta_p[iter_pr,,4])
#   
#   beta_r1 <- as.vector(rs$beta_r[iter_rs,,1])
#   beta_r2 <- as.vector(rs$beta_r[iter_rs,,2])
#   beta_r3 <- as.vector(rs$beta_r[iter_rs,,3])
#   beta_r4 <- as.vector(rs$beta_r[iter_rs,,4])
#   
#   sig_y_p <- pr$sig_y_p[iter_pr,]
#   sig_a_p <- rep(sqrt(pr$sig_s_p[iter_pr]^2 + pr$sig_o_p[iter_pr]^2),each=nj)
#     # overall variation in pr
#     # doesn't vary by species
#   
#   sig_y_r <- rs$sig_y_r[iter_rs,]
#   sig_s_r <- rep(rs$sig_s_r[iter_rs],each=nj)  # doesn't vary by species
#   phi_r <- rep(rs$phi[iter_rs],each=nj)        # doesn't vary by species
#   
#   ### RESCALE VARIANCE TERMS ###
#   # relative to first species (boin)
#   
#   sig_y_p1 <- sig_y_p[,1]
#   sig_y_r1 <- sig_y_r[,1]
#   # sds for boin
#   scale_y_p <- as.vector(sig_y_p / sig_y_p1)
#   scale_y_r <- as.vector(sig_y_r / sig_y_r1)	
#   # for each iter and species, i*j vector of sd scales relative to boin
#   # (i varies faster than j) 
#   
#   ### SIMULATE YEAR / SITE EFFECTS ###
#   
#   # SIMULATE AS VECTORS
#   # same-scale eps values for all species
#   # i.e. simulating for boin and will re-scale later
#   eps_y_p1 <- rnorm(ni*nt,0,rep(sig_y_p1,times=nt))
#   eps_y_r1 <- rnorm(ni*nt,0,rep(sig_y_r1,times=nt))
#   
#   # CONVERT TO SPECIES-SPECIFIC VALUES AS ARRAYS
#   eps_y_p <- eps_y_r <- array(NA,c(ni,nt,nj))
#   eps_y_p[] <- rep(eps_y_p1,times=nj) * rep(scale_y_p,each=nt)
#   eps_y_r[] <- rep(eps_y_r1,times=nj) * rep(scale_y_r,each=nt)
#   # relying on automatic (reverse) filling order:
#   # rows vary fastest, then cols, then gridno
#   
#   ### SIMULATE RAINFALL	###
#   
#   z <- w <- array(NA,c(ni,nt))
#   
#   zw_mu <- c(zam,wam)
#   zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
#   
#   zw <- mvrnorm(n=ni*nt, mu=zw_mu, Sigma=zw_sig)
#   
#   z[] <- zw[,1] - log(tau_p)
#   # transformed log winter rainfall
#   w[] <- zw[,2] - log(tau_p)
#   # transformed log germination rainfall
#   
#   ### CREATE DENSITY STORAGE OBJECTS ###
#   
#   ns <- ng <- nn <- array(NA,c(ni,nt,nj))
#   # summed counts - permanently saved
#   G <- m0 <- m1 <- So <- Sn <- array(NA,c(ni,nt,nj))
#   # derived params - permanently saved
#   ns[,1,] <- rep(nstart,each=ni)
#   # for each i, replicate vector of starting densities for all species
#   
#   ### BEGIN CALCULATIONS ###
#     
#   for(i in 1:ni){
#     m0[i,,] <- exp(rep(alpha_m[i,],each=nt))
#     m1[i,,] <- exp(rep(beta_m[i,],each=nt))
#     So[i,,] <- exp(-m0[i,,])
#     # annual DI seed mortality (T1+T2+T3)
#     G[i,,] <- plogis(
#       matrix(rep(alpha_G[i,],nt),nr=nt,nc=nj,byrow=T)
#       + outer(w[i,],beta_Gz[i,],"*")	
#     )  	  
#     # for each i
#     #   make two [t,j] matrices and sum, by:
#     #     multiplying each t of w by each j of beta_Gz
#     #     repeating each j of alpha_G t times
#     # G not DD, so can calculate all beforehand
#     # tau-adjusted rainfall included in Gmod
#   }
#   
#   for(t in 1:nt){
#     
#     ### GERMINATION ###
#     ng[,t,] <- G[,t,] * ns[,t,]
#       
#     ### REPRODUCTION ###
#     lgmu <- log(ng[,t,]) - (sig_s_g^2 / 2)
#     # arithmetic mean = ng[i,t,]
#     # logarithmic sd = sig_s_g[i,,j]
#     # mean of lognormal distribution = log(am) - sig^2 / 2
#       
#     lgk <- as.vector(lgmu) > lgmin
#     
#     nn[,t,][!lgk] <- 0
#     
#     nn[,t,][lgk] <- nn_int(lgmu=lgmu[lgk],lgsd=sig_s_g[lgk],zt=z[,t],
#       beta_p1[lgk],beta_r1[lgk],beta_p2[lgk],beta_r2[lgk],
#       beta_p3[lgk],beta_r3[lgk],beta_p4[lgk],beta_r4[lgk],
#       as.vector(eps_y_p[,t,])[lgk],as.vector(eps_y_r[,t,])[lgk],
#       sig_a_p[lgk],sig_s_r[lgk],phi_r[lgk]
#     )
#       # running all this within a loop because integration has to be run 
#       # one-at-a-time
#           
#       # setting range to +/-intsd to improve convergence
#       # (outside this range, ng=0 -> nn=0)  
#       # reproduction terms get first sd in g, then sd in themselves
#           
#     Sn[,t,] <- BHS(nn[,t,],m0[,t,],m1[,t,])
#       
#     ### BURIED SEED DENSITY ###
#       
#     if(t<nt) ns[,t+1,] <- So[,t,]*(ns[,t,]-ng[,t,]) + Sn[,t,]*nn[,t,]
#       
#   } # t loop
# 
#   Y <- nn/ng
#   
#   outlist <- list(
#     zam=zam,zsd=zsd,
#     ni=ni,nt=nt,nj=nj,
#     z=z,w=w,
#     G=G,So=So,Sn=Sn,Y=Y,
#     ns=ns
#   )
#   
#   if(is.null(savefile)){
#     return(outlist)
#   }
#   
#   if(!is.null(savefile)){
#     saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
#   }
#   
# }

popinv <- function(
  # pli,plr,iterseti,itersetr,
  zwy,alpha_G,beta_Gz,
  Y,Sn,So,
  ni,nt,nj=22,tmin=1,
  savefile=NULL
  ){

  cur_date <- format(Sys.Date(),"%d%b%Y")
  
  z <- matrix(zwy$zn*zsd + zam - log(tau_p),nr=ni,nc=nt)
  w <- matrix(zwy$wn*wsd + wam - log(tau_p),nr=ni,nc=nt)
  
  G <- array(dim=c(ni,nt,nj))
  for(i in 1:ni){
    G[i,,] <- plogis(
      matrix(rep(alpha_G[i,],nt),nr=nt,nc=nj,byrow=T)
      + outer(w[i,],beta_Gz[i,],"*")
      )  	  
  }

  ri <- log(G*Y*Sn + (1-G)*So)
  
  outlist <- list(rbari=apply(ri[,tmin:nt,],c(1,3),mean))
    # can calculate Y, Ye, and So from input (plr parameters)

  if(is.null(savefile)){
    return(outlist)
  }
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
}
  # assumes that both temporal and spatial terms are same for both strategies
  # e.g. invader has x times lower density than resident in every plot
  
