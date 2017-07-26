##########################################################################
### Functions for simulations of population dynamics for infinite area ###
##########################################################################

# TODO
# - check predictions
# - tau adjustments

# check <- rtrunc(n=10^6,spec="nbinom",
#        mu=exp(rnorm(10^6,eta_bar_t[l,j]+eps_y_r[i,t,j],sig_s_r[i])),
#        size=rep(phi[i],10^6),
#        a=0)

ni <- 2
nt <- 20
nj <- 22
zam <- zamo
zsd <- zsdo
wam <- wamo
wsd <- wsdo
rho <- 0.82
nstart <- 1
iterset <- NULL
savefile <- NULL

set.seed(1)
hi1 <- popsim(pl,ni,nt,nj=22,nstart,zam,zsd,wam,wsd,rho=0.82,Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,iterset=NULL,savefile=NULL,abs.tol=.Machine$double.eps^0.25)
set.seed(1)
hi2 <- popsim(pl,ni,nt,nj=22,nstart,zam,zsd,wam,wsd,rho=0.82,Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,iterset=NULL,savefile=NULL,abs.tol=.Machine$double.eps^0.25/10000)

popsim <- function(pl,ni,nt,nj=22,nstart,
	zam,zsd,wam,wsd,rho=0.82,
	Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
	iterset=NULL,savefile=NULL,abs.tol=.Machine$double.eps^0.25
	){
	# dstart = vector of starting seed totals
	# (can adjust area or density independently)
	# ni = runs; nt = years; nj = species; nk = 0.1 m^2 sites 
	# nk must be >1

	if(ni*nt*nj > 10^9 | ni*nt*nj > 10^9) stop("matrices are too large")

	require(MASS)

  cur_date <- format(Sys.Date(),"%d%b%Y")
  
	T1 <- with(Tvalues,duration[period=="T1"])
	T2 <- with(Tvalues,duration[period=="T2"])
	T3 <- with(Tvalues,duration[period=="T3"])

	go <- pl$go
	gs <- pl$gs
	pr <- pl$pr
	rs <- pl$rs

	BHS <- function(n,m0,m1){
		exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*n/tau_s )	
	}
	  # adjustment for area made when entering n (below)
	  # (could have instead been made in BHS function)
  	# nk/10: number of 0.1m^2 plots -> number of 1 m^2 plots 
  	# tau_s on coefficients: density in 0.01 m^2 = 10 x 10 cm plots

	logitmean <- function(mu,sigma){
	  flogit <- function(x){
	    plogis(x) * dnorm(x, mean=mu, sd=sigma)
	  }
	  integrate(flogit, -Inf, Inf, abs.tol=abs.tol)$value
	}
	  # code borrowed from logitnorm package
	
	nbtmean <- function(mu,phi){
	  mu / ( 1 - (phi/(mu+phi))^phi )
	} 
	  # mean for hurdle model
	
	nbtlnmean <- function(eta,sigma,phi){
	  fnbt <- function(x){
	    nbtmean(mu=exp(x),phi) * dnorm(x, mean=eta, sd=sigma)
	  }
	  integrate(fnbt, -20, 20, abs.tol=abs.tol)$value
	}
  	# finite limits required to stop integration from crashing
  	
  	# - calculate probability of each value from lognormal distribution
  	# - each of these values produces a mean from a trunc negbin distribution
  	# - then integrate to calculate the mean of these means
  	# - can be done because:
  	# sum(negbin(lognormal(mu,sig),phi)) 
  	# = sum( negbin(exp(mu+sig),phi) + negbin(exp(mu,-sig),phi) )
  	# Ref: Econometric Analysis of Count Data - Rainer Winkelmann 
	
	### SAMPLE ITERATIONS ###

	if(is.null(iterset)==T){
	  iter_go <- sample(1:length(go$alpha_G_mu),ni)
	  iter_gs <- sample(1:length(gs$phi_g),ni)
	  iter_pr <- sample(1:length(pr$sig_o),ni)
	  iter_rs <- sample(1:length(rs$phi_r),ni) # change to phi_g later
	  }
	if(is.null(iterset)==F){
	  iter_go <- iter_gs <- iter_pr <- iter_rs <- iterset
	  }
	  
	### LOAD PARAMS ###

	alpha_G <- go$alpha_G[iter_go,]
	beta_Gz <- go$beta_Gz[iter_go,]
	alpha_m <- go$alpha_m[iter_go,]
	beta_m <- go$alpha_m[iter_go,]
	
	sig_s_g <- gs$sig_s_g[iter_gs,]
	# theta_g <- gs$theta_g[iter_gs,]
		# zero-inflated, NOT hurdle
		# theta_g is prob(0)
	  # but ignoring here because zero-germinant sites don't contribute
	  # (i.e. lognormal portion of G distribution can be considered in isolation)

	beta_p <- pr$beta_p[iter_pr,,]
	sig_y_p <- pr$sig_y_p[iter_pr,]
	sig_s_p <- pr$sig_s_p[iter_pr]	 # doesn't vary by species
	sig_o_p <- pr$sig_o_p[iter_pr]   # doesn't vary by species

	beta_r <- rs$beta_r[iter_rs,,]
	sig_y_r <- rs$sig_y_r[iter_rs,]
	sig_s_r <- rs$sig_s_r[iter_rs] 	# now doesn't vary by species
	phi <- rs$phi[iter_rs]

	### RESCALE VARIANCE TERMS ###
	# relative to first species (boin)
	
	sig_y_p1 <- sig_y_p[,1]
	sig_y_r1 <- sig_y_r[,1]
		# sds for boin
	scale_y_p <- as.vector(sig_y_p / sig_y_p1)
	scale_y_r <- as.vector(sig_y_r / sig_y_r1)	
	# for each iter and species, i*j vector of sd scales relative to boin
	# (i varies faster than j) 

	### SIMULATE YEAR / SITE EFFECTS ###

	# SIMULATE AS VECTORS
		# same-scale eps values for all species
		# i.e. simulating for boin and will re-scale later
	eps_y_p1 <- rnorm(ni*nt,0,rep(sig_y_p1,times=nt))
	eps_y_r1 <- rnorm(ni*nt,0,rep(sig_y_r1,times=nt))
	# CONVERT TO SPECIES-SPECIFIC VALUES AS ARRAYS
	eps_y_p <- eps_y_r <- array(NA,c(ni,nt,nj))
	eps_y_p[] <- rep(eps_y_p1,times=nj) * rep(scale_y_p,each=nt)
	eps_y_r[] <- rep(eps_y_r1,times=nj) * rep(scale_y_r,each=nt)
		# relying on automatic (reverse) filling order:
		# rows vary fastest, then cols, then gridno

	### SIMULATE RAINFALL	###

	z <- w <- array(NA,c(ni,nt))

	zw_mu <- c(zam,wam)
	zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)

	zw <- mvrnorm(n=ni*nt, mu=zw_mu, Sigma=zw_sig)

	z[] <- zw[,1] - log(tau_p)
		# transformed log winter rainfall
	w[] <- zw[,2] - log(tau_p)
		# transformed log germination rainfall

	### CREATE DENSITY STORAGE OBJECTS ###

	ns <- ng <- no <- nn <- nnb <- array(NA,c(ni,nt,nj))
		# summed counts - permanently saved
	G <- m0 <- m1 <- So <- Sn <- array(NA,c(ni,nt,nj))
		# derived params - permanently saved

	# pi_bar_t <- pi_t <- eta_bar_t <- eta_t <- eps_o_p_t <- array(NA,c(nl,nj))
	# 	# site-level counts and params - replaced every new (i,t)

	ns[,1,] <- rep(nstart,each=ni)
		# for each i, replicate vector of starting densities for all species

	### BEGIN CALCULATIONS ###

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
	  
		for(t in 1:nt){
		  
			### GERMINATION ###
		
			ng[i,t,] <- G[i,t,]*ns[i,t,]

			### OLD SEED SURVIVAL ###

			no[i,t,] <- So[i,t,]*(ns[i,t,]-ng[i,t,])
				# no defined at START of year, so same for initial year too

			### REPRODUCTION ###

			xvec <- c(1,z[i,t],z[i,t]^2)
			  # climate data into model matrix
			
			for(j in 1:nj){
			
				# running all this within a loop because integration has to be run 
			  # one-at-a-time
			  
			  #### integrate wrt lg
			  
			  fnn <- function(g){
			    # g = log(N[g]) for a given plot
			    
			    lgmu <- log(ng[i,t,j]) - (sig_s_g[i,j]^2 / 2)
			    # arithmetic mean = ng[i,t,]
			    # logarithmic sd = sig_s_g[i,,j]
			    # mean of lognormal distribution = log(am) - sig^2 / 2
			    
			    dg <- dnorm(g,mean=lgmu,sd=sig_s_g[i,j])
			    
			    nl <- length(g)
			    x_t <- matrix(nr=nl,nc=4)
			    x_t[,1:3] <- rep(xvec,each=nl) 
			    x_t[,4] <- g - log(tau_d)
			    pi_bar_t <- beta_p[i,j,] %*% t(x_t)
			    eta_bar_t <- beta_r[i,j,] %*% t(x_t)
			    # each density (lng) has own associated world of sites
			    # but spatial aspects of pr(Y>0) and pr(Y|Y>0) considered independent,
			    # so can be simply added together
			    # can't average across years in this way because non-independent
			   
			    # ***DOES lg NEED TO BE ADJUSTED FOR TAU?***
			    
			    pr_t <- rs_t <- rep(NA,nl)
			    for(l in 1:nl){
			      pr_t[l] <- logitmean(
			        mu = pi_bar_t[l] + eps_y_p[i,t,j], 
			        sigma = sig_s_p[i] + sig_o_p[i]
			      )
			      rs_t[l] <- nbtlnmean(
			        eta = eta_bar_t[l] + eps_y_r[i,t,j], 
			        sigma = sig_s_r[i], 
			        phi = phi[i]
			      )
			    }
		
			    nY_t <- ng[i,t,j] * pr_t * rs_t
			    # expected final density of seeds for each possible germinant density
			    
			    return(dg * nY_t) 
			    # expected overall mean density of seeds
			  }
			    
			  nn[i,t,j] <- integrate(fnn, -Inf, Inf, abs.tol=abs.tol)$value

				} # close j loop

			### NEW SEED SURVIVAL ###

			Sn[i,t,] <- BHS(nn[i,t,],m0[i,t,],m1[i,t,])
			nnb[i,t,] <- Sn[i,t,]*nn[i,t,]
		
			### BURIED SEED DENSITY ###

			if(t<nt) ns[i,t+1,] <- nnb[i,t,] + no[i,t,]

			} # t loop

		} # i loop

	outlist <- list(
	    zam=zam,zsd=zsd,
			ni=ni,nt=nt,nj=nj,
	    z=z,w=w,
			G=G,So=So,Sn=Sn,
			ns=ns,ng=ng,no=no,nn=nn,nnb=nnb
			)

	if(is.null(savefile)==T){
		return(outlist)
		}
	if(is.null(savefile)==F){
		saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
		}

	}
