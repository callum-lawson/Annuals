##########################################################################
### Functions for simulations of population dynamics for infinite area ###
##########################################################################

# treat as different species (with different parameters), but link densities?
# (i.e. j indexs competitors, not species)

# every run uses new param samples AND year/site/obs effects
# make sites bigger instead of using larger number of sites...?

ni <- 2
nt <- 10
nj <- 22
nk <- 5
zam <- zam
zcv <- zcv
wam <- wam
wcv <- wcv
rho <- 0.82
nstart <- 1
iterset <- NULL
savefile <- NULL
progress <- F

# hi <- popsim(pl,ni,nt,nj=22,nk,nstart,
#   zam,zcv,wam,wcv,rho,
#   Tvalues,tau_p,tau_d,tau_s,
#   iterset,savefile,progress)

# Modelling densities per 10 x 10cm, or 1 x 1m?

popsim <- function(pl,ni,nt,nj=22,nstart,
	zam,zsd,wam,wsd,rho=0.82,
	Tvalues,tau_p=10^2,tau_d=10^2,tau_s=10^2,
	iterset=NULL,savefile=NULL,progress=F
	){
	# dstart = vector of starting seed totals
	# (can adjust area or density independently)
	# ni = runs; nt = years; nj = species; nk = 0.1 m^2 sites 
	# nk must be >1

	if(ni*nt*nj > 10^9 | ni*nt*nj > 10^9) stop("matrices are too large")

	require(truncdist)
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
	if(progress==T) pb <- txtProgressBar(min=0,max=ni,style=3)
  	# nk/10: number of 0.1m^2 plots -> number of 1 m^2 plots 
  	# tau_s on coefficients: density in 0.01 m^2 = 10 x 10 cm plots

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
	theta_g <- gs$theta_g[iter_gs,]
		# zero-inflated, NOT hurdle
		# theta_g is prob(0)
	  # possible to get all zero sites, but v. unlikely for large sites,
	  # so ignoring

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
	sig_s_g1 <- sig_s_g[,1]
	sig_s_p1 <- sig_s_p[]					# doesn't vary by species
	sig_s_r1 <- sig_s_r[]					# now doesn't vary by species	
		# sds for boin

	scale_y_p <- as.vector(sig_y_p / sig_y_p1)
	scale_y_r <- as.vector(sig_y_r / sig_y_r1)	
	scale_s_g <- as.vector(sig_s_g / sig_s_g1)
	# for each iter and species, i*j vector of sd scales relative to boin
	# (i varies faster than j) 
	scale_s_p <- as.vector(rep(1,ni*nj))
	scale_s_r <- as.vector(rep(1,ni*nj))
	# doesn't vary by species

	### SIMULATE YEAR / SITE EFFECTS ###

	# # SIMULATE AS VECTORS
	# 	# same-scale eps values for all species
	# 	# i.e. simulating for boin and will re-scale later
	# eps_y_p1 <- rnorm(ni*nt,0,rep(sig_y_p1,times=nt))
	# eps_y_r1 <- rnorm(ni*nt,0,rep(sig_y_r1,times=nt))
	# eps_s_g1 <- rnorm(ni*nk,0,rep(sig_s_g1,times=nk))
	# eps_s_p1 <- rnorm(ni*nk,0,rep(sig_s_p1,times=nk))
	# eps_s_r1 <- rnorm(ni*nk,0,rep(sig_s_r1,times=nk))
	# 	# using same sites, but doesn't matter
	# 
	# # CONVERT TO SPECIES-SPECIFIC VALUES AS ARRAYS
	# 
	# eps_y_p <- eps_y_r <- array(NA,c(ni,nt,nj))
	# eps_s_g <- eps_s_p <- eps_s_r <- array(NA,c(ni,nk,nj))
	# eps_y_p[] <- rep(eps_y_p1,times=nj) * rep(scale_y_p,each=nt)
	# eps_y_r[] <- rep(eps_y_r1,times=nj) * rep(scale_y_r,each=nt)
	# eps_s_g[] <- rep(eps_s_g1,times=nj) * rep(scale_s_g,each=nk)
	# eps_s_p[] <- rep(eps_s_p1,times=nj) * rep(scale_s_p,each=nk)
	# eps_s_r[] <- rep(eps_s_r1,times=nj) * rep(scale_s_r,each=nk)
	# 	# relying on automatic (reverse) filling order: 
	# 	# rows vary fastest, then cols, then gridno 
	# 
	# zsites <- array(NA,c(ni,nj,nk))
	#   # dim(zsites) != dim(eps_s_g)
	# zsites[] <- rbinom(ni*nj*nk,size=1,prob=rep(theta_g,times=nk))
	#   # theta_g has dim c(ni,nj)
	#   # so rep(as.vector(theta_g),times=nk) -> (i then j), nk times
	#   # theta = prob of zero
	# zsites <- aperm(zsites, c(1,3,2))
	#   # [i,j,k] -> [i,k,j] (to match eps_s_g)
	# eps_s_g[zsites==1] <- -Inf
	#   # p(germ) = exp(-Inf) = 0

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

	# ng_t <- nr_t <- pi_bar_t <- pi_t <- eta_bar_t <- eta_t <- eps_o_p_t <- array(NA,c(nk,nj))
	# 	# site-level counts and params - replaced every new (i,t)

	# x_t <- matrix(1,nr=nk,nc=4)
	# 	# makes intercept already; other params will be filled in later
	# 	# prcp replaced every (i,t)
	# 	# dens replaced every (i,t,j) [for each k]

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

		  # G[i,t,] <- plogis(
		  #   alpha_G[i,]
		  #   + beta_Gz[i,]*w[i,t]
		  #   + beta_Gd[i,]*log(ns[i,t,]/tottarea)
		  #   )
		  
			### GERMINATION ###
		
			ng[i,t,] <- G[i,t,]*ns[i,t,]

			### OLD SEED SURVIVAL ###

			no[i,t,] <- So[i,t,]*(ns[i,t,]-ng[i,t,])
				# no defined at START of year, so same for initial year too

			### REPRODUCTION ###

			x_t[,2] <- z[i,t]
			x_t[,3] <- z[i,t]^2
			  # climate data into model matrix
			
			#  eps_o_p_t[] <- rnorm(nk*nj,0,sig_o_p[i])
			   # MIGHT NOT NEED THIS ANY MORE
			
			sig_o_p[i]
			
			# arithmetic mean = ng[i,t,]
			# logarithmic sd = eps_s_g[i,,j]

			for(j in 1:nj){
			
				# running all this within a loop to avoid many calculations on zeroes
				# when species has gone extinct

			  # STILL NECESSARY?
			  
			  lgmu[i,t,] <- log(ng[i,t,]) - (eps_s_g[i,,j]^2)
			  # mean of lognormal distribution = log(am) - sig^2 / 2
			  
			  lngseq <- seq(-5,5,length.out=100)
			  dg <- dnorm(lngseq,mean=lgmu[i,t,],sd=eps_s_g[i,,j])
			  
			  # nseq REPLACING nk?
			  
			  x_t[,4] <- lngseq - log(tau_d) # check this - double-adjusting for tau?
			  
			  pi_bar_t[,j] <- beta_p[i,j,] %*% t(x_t)
			  pi_t[,j] <- pi_bar_t[,j] + eps_y_p[i,t,j] + eps_s_p[i,,j] + eps_o_p_t[,j]
			  
			  nbtmean <- function(mu,phi){
			    mu / ( 1 - (phi/(mu+phi))^phi )
			  }
			  eta_bar_t[,j] <- beta_r[i,j,] %*% t(x_t)
			  eta_t[,j] <- eta_bar_t[,j] + eps_y_r[i,t,j] + eps_s_r[i,,j]
			  rmu[,j] <- nbtmean(exp(eta_t[,j]),rs$phi[i])
			  # does this fully account for NLA?
			  # what to do with site variable here and above?
			  
			  nY_t[,j] <- plogis(pi_t[,j])*ng_t[,j]*rmu[,j] 
			  # expected final density of seeds for each possible germinant density
			  
			  nn[i,t,j] <- sum(nY_t[,j] * dg) / sum(dg)

				} # close j loop

			### NEW SEED SURVIVAL ###

			Sn[i,t,] <- BHS(nn[i,t,],m0[i,t,],m1[i,t,])
			nnb[i,t,] <- Sn[i,t,]*nn[i,t,]
		
			### BURIED SEED DENSITY ###

			if(t<nt) ns[i,t+1,] <- nnb[i,t,] + no[i,t,]

			} # t loop

   		if(progress==T) setTxtProgressBar(pb, i)

		} # i loop

	outlist <- list(
	    zam=zam,zsd=zsd,
			ni=ni,nt=nt,nj=nj,nk=nk,
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
