### Calculate ES alpha_G and beta_G by iteratively perturbing parameters  ###
### and attempting re-invasion                                            ###

BHS <- function(n,m0,m1,T3=0.6794521,tau_s=100){
  exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*n/(tau_s/10) )
}

RICKERS <- function(n,m0,m1,T3=0.6794521,tau_s=100) {
  exp(-(m0+m1*n/(tau_s/10))*T3)
}
# original model in m^2
# DD functions use density in 0.01 m^2 = 10 x 10 cm plots
# But we want to use 0.1m^2 plots (to match scale of quadrats)
# Therefore, tau_s set to 10 instead of 100
 
logitnorm <- function(x,mu,sigma){
  plogis(x) * dnorm(x, mean=mu, sd=sigma)
}
# code borrowed from logitnorm package

logitmean <- function(mu,sigma){
	 integrate(logitnorm, mu=mu, sigma=sigma,
	   lower=-Inf,
	   upper=Inf
	   )$value
}

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

nbtnorm <- function(x,eta,sigma,phi){
  nbtmean(mu=exp(x),phi=phi) * dnorm(x, mean=eta, sd=sigma)
}

nbtlnmean <- function(eta,sigma,phi,intsd=10){
	 integrate(nbtnorm,
	   eta=eta,sigma=sigma,phi=phi,
	   lower=eta-intsd*sigma,
	   upper=eta+intsd*sigma
	   )$value
}
  # finite limits required to stop integration from crashing

  # - calculate probability of each value from lognormal distribution
  # - each of these values produces a mean from a trunc negbin distribution
  # - then integrate to calculate the mean of these means
  # - can be done because:
  # sum(negbin(lognormal(mu,sig),phi))
  # = sum( negbin(exp(mu+sig),phi) + negbin(exp(mu,-sig),phi) )
  # Ref: Econometric Analysis of Count Data - Rainer Winkelmann
	# (checked by simulation that still works with zero-truncation)

fnn <- function(g,
  lgmu,x_z_t,
  beta_p,beta_r,
  eps_y_p_t,eps_y_r_t,
  sig_s_g,sig_s_p,sig_s_r,
  sig_o_p,phi_r,
  tau_d=100
  ){
  # g = log(N[g]) for a given plot
  
  dg <- dnorm(g,mean=lgmu,sd=sig_s_g)
  
  nl <- length(g)
  x_t <- matrix(nr=nl,nc=4)
  x_t[,1:3] <- rep(x_z_t,each=nl)
  x_t[,4] <- g - log(tau_d/10)
  # tau_d/10 density adjustment explained above
  
  pi_bar_t  <- beta_p %*% t(x_t) + eps_y_p_t
  eta_bar_t <- beta_r %*% t(x_t) + eps_y_r_t
  # each density (lng) has own associated world of sites
  # but spatial aspects of pr(Y>0) and pr(Y|Y>0) considered independent,
  # so can be simply added together
  # can't average across years in this way because non-independent
  
  pr_t <- rs_t <- rep(NA,nl)
  
  for(l in 1:nl){
    pr_t[l] <- logitmean(
      mu = pi_bar_t[l],
      sigma = sqrt(sig_s_p^2 + sig_o_p^2)
    )
    
    eta_t <- eta_bar_t[l]
    rs_t[l] <- nbtlnmean(
      eta = eta_t,
      sigma = sig_s_r,
      phi = phi_r
    )
  } # close l loop
  
  lnY_t <- g + log(pr_t) + log(rs_t)
  # expected log density of new seeds for each possible germinant density
  # log-transforming to try and improve numerical stability
  
  return(exp(log(dg) + lnY_t))
  # expected overall mean density of seeds
  
} # close g function

fixG <- function(w,a,b){
  plogis(a+b*w)
}

ressim <- function(w,x_z,am,bm,# as,bs,abr,
                   beta_p,beta_r,
                   eps_y_p,eps_y_r,
                   sig_s_g,sig_s_p,sig_s_r,
                   sig_o_p,phi_r,
                   So,m0,m1,
                   nt,nk,nsmin,ngmin,
                   DDFUN,
                   Sg,
                   nc=5,   # n consecutive t that ns must be < nsmin
                   nstart=1,
                   intsd=10,
                   tau_d=100
                   # full=FALSE
                   ){

  ns <- ng <- nn <- Ye <- rep(NA,nt)
  #Gres <- coaG(w,am,bm,as,bs,abr)
  Gres <- fixG(w,am,bm)
  ns[1] <- nstart
  t <- 1
  
  while(t <= nt
    & ifelse(t < nc, TRUE, FALSE %in% (ns[(t-(nc-1)):t] < nsmin))
    ){
    
    ng[t] <- Sg * Gres[t] * ns[t]
    
    if(ng[t] >= ngmin){
      
      if(nk==0){
        
        x_t <- c(x_z[t,],log(ng[t])-log(tau_d/10))
        pi_bar_t <- sum(beta_p * x_t) + eps_y_p[t]
        mu_bar_t <- sum(beta_r * x_t) + eps_y_r[t]
        pr_t <- logitnormint(mu=pi_bar_t,sigma=sig_o_p)
        rs_t <- nbtmean(exp(mu_bar_t),phi_r)
        nn[t] <- ng[t] * pr_t * rs_t
        
      } # nk==0
      
      if(nk==Inf){
        
        lgmu <- log(ng[t]) - (sig_s_g^2 / 2)
        # arithmetic mean = ng[i,t,]
        # logarithmic sd = sig_s_g[i,,j]
        # mean of lognormal distribution = log(am) - sig^2 / 2
        
        intlo <- lgmu - intsd * sig_s_g
        inthi <- lgmu + intsd * sig_s_g
        # setting range to 10 sds to improve convergence
        # (outside this range, ng=0 -> nn=0)
        
        nn[t] <- integrate(fnn, lower=intlo, upper=inthi,
          lgmu=lgmu,x_z_t=x_z[t,],
          beta_p=beta_p,beta_r=beta_r,
          eps_y_p_t=eps_y_p[t],eps_y_r_t=eps_y_r[t],
          sig_s_g=sig_s_g,sig_s_p=sig_s_p,sig_s_r=sig_s_r,
          sig_o_p=sig_o_p,phi_r=phi_r
          )$value
        
      } # nk==Inf
      
    } # close if function
    
    if(ng[t] < ngmin){
      nn[t] <- 0
    }
    
    Ye[t] <- nn[t] * DDFUN(nn[t],m0,m1) / ng[t]
    if(t<nt) ns[t+1] <- ns[t] * ( (1-Gres[t])*So + Gres[t]*Ye[t] )
    t <- t + 1
  } # close t loop
  
  return(data.frame(Gres=Gres,Ye=Ye))
  # if(full==FALSE) return(data.frame(Gres=Gres,Ye=Ye))
  # if(full==TRUE) return(data.frame(Gres=Gres,Ye=Ye,ns=ns))
}

invade <- function(w,ami,bmi,Gres,Ye,So,nt,nb){ # asi,bsi,abri,
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
  nr,nt,nb,nk,
  zam,wam,zsd,wsd,rho,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_s_g,sig_s_p,sig_s_r,
  sig_o_p,phi_r,
  m0,m1,
  am0,bm0,
  # as0,bs0,
  # abr0,
  DDFUN,
  Sg,
  smut_m,# smut_s=0.1,smut_r=0.1,
  nsmin,
  ngmin,
  lastonly,
  tau_p=100
  ){
  
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
  
  es <- data.frame(am=rep(NA,times=nr),
                   bm=rep(NA,times=nr)
                   # as=rep(NA,times=nr),
                   # bs=rep(NA,times=nr),
                   # abr=rep(NA,times=nr)
                   )
  es[1,] <- c(am0,bm0) # ,as0,bs0,abr0)

  for(i in 1:nr){
    
    if(i==1){
      rd <- with(es[i,], ressim(zw[,2],x_z,am,bm,# as,bs,abr,
                                beta_p,beta_r,
                                eps_y_p,eps_y_r,
                                sig_s_g,sig_s_p,sig_s_r,
                                sig_o_p,phi_r,
                                So,m0,m1,
                                nt,nk,nsmin,ngmin,
                                DDFUN,
                                Sg
                                ) )
        # simulate starting resident dynamics
    }
    ami <- es$am[i] + rnorm(1,0,smut_m)
    bmi <- es$bm[i] + rnorm(1,0,smut_m)
    # asi <- es$as[i] * exp(rnorm(1,0,smut_s))
    # bsi <- es$bs[i] * exp(rnorm(1,0,smut_s))
    # abri <- plogis( (qlogis((es$abr[i]+1)/2) + rnorm(1,0,smut_r)) )*2 - 1 
      # transform [0,1] to correlation range of [-1,1]
    # if(abri==-1) abri <- -0.99
    # if(abri==+1) abri <- +0.99
    invaded <- invade(zw[,2],ami,bmi,rd$Gres,rd$Ye,So,nt,nb) # asi,bsi,abri,
    if(i < nr){
      if(invaded==TRUE){
        es[i+1,] <- c(ami,bmi) # asi,bsi,abri
        rd <- ressim(zw[,2],x_z,ami,bmi,# asi,bsi,abri,
                     beta_p,beta_r,
                     eps_y_p,eps_y_r,
                     sig_s_g,sig_s_p,sig_s_r,
                     sig_o_p,phi_r,
                     So,m0,m1,
                     nt,nk,nsmin,ngmin,
                     DDFUN,
                     Sg,
                     ) # simulate new resident dynamics
      } 
      if(invaded==FALSE){
        es[i+1,] <- es[i,]
      } 
    }
  } # close i loop
  
  if(lastonly==T){
    return(es[nr,])
  }
  
  if(lastonly==F){
    outlist <- list(zw=zw,es=es)  
    return(outlist)
  }

}

multievolve <- function(  
  ni,nj,nr,nt,nb,nk,
  zam,wam,zsd,wsd,rho=0.82,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_s_g=NULL,sig_s_p=NULL,sig_s_r=NULL,
  sig_o_p,phi_r,
  m0,m1,
  am0,bm0,
  DDFUN=BHS,
  Sg=1,
  smut_m=0.5,# smut_s=0.1,smut_r=0.1,
  nsmin=10^-10,
  ngmin=10^-50,
  lastonly=T,
  savefile=NULL
  ){
  
  cur_date <- format(Sys.Date(),"%d%b%Y")
  
  amm <- bmm <- matrix(nr=ni,nc=nj)
  
  for(i in 1:ni){
    
    for(j in 1:nj){ 
      
    set.seed(i) 
      # climate and random years effects the same for all species
      # but parameter drift isn't (disrupted by differences in timing
      # of invasions, i.e. population dynamics simulations)
    
    ESS <- evolve(
      nr=nr,nt=nt,nb=nb,nk,
      zam=zam,wam=wam,zsd=zsd,wsd=wsd,rho=0.82,
      beta_p=beta_p[i,j,],beta_r=beta_r[i,j,],
      sig_y_p=sig_y_p[i,j],sig_y_r=sig_y_r[i,j],
      sig_s_g[i,j],sig_s_p[i],sig_s_r[i],
      sig_o_p=sig_o_p[i],phi_r=phi_r[i],
      m0=m0[i,j],m1=m1[i,j],
      am0=am0[i],bm0=bm0[i],
      DDFUN,
      Sg,
      smut_m,
      nsmin,
      ngmin,
      lastonly
    )
  
    amm[i,j] <- ESS$am
    bmm[i,j] <- ESS$bm
    
    } # j loop
  
  } # i loop
    
  outlist <- list(
    zam=zam,zsd=zsd,
    ni=ni,nj=nj,nr=nr,nt=nt,nt=nt,nb=nb,
    am=amm,bm=bmm
    )
  
  if(is.null(savefile)){
    return(outlist)
  }
  
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
    
}

simcombine <- function(insiml){
  
  varl <- insiml[[1]]
  nvar <- length(varl)
  dimvar <- lapply(varl,dim)
  ndimvar <- sapply(dimvar,length)
  firstmpos <- match(1:nclim,mpos)
  
  # Combine sim matrices by sim type
  outsiml <- vector("list", nvar)
  names(outsiml) <- names(varl)
  for(i in 1:nvar){
    if(ndimvar[i]>0){
      for(j in 1:nclim){
        for(k in 1:cpc){
          if(k==1) ast <- insiml[[firstmpos[j]]][[i]]
          else ast <- abind(ast,insiml[[firstmpos[j]+(k-1)]][[i]],along=1)
        }
        if(j==1){
          outsiml[[i]] <- array(ast,dim=c(dim(ast),1))
        }
        else{
          outsiml[[i]] <- abind(outsiml[[i]],ast,along=ndimvar[i]+1)
        }
      }
      if(nclim>1) dimnames(outsiml[[i]])[[ndimvar[i]+1]] <- cnames_unique
    }
  }
  outsiml <- outsiml[ndimvar>0]
  return(outsiml)
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