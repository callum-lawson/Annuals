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
# In infinite-area models, all that matters for ES G is that same areas
# are used for plant and seed densities (scaling irrelevant - just alters 
# densities, not ES G)
 
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
  denom <- 1 - (phi/(mu+phi))^phi
  ifelse(denom==0 | mu==0, 0, mu/denom)
}
  # mean for hurdle model
  # ifelse prevents calculation failing when expected reproduction = 0

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
  # - can be done because just calculating mean of means for each plot type

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

pradj <- function(pr,mu,phi){
  q <- dnbinom(0, mu=mu, size=phi) # Pr(Y>0)
  return(pr / (1-q)) # zero-inflated
}

sprinkle <- function(x,kseq,probs){
  rmultinom(n=kseq, size=x, prob=probs)
}

ressim <- function(w,x_z,am,bm,# as,bs,abr,
                   beta_p,beta_r,
                   eps_y_p,eps_y_r,
                   sig_s_g,sig_s_p,sig_s_r,
                   sig_o_p,phi_r,theta_g,
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
        eta_bar_t <- sum(beta_r * x_t) + eps_y_r[t]
        pr_t <- logitnormint(mu=pi_bar_t,sigma=sig_o_p)
        rs_t <- nbtmean(exp(eta_bar_t),phi_r)
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
    
    Ye[t] <- ifelse(nn[t]==0, 0, nn[t] * DDFUN(nn[t]*(1-theta_g),m0,m1) / ng[t])
    if(t<nt) ns[t+1] <- ns[t] * ( (1-Gres[t])*So + Gres[t]*Ye[t] )
    t <- t + 1
    
  } # close t loop
  
  return(data.frame(Gres=Gres,Ye=Ye))
  # if(full==FALSE) return(data.frame(Gres=Gres,Ye=Ye))
  # if(full==TRUE) return(data.frame(Gres=Gres,Ye=Ye,ns=ns))
}

invade_infinite <- function(w,ami,bmi,Gres,Ye,So,nt,nb){ # asi,bsi,abri,
  
  if(NA %in% Ye){   
    invaded <- TRUE 
    # if resident goes extinct, invader establishes immediately
  }
  
  if(!NA %in% Ye){ # t = final value at which loop stopped 
    #Ginv <- coaG(w,ami,bmi,asi,bsi,abri)
    Ginv <- fixG(w,ami,bmi)
    delta_r <- log((1-Ginv)*So + Ginv*Ye) - log((1-Gres)*So + Gres*Ye)
    invaded <- mean(delta_r[(nb+1):nt]) > 0
  }

  return(invaded)
  
}

invade_finite <- function(w,x_z,am,bm,ami,bmi,
                          beta_p,beta_r,
                          eps_y_p,eps_y_r,
                          eps_s_g,eps_s_p,eps_s_r,
                          sig_o_p,phi_r,
                          So,m0,m1,
                          nt,nb,nk,nsmin,ngmin,
                          DDFUN,
                          Sg,
                          nc=5,   # n consecutive t that ns must be < nsmin
                          nstart=1,
                          intsd=10,
                          tau_d=100,
                          mumax=10^6){
  
  # require(truncdist)
  require(countreg)
  require(MASS)
  
  ns <- array(dim=c(nt,2)) # 2 = res and inv
  Gres <- fixG(w,am,bm)
  Ginv <- fixG(w,ami,bmi)

  ns[1,] <- c(round(nstart*nk/10,0),0)
  ns[1,][ns[1,]==0] <- 1
    # if starting density is 0, set to 1 instead
  t <- 1
  kseq <- 1:nk
  ng_k <- matrix(nr=nk,nc=2)
  
  while(t < nt
    & ns[t,1] > 0
    & (t <= (nb+1) | ns[t,2] > 0)
  ){
    
    if(t==(nb+1)){
      ns[t,2] <- 1
      # 1 invader introduced at t = nb + 1
    }
    
    ng <- rbinom(2,prob=c(Gres[t],Ginv[t]),size=ns[t,])
    no <- rbinom(2,prob=So,size=ns[t,]-ng)
    
    if(sum(ng)==0){
      nnb <- rep(0,2)
    }
    
    if(sum(ng)>0){
      ng_k[] <- sapply(ng,sprinkle,kseq,probs=eps_s_g)
      # using spatial terms as weights 
      # (normalised within function to sum to 1)

      ngt_k <- rowSums(ng_k)  # total germinants (residents and invaders)
      isg_k <- ngt_k>0        # binary: are there any germinants in plot?
      nkg <- sum(isg_k)       # total number of *plots* with germinants
      
      x_k <- array(dim=c(nkg,4))
      x_k[,1:3] <- rep(x_z[t,],each=nkg)
      x_k[,4] <- log(ngt_k[isg_k]) - log(tau_d/10)
      
      eps_o_p_k <- rnorm(nkg,0,sig_o_p)
      
      pr_k <- plogis(beta_p %*% t(x_k) + eps_y_p[t] + eps_s_p[isg_k] + eps_o_p_k)
      mu_k <- exp(beta_r %*% t(x_k) + eps_y_r[t] + eps_s_r[isg_k])
      mu_k[mu_k>mumax] <- mumax
        # maximum expected per-capita reproduction
        # prevents negative binomial distribution with p=0
      
      pradj_k <- pradj(pr_k,mu_k,phi_r)
      
      qu <- pradj_k <= 1 # quick plots
      nqu <- sum(qu)
      nsl <- nkg-nqu
      
      nr_kq <- array(dim=c(nqu,2))
      nr_ks <- array(dim=c(nkg-nqu,2))
      nr_kq[] <- rbinom(nqu*2,prob=rep(pradj_k[qu],2),size=ng_k[isg_k][qu])
      nr_ks[] <- rbinom(nsl*2,prob=rep(pr_k[!qu],2),size=ng_k[isg_k][!qu])
      
      nr_all <- sum(nr_kq) + sum(nr_ks)
      
      if(nr_all==0){
        nnb <- c(0,0)
      }

      if(nr_all>0){
        
        qp <- nr_kq > 0
        nn1m <- matrix(0,nr=nrow(nr_kq),nc=ncol(nr_kq))
        nn1m[qp] <- rnbinom(sum(qp), 
                  prob=phi_r/(phi_r+rep(mu_k[qu],2)[qp]), 
                  size=phi_r * nr_kq[qp]
                  )
        nn1 <- colSums(nn1m)
        
        whichpos <- which(nr_ks > 0)
        isinv_k <- factor(whichpos > nsl,levels=c(FALSE,TRUE)) 
        # if TRUE, then in 2nd column
        isinv_r <- rep(isinv_k,nr_ks[whichpos])
        mus <- rep(rep(mu_k[!qu],2)[whichpos],nr_ks[whichpos])
        y_all <- rztnbinom(n=sum(nr_ks), mu=mus, size=phi_r)
        # y_all <- rtrunc(n=nr_all,spec="nbinom",mu=mus,size=phi_r,a=0)
        nn2 <- tapply(y_all,isinv_r,sum) 
        nn2[is.na(nn2)] <- 0 # required when no reproducers in resident / invader
        nnt <- nn1 + nn2
        Sn <- ifelse(nnt==0,0,DDFUN(nnt/nk,m0,m1))
        # ifelse needed because doing separately for res and inv
        # division by 10 (i.e. scaling up to m^2) occurs within DDFUN
        nnb <- rbinom(2,prob=Sn,size=nnt)
        
      }
 
    }
    
    ns[t+1,] <- nnb + no
    t <- t + 1
    
  } # close t loop
    # simulation stops when either resident or invader has gone extinct
    # (or when maximum time limit has been reached -> coalition)
  
  if(ns[t,2] > 0){
    invaded <- TRUE
  }
  if(ns[t,2] == 0){
    invaded <- FALSE
  }
  
  return(invaded)
  
}

evolve <- function(
  nr,nt,nb,nk,
  zam,wam,zsd,wsd,rho,
  beta_p,beta_r,
  sig_y_p,sig_y_r,
  sig_s_g,sig_s_p,sig_s_r,
  sig_o_p,phi_r,theta_g,
  m0,m1,
  am0,bm0,
  # as0,bs0,
  # abr0,
  DDFUN,
  Sg,
  smut_a,smut_b,# smut_s=0.1,smut_r=0.1,
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
  
  es <- data.frame(am=rep(NA,times=nr),bm=rep(NA,times=nr))
    # as=rep(NA,times=nr), bs=rep(NA,times=nr),abr=rep(NA,times=nr)

  es[1,] <- c(am0,bm0) # ,as0,bs0,abr0)

  for(r in 1:nr){
    
    ami <- es$am[r] + rnorm(1,0,smut_a)
    bmi <- es$bm[r] + rnorm(1,0,smut_b)
    # asi <- es$as[r] * exp(rnorm(1,0,smut_s))
    # bsi <- es$bs[r] * exp(rnorm(1,0,smut_s))
    # abri <- plogis( (qlogis((es$abr[r]+1)/2) + rnorm(1,0,smut_r)) )*2 - 1 
    # transform [0,1] to correlation range of [-1,1]
    # if(abri==-1) abri <- -0.99
    # if(abri==+1) abri <- +0.99
    
    if(nk %in% c(0,Inf)){
      rd <- with(es[r,], ressim(zw[,2],x_z,am,bm,# as,bs,abr,
        beta_p,beta_r,
        eps_y_p,eps_y_r,
        sig_s_g,sig_s_p,sig_s_r,
        sig_o_p,phi_r,theta_g,
        So,m0,m1,
        nt,nk,nsmin,ngmin,
        DDFUN,
        Sg
      ) )
      invaded <- invade_infinite(zw[,2],ami,bmi,rd$Gres,rd$Ye,So,nt,nb) 
        # asi,bsi,abri,
    }
    
    if(nk>0 & nk<Inf){

      eps_s_p <- rnorm(nk,0,sig_s_p)
      eps_s_r <- rnorm(nk,0,sig_s_r)
        # done before g because want to match set.seed
      
      eps_s_g <- exp(rnorm(nk,0,sig_s_g))
      zsites <- rbinom(nk,size=1,prob=theta_g) 
      while(sum(zsites)==nk){
        zsites <- rbinom(nk,size=1,prob=theta_g) 
      }
        # theta = prob of zero
        # redraw until have at least one non-zero site
      eps_s_g[zsites==1] <- 0
      
      invaded <- with(es[r,], invade_finite(w=zw[,2],x_z,am,bm,ami,bmi,
                               beta_p,beta_r,
                               eps_y_p,eps_y_r,
                               eps_s_g,eps_s_p,eps_s_r,
                               sig_o_p,phi_r,
                               So,m0,m1,
                               nt,nb,nk,nsmin,ngmin,
                               DDFUN,
                               Sg
                               ))
    }
    
    if(r < nr){
      if(invaded==TRUE){
        es[r+1,] <- c(ami,bmi) # asi,bsi,abri
      } 
      if(invaded==FALSE){
        es[r+1,] <- es[r,]
      } 
    }
    
  } # close r loop
  
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
  sig_o_p,phi_r,theta_g=NULL,
  m0,m1,
  am0,bm0,
  DDFUN=BHS,
  Sg=1,
  smut_a=2,smut_b=2,# smut_s=0.1,smut_r=0.1,
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
    
      if(ni>1){
      
      ESS <- evolve(
        nr=nr,nt=nt,nb=nb,nk=nk,
        zam=zam,wam=wam,zsd=zsd,wsd=wsd,rho=0.82,
        beta_p=beta_p[i,j,],beta_r=beta_r[i,j,],
        sig_y_p=sig_y_p[i,j],sig_y_r=sig_y_r[i,j],
        sig_s_g=sig_s_g[i,j],sig_s_p=sig_s_p[i],sig_s_r=sig_s_r[i],
        sig_o_p=sig_o_p[i],phi_r=phi_r[i],
        theta_g=ifelse(is.null(theta_g),NULL,theta_g[i,j]),
        m0=m0[i,j],m1=m1[i,j],
        am0=am0[i],bm0=bm0[i],
        DDFUN,
        Sg,
        smut_a,smut_b,
        nsmin,
        ngmin,
        lastonly
      )
  
    }
      
    if(ni==1){
        
      ESS <- evolve(
        nr=nr,nt=nt,nb=nb,nk=nk,
        zam=zam,wam=wam,zsd=zsd,wsd=wsd,rho=0.82,
        beta_p=beta_p[j,],beta_r=beta_r[j,],
        sig_y_p=sig_y_p[j],sig_y_r=sig_y_r[j],
        sig_s_g=sig_s_g[j],sig_s_p=sig_s_p,sig_s_r=sig_s_r,
        sig_o_p=sig_o_p,phi_r=phi_r,
        theta_g=ifelse(is.null(theta_g),NULL,theta_g[j]),
        m0=m0[j],m1=m1[j],
        am0=am0,bm0=bm0,
        DDFUN,
        Sg,
        smut_a,smut_b,
        nsmin,
        ngmin,
        lastonly
      )
        
    }
      
    amm[i,j] <- ESS$am
    bmm[i,j] <- ESS$bm
    
    } # j loop
  
  } # i loop
    
  outlist <- list(
    zam=zam,zsd=zsd,
    ni=ni,nj=nj,nr=nr,nt=nt,nb=nb,
    am=amm,bm=bmm
    )
  
  if(is.null(savefile)){
    return(outlist)
  }
  
  if(!is.null(savefile)){
    saveRDS(outlist,paste0("Sims/",savefile,"_",cur_date,".rds"))
  }
    
}

simcombine <- function(insiml,nclim,cpc){
  
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