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
smut_a=5,smut_b=5,# smut_s=0.1,smut_r=0.1,
nsmin=10^-10,
ngmin=10^-50,
lastonly=T,
savefile=NULL


# Inputs ------------------------------------------------------------------

# read-in file path (specified in submit script)
# ij <- task number

# Read in data ------------------------------------------------------------

ind <- readRDS("Sims/invasion_inputs.rds")
params <- ind$params[ij,]
series <- ind$series[,,ij]

# Calculate and save ESS --------------------------------------------------

source("Source/invasion_helper.R")
ess <- evolve(params=pd[NI,],series) # normalised climate + strategy evolution
saveRDS(ess,paste0("Sims/",savefile,"_",cur_date,".rds"))

  ){
  
  require(MASS)
  set.seed(NI) # keeps same year random effects structure for a given iteration
      
  zam=zamo+mam*zsdo,zsd=zsdo*msd,
  wam=wamo+mam*wsdo,wsd=wsdo*msd,
  zw_mu <- c(zam,wam) - log(tau_p)
  zw_sig <- matrix(c(zsd^2,rep(rho*zsd*wsd,2),wsd^2),nr=2,nc=2)
  zw <- mvrnorm(n=nt, mu=zw_mu, Sigma=zw_sig)
  z <- zw[,1]
  w <- zw[,2]
  
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
  
  dout <- cbind(zw=zw,es=es)
  
  
  if(lastonly==T){
    return(
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
  smut_a=5,smut_b=5,# smut_s=0.1,smut_r=0.1,
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
