#############################################################
### Functions for calculating pr, rs, and pcr predictions ###
#############################################################

### LOAD PARAMS

go <- readRDS("go_pars_binomialG_tdistbpar_normndat_16Feb2017.rds")
gs <- readRDS("gnzhh_onhh_pars_medians_26Oct2015.rds")
pr <- readRDS("pr_pars_yearhet_squared_pc_02Mar2016.rds")
rs <- readRDS("rs_pars_yearhet_squared_pc_trunc_05Mar2016.rds")

### GO

# Gcalc <- function(z,alpha,beta){
#   plogis(alpha + beta*z)
#   }

Gcalc <- function(z,d,alpha,beta_z){
	plogis(alpha + beta_z*z) # d = log(dens)
  }

Scalc <- function(olsd,nsd,m0,m1,T2=dl$T2,T3=dl$T3){
	nd <- nsd * exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*nsd/tau_s )
		# tau_s only in second part
	  # (because only have to adjust survival rate, not density going in)
	od <- olsd * exp(-m0*(T2+T3))
	cd <- nd + od
	return(cd/(nsd+olsd))
  }

Sncalc <- function(nd,m0,m1,T3){
  exp(-m0*T3) / ( 1 + (m1/m0)*(1-exp(-m0*T3))*nd/tau_s )
  }

### RS 

maxiter <- max(length(pr$sig_o_p),length(rs$phi))

makeseq <- function(x){
	seq(min(x,na.rm=T),max(x,na.rm=T),length.out=nseq)
	}

makeav <- function(x,FUN){
	rep(FUN(x,na.rm=T),times=nseq)
	}

nbtmean <- function(mu,phi){
 	mu / ( 1 - (phi/(mu+phi))^phi )
	}

# below function can't be used for calibration because missing eps terms

r_lp <- function(betas,ytype,xtype,xarr=NULL,avmeth=mean,maxiter=niter){

		# mean of prcp/germd weighted by number of obs 
		# beta: [niter,nspecies,3]

	if(xtype!="data"){
		xarr <- array(dim=c(nseq,4,nspecies))
		pred_lp <- pred_lp_notrunc <- array(dim=c(maxiter,nseq,nspecies))
		}

	if(xtype=="data"){
		xarr <- xarr
		pred_lp <- pred_lp_notrunc <- array(dim=c(maxiter,dim(xarr)[1],nspecies)) 
			# changed "nyear" to "nseq"
		}

	for(i in 1:nspecies){

		# if(ytype=="pr") subdat <- subset(prdat,species==i)
		# if(ytype=="rs") subdat <- subset(rsdat,species==i) 
	  subdat <- subset(prdat,species==i)
    	# using prdat for both pr and rs predictions because want
	    # prcp and dens to have same range (e.g. for calculating pcr)
			# species is numeric in curdat

		if(xtype=="prcp"){
			xarr[,,i] <- c(
				rep(1,nseq),
				makeseq(log(subdat$prcp/tau_p)), # change if change variable
				makeseq(log(subdat$prcp/tau_p))^2, # change if change variable
				makeav(log(subdat$germd/tau_d),FUN=avmeth) # change if change variable
				)
			}			

		if(xtype=="dens"){
			xarr[,,i] <- c(
				rep(1,nseq),
				makeav(log(subdat$prcp/tau_p),FUN=avmeth), # change if change variable
				makeav(log(subdat$prcp/tau_p),FUN=avmeth)^2, # change if change variable
				makeseq(log(subdat$germd/tau_d)) # change if change variable
				)
			}

		if(ytype=="pr"){
			pred_lp[,,i] <- betas[rep_len(1:niter,maxiter),i,] %*% t(xarr[,,i])
			}
			# cols = seqno; rows = iter
			# rep_len -> resamples sequence until reaches a length of maxiter
			# so drawing from sampled parameters a set number (maxiter) of times

		if(ytype=="rs"){
			pred_lp_notrunc[,,i] <- betas[rep_len(1:niter,maxiter),i,] %*% t(xarr[,,i])
				# no squared: sqpos <- 3; xarr[,-sqpos,i]
			pred_lp[,,i] <- log(nbtmean(exp(pred_lp_notrunc[,,i]),rs$phi[i]))
			}
			# accounting for zero-truncation
			# don't need to account for variance because 
			# none within exp() function (site/year already accounted for)

		} # close i loop

	pred_lp[is.finite(pred_lp)==F] <- NA

	list(xarr=xarr,pred_lp=pred_lp)

	}

r_quant <- function(pred){

	nseq <- dim(pred$xarr)[1]	
	quant <- array(dim=c(nseq,nquant,nspecies))
	for(i in 1:nspecies){
		quant[,,i] <- aaply(pred$pred_lp[,,i],.margins=2,.fun=quantile,probs=myprobs,na.rm=T)
		}

	return(quant)
  	
  }

### OTHER FUNCTIONS

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x,y)
  r_abs <- abs(r)
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r_abs)
  }
