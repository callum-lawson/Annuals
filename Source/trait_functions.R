#############################################################################
# Functions for calculating derived trait values from popualtion parameters #
#############################################################################

# Reproduction ------------------------------------------------------------

zopt_f <- function(y,z){
  z[which(y==max(y))[1]]
}

minabs_f <- function(u,v){
  which(abs(u-v)==min(abs(u-v)))[1]
}

zwid_f <- function(y,z,qvals=c(0.25,0.75)){
  sumy <- cumsum(exp(y))
  qy <- sumy/max(sumy)
  qpos <- sapply(qvals,minabs_f,v=qy)
  zvals <- z[qpos]
  zwid <- diff(zvals)
}

zseq <- makeseq(log(seq(minprcp,maxprcp,length.out=2)/tau_p)) 
# was used in calculations but had not yet been saved as an object

rz <- pcr_prcp_lp$pred_lp
lz <- exp(rz)

zopt <- apply(rz,c(1,3),zopt_f,z=zseq)
zwid <- apply(rz,c(1,3),zwid_f,z=zseq)
rmax <- apply(rz,c(1,3),max,na.rm=T)
lsum <- apply(lz,c(1,3),sum)

rzrel <- aaply(exp(rz),c(1,3),function(x) x/sum(x))
rzrel <- log(aperm(rzrel,c(1,3,2)))

zopt_med <- apply(zopt,2,median,na.rm=T) # but check why getting NAs
zwid_med <- apply(zwid,2,median,na.rm=T)
rmax_med <- apply(rmax,2,median,na.rm=T)
lsum_med <- apply(lsum,2,median,na.rm=T)

rmed <- apply(rz,c(2,3),median)
rrelmed <- apply(rzrel,c(2,3),median)

# Mortality ---------------------------------------------------------------

Kncalc <- function(m0,m1,T3){
  exp(-m0*T3) / ( (m1/m0)*(1-exp(-m0*T3))/tau_s )
}
# Yodzis, fig on p53
# density effectively in units of 0.01m^2 (10cm x 10cm) 

hncalc <- function(m0,m1,T3){
  1 / ( (m1/m0)*(1-exp(-m0*T3))/tau_s )
}
# half-saturation density
# calculated from c1/(1+c2N) = (1/2)(c1/c2)
# based on Yodzis, p53

godmean_f <- function(alpha,beta){
  -alpha/beta
}

godvar_f <- function(beta){
  pi^2/(3*beta^2)
}
# from Godfray & Rees 2002
