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
