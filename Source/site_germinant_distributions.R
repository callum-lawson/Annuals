###############################################################################
# How does distributing germinants among discrete sites affect the among-site #
# variance in the number of germinants?                                       #
###############################################################################

obskvar <- function(nk,ng,ssd){
  kseq <- 1:nk
  epss <- rnorm(nk,0,ssd)
  obss <- as.numeric(table(
    factor(
      sample(kseq, 
        ng, 
        replace=T, 
        prob=exp(epss)
      ),
      levels=kseq
    )
  ))

  return(var(obss))
}
  
obskvar(10,10,0.1)

ssd <- 0.1
ng <- 10^6
nkseq <- 10^(1:6)
nnk <- length(nkseq)
nsim <- 100

ok <- tk <- matrix(nr=nsim,nc=nnk)
for(i in 1:nnk){
  for(j in 1:nsim){
    ok[j,i] <- obskvar(ng,nkseq[i],ssd)
  }
}

tk[] <- rep( (exp(ssd^2)-1)*exp(2*(log(ng/nkseq) - (ssd^2 / 2)) + ssd^2), each=nsim)
ok/tk






