############################################################
# Population size and germination graphs for presentations #
############################################################

# Oslo talk plots ---------------------------------------------------------

tpos <- 15
myylim <- c(-5,0)

# MEANS AND VARIANCES

cclimpos <- which(cnames_unique=="mu1_cv1")
reld <- function(m){
  m[,2,tpos,] - rep(m[cclimpos,2,tpos,],each=dim(m)[1])
}
md_nsf <- reld(q_nsf)
# cclimpos = constant climate position (here mu1_cv1, later mu1_cv0)
# 2 = median

keepnames <- c("mu08_cv1","mu1_cv12","mu08_cv12")
# keepnames <- c("mu09_cv1","mu1_cv11","mu09_cv11")
keepclim <- cnames_unique %in% keepnames
keepsp <- (spvals %in% c("plpa","vuoc"))==F
md_nsf2 <- md_nsf[keepclim,keepsp]
md_nsd <- melt(md_nsf2)
names(md_nsd) <- c("scenario","species","dlN")
# dlN = difference relative to mu1_cv1

msy$olsdhat <- with(msy,totlive/area_o)
msy$godhat <- with(msy,germdhat+olsdhat)
G_med <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,median,na.rm=T)
)
md_nsd$G_med <- G_med[md_nsd$species]

medvarbet <- function(G,...){
  median(G*(1-G),...)
}
G_var <- with(subset(msy,year>2000),
  tapply(germdhat/godhat,species,medvarbet,na.rm=T)
)
md_nsd$G_var <- G_var[md_nsd$species]

myteal <- rgb(190,235,159,maxColorValue=255)

pdf(paste0("Output/mean_variance_comparison_",
  paste(keepnames,collapse="_"),
  "_t",tpos,"_",
  format(Sys.Date(),"%d%b%Y"),".pdf"
),
  width=4.5,height=4.5)

par(mfrow=c(1,1),mar=c(5,5,2,2),las=1,bty="l")
plot(dlN ~ scenario, 
  data=md_nsd,
  ylim=myylim,
  range=0,
  medlwd=2,
  boxwex=0.35,
  lty=1,
  ylab=expression(Delta~ln~population~size),
  xlab="Rainfall change scenario",
  names=c("mean","CV","both"),
  border="white" # makes invisible
)

abline(h=0,lty=2)

plot(dlN ~ scenario, 
  data=md_nsd,
  ylim=myylim,
  range=0,
  col=myteal,
  medlwd=2,
  boxwex=0.35,
  lty=1,
  ylab=expression(Delta~ln~population~size),
  xlab="Rainfall change scenario",
  names=c("mean","CV","both")
)

abline(h=0,lty=2)

# plot(dlN ~ qlogis(G_med),
#   data=subset(md_nsd,scenario==keepnames[3]),
#   ylim=myylim,
#   bg=myteal,
#   col="black",
#   pch=21,
#   cex=1.5,
#   ylab=expression(Delta~ln~population~size),
#   xlab="Species germination probability\n(logit scale)"
#   )
# 
# abline(h=0,lty=2)

dev.off()

# BES talk plots ----------------------------------------------------------

cclimpos <- which(cnames_unique=="mu1_cv0")
# changing reference climate to mu1_cv0
md_nsf3 <- reld(q_nsf)
md_nsd3 <- melt(md_nsf3[,keepsp])
names(md_nsd3) <- c("scenario","species","dlN")
G_mod <- plogis(q_Gf["mu1_cv1",2,tpos,])[keepsp]
consdiff <- md_nsd3$dlN[md_nsd3$scenario=="mu1_cv1"]

plot(dlN ~ qlogis(G_med[keepsp]), data=subset(md_nsd3,scenario=="mu1_cv1"))

iota_mu <- -go$alpha_G/go$beta_Gz
iota_mu_med <- apply(iota_mu,2,median)
iota_sig <- pi^2/(3*go$beta_Gz^2)
iota_sig_med <- apply(beta_G_var,2,median)
# from Godfray & Rees 2002

pairs(cbind(G_mod,mu=iota_mu_med[keepsp],lsig=log(iota_sig_med[keepsp])))
# lower G mostly reflects higher threshold
# (i.e. not necessarily bet-hedging - take longer to germinate but do so at
# same time)

plot(G_mod,consdiff)
plot(iota_mu_med[keepsp],consdiff)
plot(log(iota_sig_med[keepsp]),consdiff)
# nothing predicts whether a bet-hedger (cur env definition) 

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])
# less var -> seem to do slightly better under mean change
# note change from md_nsd3 to md_nsd

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu1_cv12"])
# less var -> do better under variance change

plot(G_mod[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
plot(iota_mu_med[keepsp],md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
plot(log(iota_sig_med[keepsp]),md_nsd$dlN[md_nsd$scenario=="mu08_cv12"])
# less var -> do slightly better under overall change

matplot(qlogis(G_med),t(q_Gf[,2,tpos,]),pch=16)

diffplot <- function(x,y,xname,...){
  plot(y ~ x,bg=myteal,col="black",pch=21,ann=F,cex=1.5,...)
  addxlab(xname)
}
# custlettlab <- function(i,label){
#   mtext(bquote(bold((.(letters[i])))~.(label)),
#     side=3,line=0.5,xpd=T,adj=0.05)
#   }
# addtoplab <- function(tname){
#   mtext(bquote(bold(tname)),side=3,line=2,las=1,xpd=T,cex=1.1)
#   }

diffplot_all <- function(...){
  xlabcur <- expression(Delta~ln~N~current~variability)
  par(mfrow=c(2,3),mar=c(5,2.5,5,2.5),oma=c(4,4,2,2),las=1,bty="l")
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu08_cv1"],xlabcur,...)
  mtext(bquote(bold(mean)),side=3,line=3,las=1,xpd=T,cex=1.3)
  addylab(expression(Delta~ln~population~size~(N)))
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu1_cv12"],xlabcur,...)
  mtext(bquote(bold(cv)),side=3,line=3,las=1,xpd=T,cex=1.3)
  diffplot(consdiff,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],xlabcur,...)
  mtext(bquote(bold(both)),side=3,line=3,las=1,xpd=T,cex=1.3)
  
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu08_cv1"],"germination probability",...)
  addylab(expression(Delta~ln~population~size~(N)))
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu1_cv12"],"germination probability",...)
  diffplot(G_mod,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],"germination probability",...)
}

pdf(paste0("Plots/bet_hedging_prediction_",
  paste(keepnames,collapse="_"),
  "_t",tpos,"_",
  format(Sys.Date(),"%d%b%Y"),".pdf"
),
  width=10,height=8
)
diffplot_all(type="n")
diffplot_all()
dev.off()

pdf(paste0("Plots/bet_hedging_prediction_K_So",
  "_t",tpos,"_",
  format(Sys.Date(),"%d%b%Y"),".pdf"
),width=9,height=4.5)
par(mfrow=c(1,2),mar=c(2,2,2,2),oma=c(4,4,1,1),las=1,bty="l")

diffplot(mta$lKnmed,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="ln max seed density",type="n")
addylab(expression(Delta~ln~population~size))
diffplot(mta$So,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="logit dormant seed survival",type="n")

diffplot(mta$lKnmed,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="ln max seed density")
addylab(expression(Delta~ln~population~size))
diffplot(mta$So,md_nsd$dlN[md_nsd$scenario=="mu08_cv12"],
  xname="logit dormant seed survival")

dev.off()

summary(glm(md_nsd$dlN[md_nsd$scenario=="mu08_cv12"]~consdiff))
summary(glm(md_nsd$dlN[md_nsd$scenario=="mu08_cv12"]~G_mod))

plot(md_nsd$dlN[md_nsd$scenario=="mu1_cv12"] ~ md_nsd$dlN[md_nsd$scenario=="mu08_cv1"])