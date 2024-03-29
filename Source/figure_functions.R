########################################################
### Functions to plot figures from Stan model output ###
########################################################

### SMALL FUNCTIONS AND GRAPHICAL PARAMS

myjitter <- function(x) jitter(x,factor=10)
mysupsmu <- function(x,y) supsmu(x,y,bass=0)

plotsetup <- function(){
	par(mfrow=c(6,4),mar=c(2.5,2.5,4,2.5),oma=c(4,4,1,1),las=1,bty="l")
	}

lettlab <- function(i){
	mtext(bquote(bold((.(letters[i])))~.(spvals[i])),side=3,line=0.5,xpd=T,adj=0.05)
	}

addxlab <- function(xname){
	mtext(xname,side=1,line=4,las=1,xpd=T)
	}

addylab <- function(yname){
	mtext(yname,side=2,line=3.5,las=3,xpd=T)
	}

addledge <- function(ltext,...){
	plot(1,1,type="n",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
  	legend("center",legend=ltext,bty="n",...)
	}

plotwidth <- 10
plotheight <- 16

myblue <- rgb(red=0, green=0, blue=1, 0.25, maxColorValue = 1)
myred <- rgb(red=1, green=0, blue=0, 0.25, maxColorValue = 1)
black_rgb <- col2rgb("black")
myblack <- rgb(red=black_rgb[1],green=black_rgb[2],blue=black_rgb[3],alpha=50,maxColorValue = 255)

### PREDICTIONS AND DATA

preddatplot <- function(plotname,xdat,ydat,sind,xpred,ypred,xptype="vector",pcol,xname,yname,oneoneline=F,mycex=NULL,...){

	# sind = species index
	# pcol = point colour

	pdf(paste0("Plots/",plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
		width=plotwidth,height=plotheight)
	
	plotsetup()

	if(is.null(mycex)) mycex <- rep(1,length(ydat))

	for(i in 1:nspecies){
		pdat <- is.na(ydat)==F & sind==spvals[i] # data to plot
		y <- ydat[pdat]
		x <- xdat[pdat]
		plot(y~x,pch=21,col="black",bg=pcol,cex=mycex[pdat],...)
		# lines(mysupsmu(x,y),col="gray")
		if(xptype=="vector"){
			matplot(xpred,ypred[,,i],type="l",col="black",lty=c(2,1,2),add=T)
			}
		if(xptype=="matrix"){
			matplot(xpred[,i],ypred[,,i],type="l",col="black",lty=c(2,1,2),add=T)
			}
		if(oneoneline==T){
			abline(0,1,lty=3)
			}
	
		lettlab(i)

		if(i %in% 19:23) addxlab(xname) 
		if(i %in% seq(1,23,4)) addylab(yname) 

		}

	dev.off()

	}

### CALIBRATION

calplotall <- function(xdat,ydat,...){
	
	plot(ydat~xdat,
		pch=16,
		col=myblack,
		...
		)
	abline(0,1,col="red",lty=2)
	lines(mysupsmu(xdat,ydat),col="blue")

	}

	# add r-squared bit

calplot <- function(plotname,ym,pcol,xname,yname,...){
	
	# sind = species index
	# pcol = point colour
	# require(Hmisc)

	pdf("Plots/",paste0(plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
		width=plotwidth,height=plotheight)
	
	plotsetup()

	counter <- 0

	for(i in 1:nspecies){
		if( sum(is.na(ym[,3,i])==F & is.nan(ym[,3,i])==F)>=2 ){

			counter <- counter + 1

			# errbar(ym[,1,i],ym[,3,i],ym[,2,i],ym[,4,i],pch=16,col=pcol,xlab="",ylab="",...)
			plot(ym[,3,i],ym[,1,i],pch=16,col=pcol,xlab="",ylab="",...)
				# pred(x) vs obs(y)

			lettlab(i)

			if(i %in% 20:23) addxlab(xname) 
			if(counter %in% seq(1,23,4)) addylab(yname) 

			abline(0,1,col="red",lty=2)

			}
		}

	dev.off()

	}

### DATA HISTOGRAM

comphist <- function(plotname,x,y1,y2,xname,...){
	pdf(paste0("Plots/",plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
		width=plotwidth,height=plotheight)

	plotsetup()

	for(i in 1:nspecies){
		mybreaks <- hist(y1[x==i],
			col=myred,
		  density=20,
			main=spvals[i],
			prob=T,	
		  breaks=50,
			...
			)$breaks
		hist(y2[x==i],
			col=myblue,
			add=T,
			prob=T,
		  breaks=c(-Inf,mybreaks,Inf),
			...
			)

		if(i %in% 20:23) addxlab(xname) 
		if(i %in% seq(1,23,4)) addylab("Probability density") 
		}

	addledge(ltext=c("data","predictions"),fill=c(myred,myblue))

	dev.off()

	}

### CORRELATION HISTORGRAM

flexxlab <- function(xname,myline=5.25){
	mtext(xname,side=1,las=1,xpd=T,line=myline)
	}

corhist <- function(x,...){
	require(MASS)
	hist(x,prob=T,main="",breaks=30,xlim=c(-1,1),...)
	abline(v=0,col="red",lty=2)
	}

lettlab2 <- function(i,myline=-0.5,...){
	mtext(bquote(bold((.(letters[i])))),side=3,line=myline,xpd=T,adj=0.05,...)
	}

### SIMULATION OUTPUT FIGURES

seriesplot <- function(qa,varname,yname,quantiles=T){

  ncols <- length(cols)
  ltys <- c(3,1,3)
  nltys <- length(ltys)

  if(quantiles==T){
    ca <- acast(melt(qa),Var2 ~ Var4 + Var1 ~ Var3)
    colvec <- rep(cols,each=nltys)
    ltyvec <- rep(ltys,times=ncols)
  }
  if(quantiles==F) {
    ca <- acast(melt(qa),Var1 ~ Var3 ~ Var2)
    colvec <- cols
    ltyvec <- rep(1,times=ncols)
  }
    
	pdf(paste0("Plots/timeseries_",cnames_merged,"_",varname,"_",
	  format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
	
		plotsetup()
		
		for(j in 1:nspecies){
		
			matplot(ca[,,j],type="l",
				col=colvec,
				lty=ltyvec,
				xlab="",ylab=""
				)

			lettlab(j)

			if(j %in% 19:23) addxlab("Time (yr)") 
			if(j %in% seq(1,23,4)) addylab(yname) 

	  	}
	
		addledge(ltext=colledgetext,col=unique(colvec),lty=1)
		addledge(ltext=detledgetext)

	dev.off()
	
	}

relchange <- function(a,scenbase,scennew,tpos=15,keepsp){
  sposbase <- match(scenbase,cnames_unique)
  sposnew <- match(scennew,cnames_unique)
  nnew <- length(scennew)
  nkspecies <- sum(keepsp)
  anew <- matrix(nr=nnew,nc=nkspecies,
    dimnames=list(cnames_unique[sposnew],1:nkspecies)
  )
  for (i in 1:nnew){
    anew[i,] <- a[2,tpos,keepsp,sposnew[i]] - a[2,tpos,keepsp,sposbase]
  }
  return(t(anew))
  # 2 because median (ignoring other quantiles)
  # gives matrix; could add functions to convert to dataframe
}

pairplot <- function(plotname,a,npdim,w=8,h=8){
  pdf(paste0("Plots/",plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),width=w,height=h)
  # npdim = new page dim = whcih dim gets new page each time?
  npdim_n <- dim(a)[npdim]
  for(i in 1:npdim_n){
    if(npdim==2){
      pairs(a[,i,],lower.panel=panel.cor,pch="+",main=dimnames(a)[[npdim]][i])
    }
    if(npdim==3){
      pairs(a[,,i],lower.panel=panel.cor,pch="+",main=dimnames(a)[[npdim]][i])
    }
  }
  dev.off()
}

# withinplot <- function(parlist,simlist,parname,simname,
#   partrans_fun=NULL,simtrans_fun=NULL,
#   agg_fun=NULL,smooth=T,...){
#   
#   pdf(paste0("Plots/",parname,"_",simname,"_",format(Sys.Date(),"%d%b%Y"),".pdf"),
#     width=plotwidth,height=plotheight)
#   
#   iters <- as.vector(unlist(itersetl))
#     # same parameter sets for all clims
#   
#   for(i in 1:nclim){
#     
#     # may need to change iters extracted when change sim code
#     if(is.null(partrans_fun)) xmat <- parlist[[parname]][iters,]
#     else xmat <- partrans_fun( parlist[[parname]][iters,] )
#     
#     if(is.null(simtrans_fun)) yarr <- simlist[[i]][[simname]]
#     else yarr <- simtrans_fun(simlist[[i]][[simname]])
#     
#     if(is.null(agg_fun)) ymat <- yarr[,tpos,]
#     else ymat <- apply(yarr,c(1,3),agg_fun)
#     
#     plotsetup()
#     
#     for(j in 1:nspecies){
#       
#       plot(ymat[,j]~xmat[,j],...)
#       if(smooth==T) lines(mysupsmu(xmat[,j],ymat[,j]),col="red")
#       
#       lettlab(j)
#       
#       if(j %in% 19:23) addxlab(parname) 
#       if(j %in% seq(1,23,4)) addylab(simname) 
#       
#     }
#     
#     addledge(ltext=cnames_unique[i])
#     
#   }
#   
#   dev.off()
#   
# }

densplot <- function(na,yname,t=tpos,vab=NULL,...){
    
    pdf(paste0("Plots/density_",yname,"_",
      format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
    
    plotsetup()
    
    for(j in 1:nspecies){
      
      plot(density(na[,t,j,1],na.rm=T),col=cols[1],main="",
        xlim=round(range(na[,t,j,][is.finite(na[,t,j,])]),0),
        ...
        )
      for(i in 2:nclim){
        lines(density(na[,t,j,i],na.rm=T),col=cols[i])
        }
      if(!is.null(vab)) abline(v=vab[j],lty=3)
      
      lettlab(j)
      
      if(j %in% 19:23) addxlab(yname) 
      if(j %in% seq(1,23,4)) addylab("density") 
      
    }
    
    addledge(ltext=colledgetext,col=cols,lty=1)
    addledge(ltext=detledgetext)
    
    dev.off()
    
  }

parplot <- function(x,y,xname,yname,t=NULL,tran=25,...){
  
  cols_rgb <- col2rgb(cols)
  trancols <- rgb(
    red=cols_rgb[1,],
    green=cols_rgb[2,],
    blue=cols_rgb[3,],
    alpha=tran,
    maxColorValue = 255
  )
  
  xdim <- dim(x)
  ydim <- dim(y)
  nxdim <- length(dim(x))
  nydim <- length(dim(y))
  
  if(!is.null(t)){
    pdf(paste0("Plots/parplot_",xname,"_",yname,"_t",t,"_",format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
  }
  
  if(is.null(t)){
    pdf(paste0("Plots/parplot_",xname,"_",yname,"_",format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
  }

  plotsetup()
  
  for(j in 1:nspecies){
    
    if(!is.null(t)){
      ys <- y[,t,j,]
      if(nxdim==4){
        xs <- x[,t,j,]
      }
      if(nxdim==3){
        xs <- x[,t,]
      }
      if(nxdim==2){
        xs <- array(dim=dim(ys))
        xs[] <- rep(x[,j],nclim) # not necessary - matplot
      }
    }
    if(is.null(t)){
      if(nydim==4){
        ys <- acast(melt(y[,,j,]), Var1 + Var2 ~ Var3)
        if(nxdim==4){
          xs <- acast(melt(x[,,j,]), Var1 + Var2 ~ Var3)
        }
        if(nxdim==3){
          xs <- acast(melt(x), Var1 + Var2 ~ Var3)
        }
      }
      if(nydim==3){
        ys <- y[,j,]
        if(nxdim==3){
          xs <- x[,j,]
        }
        if(nxdim==2){
        xs <- array(dim=dim(ys))
        xs[] <- rep(x[,j],nclim) # not necessary - matplot
        }
      }
    }
    
    matplot(xs,ys,col=trancols,pch=16,...)
    for (k in 1:nclim) {
      lines(mysupsmu(xs[,k],ys[,k]),col=cols[k])
    }  
    
    lettlab(j)
    
    if(j %in% 19:22) addxlab(xname) 
    if(j %in% seq(1,23,4)) addylab(yname) 
    
  }
  
  addledge(ltext=colledgetext,col=cols,lty=1)
  addledge(ltext=detledgetext)
  
  dev.off()
  
}


# Parplot with density overlay --------------------------------------------

pardensplot <- function(x,y,xname,yname,t=NULL,tran=25,...){
  
  cols_rgb <- col2rgb(cols)
  trancols <- rgb(
    red=cols_rgb[1,],
    green=cols_rgb[2,],
    blue=cols_rgb[3,],
    alpha=tran,
    maxColorValue = 255
  )
  
  xdim <- dim(x)
  ydim <- dim(y)
  nxdim <- length(dim(x))
  nydim <- length(dim(y))
  
  if(!is.null(t)){
    pdf(paste0("Plots/pardensplot_",xname,"_",yname,"_t",t,"_",format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
  }
  
  if(is.null(t)){
    pdf(paste0("Plots/pardensplot_",xname,"_",yname,"_",format(Sys.Date(),"%d%b%Y"),".pdf"), width=plotwidth,height=plotheight)
  }
  
  plotsetup()

  dcur <- round(relchange(q_ns,scenbase="mu1_sd0",scennew="mu1_sd1",keepsp=rep(T,nj)),2)
  dfut <- round(relchange(q_ns,scenbase="mu1_sd1",scennew="mu095_sd12",keepsp=rep(T,nj)),2)
  # build flexibility in later if needed
  
  for(j in 1:nspecies){
    
    if(!is.null(t)){
      ys <- y[,t,j,]
      if(nxdim==4){
        xs <- x[,t,j,]
      }
      if(nxdim==3){
        xs <- x[,t,]
      }
      if(nxdim==2){
        xs <- array(dim=dim(ys))
        xs[] <- rep(x[,j],nclim) # not necessary - matplot
      }
    }
    if(is.null(t)){
      if(nydim==4){
        ys <- acast(melt(y[,,j,]), Var1 + Var2 ~ Var3)
        if(nxdim==4){
          xs <- acast(melt(x[,,j,]), Var1 + Var2 ~ Var3)
        }
        if(nxdim==3){
          xs <- acast(melt(x), Var1 + Var2 ~ Var3)
        }
      }
      if(nydim==3){
        ys <- y[,j,]
        if(nxdim==3){
          xs <- x[,j,]
        }
        if(nxdim==2){
          xs <- array(dim=dim(ys))
          xs[] <- rep(x[,j],nclim) # not necessary - matplot
        }
      }
    }
    
    propheight <- 0.25

    xd <- apply(xs,2,density,na.rm=T)
    xdymax <- sapply(xd,function(a) max(a$y))
    xdyscale <- propheight * diff(range(ys)) / max(xdymax)
    
    matplot(xs,ys,col=trancols,pch=16,...)
    for (k in 1:nclim) {
      lines(mysupsmu(xs[,k],ys[,k]),col=cols[k])
      xd[[k]]$y <- xd[[k]]$y * xdyscale - abs(min(ys))
      lines(xd[[k]],col=cols[k])
    }  

    lettlab(j)
    
    if(j %in% 19:22) addxlab(xname) 
    if(j %in% seq(1,23,4)) addylab(yname) 
    
    legend("topright",bty="n", 
      legend=c(
        paste("dcur",dcur[j]),
        paste(" dfut",dfut[j])
        )
      )
    
  }
  
  addledge(ltext=colledgetext,col=cols,lty=1)
  addledge(ltext=detledgetext)
  
  dev.off()
  
}

# Difference between climate scenarios ------------------------------------

diffplot <- function(xarr,yarr,refcl,t=NULL,xname,yname,xdiff=F,tran=50,...){
  
  if(!is.null(t)){
    plotname <- paste0("Plots/diffplot_",xname,"_",yname,"_",refcl,"_t",t,"_",
      format(Sys.Date(),"%d%b%Y"),".pdf")
  }
  else{
    plotname <- paste0("Plots/diffplot_",xname,"_",yname,"_",refcl,"_median_",
      format(Sys.Date(),"%d%b%Y"),".pdf")
  }
    
  pdf(plotname, width=plotwidth,height=plotheight)
  
  cols_rgb <- col2rgb(cols)
  trancols <- rgb(
    red=cols_rgb[1,],
    green=cols_rgb[2,],
    blue=cols_rgb[3,],
    alpha=tran,
    maxColorValue = 255
  )
  
  inclcl <- which(cnames_unique!=refcl) # drop reference climate
  for(m in inclcl){
    plotsetup()
    for(j in 1:nspecies){
      
      if(!is.null(t)){
        dyj <- yarr[,t,j,m] - yarr[,t,j,refcl]
        if(xdiff==F) xj <- xarr[,j]
        if(xdiff==T) xj <- xarr[,t,j,m] - xarr[,t,j,refcl]
      }
      else{
        dyj <- apply(yarr[,,j,m],1,median) - apply(yarr[,,j,refcl],1,median)
        if(xdiff==F) xj <- xarr[,j]
        if(xdiff==T) xj <- apply(xarr[,,j,m],1,median) - apply(xarr[,,j,refcl],1,median)
      }
      
      plot(xj,dyj,pch=16,col=trancols[m],...)
      lines(mysupsmu(xj,dyj),col="black")
      
      lettlab(j)
      
      if(j %in% 19:23) addxlab(xname) 
      if(j %in% seq(1,23,4)) addylab(yname) 
      
    }
    addledge(ltext=paste0(cnames_unique[m],"-",refcl))
  }
  
  dev.off()
}


# PIPs --------------------------------------------------------------------

pipplot <- function(z,xname,yname,pointvals=NULL,...){

  require(fields) # for tim.colors
  
  pdf(paste0("Plots/pipplot_",format(Sys.Date(),"%d%b%Y"),".pdf"), 
    width=plotwidth,height=plotheight)
  
  for(m in 1:nclim){
    
    plotsetup()
    
    for(j in 1:nspecies){
      
      image(z=pip[,,j,m],col=tim.colors(100),...)
      if(!is.null(pointvals)){
        points(pointvals[[1]][j],pointvals[[1]][j],pch="+",cex=1.5)
        points(pointvals[[2]][j],pointvals[[2]][j],pch=16,cex=1.5)
      }
      
      lettlab(j)
      
      if(j %in% 19:22) addxlab(xname) 
      if(j %in% seq(1,23,4)) addylab(yname) 
      
    }
    
  # addledge(ltext=colledgetext,col=cols,lty=1)
  # addledge(ltext=detledgetext)
  
  }
  
  dev.off()
  
}


### BEN BOLKER OVERDISPERSION

overdisp_fun <- function(model) {
  ## number of variance parameters in 
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
# but overdisp_fun apparently doens't work in Poisson not >5

### PANEL PLOTS

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
