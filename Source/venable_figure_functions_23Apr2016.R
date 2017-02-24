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

	pdf(paste0(plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
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

	pdf(paste0(plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
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
	pdf(paste0(plotname,format(Sys.Date(),"%d%b%Y"),".pdf"),
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

seriesplot <- function(simname,a,yname,cols,ltys,colledgetext,detledgetext){

	pdf(paste0("timeseries_",simname,"_",format(Sys.Date(),"%d%b%Y"),".pdf"),
			width=plotwidth,height=plotheight)
	
		plotsetup()
		
		colvec <- rep(cols,each=nltys)
		ltyvec <- rep(ltys,times=ncols)

		for(j in 1:nspecies){
		
			matplot(t(a[,,j]),type="l",
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



