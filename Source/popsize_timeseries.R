######################################################################
### Plotting population sizes of seedlings and seeds through time  ###
######################################################################

### LOAD DATA 

msy <- read.csv("Output/msy_15Jan2016.csv",header=T)

#####################
### DEFINE PARAMS ###
#####################

### DATA PARAMS
nspecies <- nlevels(msy$species) # same for all other datasets
spvals <- levels(msy$species)
nyear <- length(unique(msy$year))

###################
### PLOT GRAPHS ###
###################

pdf(paste0("Plots/pop_timeseries_",format(Sys.Date(),"%d%b%Y"),".pdf"),
	width=20,height=16)
	
par(mfrow=c(6,4),mar=c(2.5,2.5,4,2.5),oma=c(4,4,1,1),las=1,bty="l")

for(i in 1:nspecies){
	newdat <- subset(msy,species==spvals[i], select=c("year","germd","nsd","olsd"))
	matplot(newdat$year, log10(newdat[,-1]), main=spvals[i], type="b", pch=16,
		col=c("darkgreen","red","black"),lty=1)
	}
plot(1,1,type="n",ann=F,xaxt="n",yaxt="n",bty="n")
legend("center",legend=c("germinants","new seed","old live seed"),
	col=c("darkgreen","red","black"),lty=1,bty="n")

dev.off()


