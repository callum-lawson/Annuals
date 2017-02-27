### Correct names of species in Census and Seed data
### and remove low-abundance species

# maindir <- "D:/Users/calluml/Dropbox/NIOO/"
# maindir <- "C:/Users/Callum/Documents/My Dropbox/NIOO/"
setwd(paste0(maindir,"Analyses/Venable"))

####################
### READ IN DATA ###
####################

### SPECIES LIST 

sl <- read.csv(paste(maindir,
	"Data/Venable/species_list.csv",sep=""),
	header=T)

names(sl)[names(sl)=="Species.code"] <- "Code"
sl <- sl[,names(sl)!="synonym"]
slsp <- levels(sl$Code)

### CENSUS DATA 

cd <- read.csv(paste(maindir,
	"Data/Venable/census_data_07May2015.csv",sep=""),
	header=T)

cd$species <- as.factor(tolower(as.character(cd$species))) 
	# converts to lower case (some had species in caps)

cdsp <- levels(cd$species)

### SEED BANK DATA

sb <- read.csv(paste(maindir,
	"Data/Venable/Seed_Bank_13May2015.csv",sep=""),
	header=T)

sb$species <- as.factor(tolower(as.character(sb$species))) 
	# converts to lower case (later measurements had species in caps)

sbsp <- levels(sb$species)

######################
### CENSUS vs LIST ###
######################

### Present in LIST but missing in CENSUS:
# brto, cimo, euch, fiar, olli, stca

setdiff(slsp,cdsp)

### Present in CENSUS but missing in LIST:
# crsp	1587
# da-SP 	1 
# eriog 	68
# lo-as 	123
# losp 	51
# phsp 	101
# plsp 	4

( cdex <- setdiff(cdsp,slsp) )
table(cd$species)[cdsp %in% cdex]

####################
### SEED vs LIST ###
####################

### Present in LIST but missing in SEED:
# asnu, chri, crco, erab, fiar, gist, libi, losq, oepr, plar

setdiff(slsp,sbsp)

### Present in SEED but missing in LIST:
# crma 	1
# crsp 	4110
# eracil	1
# ersp 	4126
# eusp 	4126
# kagr 	1
# losp 	4126
# pesp 	1796
# phsp 	4122
# plsp 	1796

( sbex <- setdiff(sbsp,slsp) )
table(sb$species)[sbsp %in% sbex]

######################
### CENSUS vs SEED ###
######################

table(cd$species)[setdiff(cdsp,sbsp)]
	# SEED data present for all important CENSUS species
table(sb$species)[setdiff(sbsp,cdsp)]
	# some good SEED data for species missing in CENSUS

shared <- intersect(cdsp,sbsp)

############################
### THREE-WAY COMPARISON ###
############################

### Present in CENSUS and SEED but missing in LIST
# crsp, losp, phsp, plsp

setdiff(shared,slsp)
	# all "sp" - ID unclear, and multiple species of that genus in data
	# therefore removed:

allshared <- intersect(shared,slsp)

###########################
### REMOVE RARE SPECIES ###
###########################

nrecords <- table(cd$species)
minindivs <- 500
bigger <- names(nrecords[nrecords>minindivs])
bigger <- bigger[bigger!="brru"]
	# brru has many germinants but basically all zeroes for 
	# number of new (cd) and old seeds (sb)
allsharedbig <- intersect(allshared,bigger)

### HISTORGRAM
# pdf(paste0("species_abundances_",format(Sys.Date(),"%d%b%Y"),".pdf"),
# 	width=20,height=7)
# barplot(log10(sort(table(cd$species),decreasing=T)),ylab=expression(log[10](Nrecords)),las=2)
# abline(h=2,col="red",lty=2)
# dev.off()

############################
### CREATE NEW CSV FILES ###
############################

### CENSUS
cdnew <- subset(cd,species %in% allsharedbig)
write.csv(cdnew,file=paste0("census_data_sharedspecies_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)

### SEED BANK
sbnew <- subset(sb,species %in% allsharedbig)
write.csv(sbnew,file=paste0("Seed_Bank_sharedspecies_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)

### SPECIES LIST

slnew <- subset(sl,Code %in% allsharedbig)
write.csv(slnew,file=paste0("species_list_",format(Sys.Date(),"%d%b%Y"),".csv"),row.names=F)
