### Assemble individual RDS ESS outputs into a single RDS file

# Parsing arguments -------------------------------------------------------

library(optparse)

parser <- OptionParser(
  usage       = "Rscript %prog datapath workpath label",
  description = "\ncombine ESS outputs into single RDS file"
)

cli <- parse_args(parser, positional_arguments = 3)
# other args: sourcepath, workpath

savepath <- paste0(cli$args[1],"Sims/")
workpath <- cli$args[2]
label    <- cli$args[3]

# Read in files -----------------------------------------------------------

pd <- readRDS(paste0(workpath,"ESS_input_",label,".rds"))
ntasks  <- nrow(pd) # = nj * ni * nm
finite <- with(pd[1,], nk>0 & nk<Inf)

zw <- with(pd[1,], array(dim=c(nt,2,nj,ni)))
if(finite==FALSE){
  es <- with(pd[1,], array(dim=c(nr,2,nj,ni)))
} 
if(finite==TRUE){
  es <- with(pd[1,], array(dim=c(nr,3,nj,ni)))
} 
  # extra column for extinction in finite populations

for(i in 1:ntasks){
  filename <- paste0(workpath,"ESS_output_",i,"_",label,".rds")
  if(file.exists(filename)){
    curlist <- readRDS(filename)
    zw[,,pd$species[i],pd$iteration[i]] <- curlist$zw
    es[,,pd$species[i],pd$iteration[i]] <- as.matrix(curlist$es)
  }
}

# Save output -------------------------------------------------------------

saveRDS(zw,paste0(savepath,"ESS_climate_",label,".rds"))
saveRDS(es,paste0(savepath,"ESS_strategies_",label,".rds"))