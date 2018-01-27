### Assemble individual RDS ESS outputs into a single RDS file

# Parsing arguments -------------------------------------------------------

args <- commandArgs(trailingOnly=TRUE)

storepath <- args[1]
label     <- args[2]

# Read in files -----------------------------------------------------------

pd <- readRDS(paste0("Sims/ESS_input_",label,".rds"))
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
  filename <- paste0(storepath,"/ESS_output_",i,"_",label,".rds")
  if(file.exists(filename)){
    curlist <- readRDS(filename)
    zw[,,pd$species[i],pd$iteration[i]] <- curlist$zw
    es[,,pd$species[i],pd$iteration[i]] <- as.matrix(curlist$es)
  }
}

# Save output -------------------------------------------------------------

saveRDS(zw,paste0("Sims/ESS_climate_",label,".rds"))
saveRDS(es,paste0("Sims/ESS_strategies_",label,".rds"))